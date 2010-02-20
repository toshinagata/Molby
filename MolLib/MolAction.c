/*
 *  MolAction.c
 *
 *  Created by Toshi Nagata on 07/06/23.
 *  Copyright 2007-2008 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "Ruby_bind/Molby.h"  /*  Support Ruby binding  */

#include "IntGroup.h"
#include "Molecule.h"
#include "MolAction.h"
#include "MD/MDCore.h"

#include <stdarg.h>
#include <ctype.h>

const char *gMolActionNone            = "none";
const char *gMolActionAddAnAtom       = "addAtom:Ai;i";
const char *gMolActionDeleteAnAtom    = "deleteAtom:i";
const char *gMolActionMergeMolecule   = "mergeMol:MG";
const char *gMolActionUnmergeMolecule = "unmergeMol:G";
const char *gMolActionAddBonds        = "addBonds:I";
const char *gMolActionDeleteBonds     = "deleteBonds:I";
const char *gMolActionAddAngles       = "addAngles:IG";
const char *gMolActionDeleteAngles    = "deleteAngles:G";
const char *gMolActionAddDihedrals    = "addDihedrals:IG";
const char *gMolActionDeleteDihedrals = "deleteDihedrals:G";
const char *gMolActionAddImpropers    = "addImpropers:IG";
const char *gMolActionDeleteImpropers = "deleteImpropers:G";
/* const char *gMolActionReplaceTables   = "replaceTables:III"; */
const char *gMolActionTranslateAtoms  = "translateAtoms:vG";
const char *gMolActionRotateAtoms     = "rotate:vdvG";
const char *gMolActionTransformAtoms  = "transform:tG";
const char *gMolActionSetAtomPositions = "atomPosition:GV";
const char *gMolActionInsertFrames    = "insertFrames:GV";
const char *gMolActionRemoveFrames    = "removeFrames:G";
const char *gMolActionSetSelection    = "selection:G";
const char *gMolActionChangeResidueNumber = "changeResSeq:Gi";
const char *gMolActionChangeResidueNumberForUndo = "changeResSeqForUndo:GIi";
const char *gMolActionChangeResidueNames = "changeResNames:IC";
const char *gMolActionOffsetResidueNumbers = "offsetResSeq:Gii";
const char *gMolActionChangeNumberOfResidues = "changeNres:i";
const char *gMolActionReorderAtoms    = "reorderAtoms:I";
const char *gMolActionExpandBySymmetry = "expandSym:Giiii";
const char *gMolActionDeleteSymmetryOperation = "deleteSymop";
const char *gMolActionAddSymmetryOperation = "addSymop:t";
const char *gMolActionSetCell         = "setCell:Di";
const char *gMolActionSetBox          = "setBox:vvvvi";
const char *gMolActionClearBox        = "clearBox";
/*const char *gMolActionSetParameterAttribute = "setParameterAttr:iirri"; */
const char *gMolActionAddParameters   = "addParameters:iGU";
const char *gMolActionDeleteParameters = "deleteParameters:iG";
const char *gMolActionCartesianToXtal  = "cartesianToXtal";
const char *gMolActionXtalToCartesian  = "xtalToCartesian";

/*  A Ruby array to retain objects used in MolActionArg  */
static VALUE sMolActionArgValues = Qfalse;

/*  Action arguments  */
/*  (Simple types)  i: Int, d: double, s: string, v: Vector, t: Transform, u: UnionPar
 (Array types)   I: array of Int, D: array of double, V: array of Vector, C: array of char, T: array of Transform, U: array of UnionPars
 (Complex types) M: Molecule, G: IntGroup, A: Atom
 (Ruby value)    b: Ruby boolean, r: Ruby object, R: an array of Ruby object (a Ruby array)
 (Return value from Ruby script) i, d, s, v, t + 0x80 */
typedef struct MolActionArg {
	Byte type;
	union {
		Int ival;
		double dval;
		struct {  /*  Array type: the size of the element is defined by the argument type  */
			Int nitems; /* Number of items */
			void *ptr;
		} arval;
		struct IntGroup *igval; /* The record is retained but not duplicated */
		struct Molecule *mval; /* The record is retained but not duplicated */
		struct Atom *aval; /* The value is duplicated, so any pointer can be passed */
		VALUE vval;        /* The value is retained in sMolActionArgValues  */
		void *pval;        /* Store address for the return value; only meaningful in performing Ruby script  */
	} u;
} MolActionArg;

/*  For debug  */
/* static int sMolActionCount = 0; */

/*  Create a new MolAction with the given arguments.
 *  Note: If the argument is an array type, it is represented by an integer (number
 *  of items) followed by the pointer.  */
MolAction *
MolActionNewArgv(const char *name, va_list ap)
{
	MolAction *action;
	const char *p;

	action = (MolAction *)malloc(sizeof(MolAction));
	if (action == NULL)
		goto low_memory;
	memset(action, 0, sizeof(MolAction));
	action->refCount = 1;
	action->frame = -1;
	action->name = strdup(name);

	/*  For debug  */
/*	fprintf(stderr, "MolAction %p created, count = %d\n", action, ++sMolActionCount); */
	
	/*  Handle arguments  */
	p = strchr(name, ':');
	if (p == NULL)
		return action;
	p++;

	while (*p) {
		MolActionArg arg;
		memset(&arg, 0, sizeof(MolActionArg));
		arg.type = *p++;
		switch (arg.type) {
			case 'i':
				arg.u.ival = va_arg(ap, Int);
				break;
			case 'd':
				arg.u.dval = va_arg(ap, double);
				break;
			case 'I':
			case 'D':
			case 'C':
			case 's':
			case 'V':
			case 'v':
			case 'T':
			case 't':
			case 'U':
			case 'u': {
				/*  Array type of some kind  */
				void *ptr;
				Int nitems, itemsize, allocsize;
				if (arg.type != 's' && arg.type != 'v' && arg.type != 't' && arg.type != 'u')
					nitems = va_arg(ap, Int);
				ptr = va_arg(ap, void *);
				switch (arg.type) {
					case 'I': itemsize = sizeof(Int); break;
					case 'D': itemsize = sizeof(double); break;
					case 'C':
					case 's': itemsize = sizeof(char); break;
					case 'V':
					case 'v': itemsize = sizeof(Vector); break;
					case 'T':
					case 't': itemsize = sizeof(Transform); break;
					case 'U':
					case 'u': itemsize = sizeof(UnionPar); break;
				}
				if (arg.type == 's')
					nitems = (ptr == NULL ? 0 : strlen((char *)ptr));
				else if (arg.type == 'v' || arg.type == 't' || arg.type == 'u')
					nitems = 1;
				arg.u.arval.nitems = nitems;
				allocsize = itemsize * nitems + (arg.type == 's'); /* +1 for string type */
				if (allocsize == 0)
					arg.u.arval.ptr == NULL;
				else {
					arg.u.arval.ptr = calloc(1, allocsize);
					if (arg.u.arval.ptr == NULL)
						goto low_memory;
					memmove(arg.u.arval.ptr, ptr, itemsize * nitems);
				}
				break;
			}
/*			case 's':
				arg.u.sval = va_arg(ap, char *);
				arg.u.sval = strdup(arg.u.sval);
				break; */
/*			case 'v': {
				Vector *vp = va_arg(ap, Vector *);
				arg.u.vval = (Vector *)malloc(sizeof(Vector));
				if (arg.u.vval == NULL)
					goto low_memory;
				*(arg.u.vval) = *vp;
				break;
			} */
			case 'G':
				arg.u.igval = va_arg(ap, IntGroup *);
				if (arg.u.igval != NULL)
					IntGroupRetain(arg.u.igval);
				break;
			case 'M':
				arg.u.mval = va_arg(ap, Molecule *);
				MoleculeRetain(arg.u.mval);
				break;
			case 'A': {
				Atom *aval = va_arg(ap, Atom *);
				arg.u.aval = (Atom *)malloc(gSizeOfAtomRecord);
				if (arg.u.aval == NULL)
					goto low_memory;
				AtomDuplicate(arg.u.aval, aval);
				break;
			}
			case 'r':
			case 'R': {
				int rtype;
				arg.u.vval = va_arg(ap, VALUE);
				rtype = TYPE(arg.u.vval);
				if (rtype != T_NIL && rtype != T_FIXNUM && rtype != T_SYMBOL && rtype != T_TRUE && rtype != T_FALSE) {
					if (sMolActionArgValues == Qfalse) {
						sMolActionArgValues = rb_ary_new();
						rb_global_variable(&sMolActionArgValues);
					}
					rb_ary_push(sMolActionArgValues, arg.u.vval);
				}
				break;
			}
			case 'b':
				arg.u.vval = (va_arg(ap, Int) ? Qtrue : Qfalse);
				break;
			case ';':
				if (*p == 'i' || *p == 'd' || *p == 's' || *p == 'v' || *p == 't' || *p == 'r') {
					arg.type = (*p | 0x80);
					p++;
					arg.u.pval = va_arg(ap, void *);
				} else {
					fprintf(stderr, "Internal error: unknown return-type \'%c\' in NewMolAction()", *p);
				}
				break;
			default:
				fprintf(stderr, "Internal error: unknown argument type \'%c\' in NewMolAction()", arg.type);
				break;
		}
		if (AssignArray(&action->args, &action->nargs, sizeof(arg), action->nargs, &arg) == NULL)
			goto low_memory;
	}
	return action;
				
  low_memory:
	if (action != NULL)
		free(action);
	Panic("Low memory during creation of action record");
	return NULL;  /* Not reached */
}

MolAction *
MolActionNew(const char *name, ...)
{
	va_list ap;
	MolAction *retval;
	va_start(ap, name);
	retval = MolActionNewArgv(name, ap);
	va_end(ap);
	return retval;
}

MolAction *
MolActionRetain(MolAction *action)
{
	if (action != NULL)
		action->refCount++;
	return action;
}

void
MolActionRelease(MolAction *action)
{
	int i;
	if (action == NULL)
		return;
	if (--action->refCount > 0)
		return;
	if (action->name != NULL)
		free((void *)action->name);
	for (i = 0; i < action->nargs; i++) {
		MolActionArg *argp = action->args + i;
		switch (argp->type) {
			case 'i':
			case 'd':
			case 'b':
				break;
			case 'I':
			case 'D':
			case 'C':
			case 's':
			case 'V':
			case 'v':
			case 'T':
			case 't':
			case 'U':
			case 'u':
				if (argp->u.arval.ptr != NULL)
					free(argp->u.arval.ptr);
				break;
			
		/*	case 'a':
				if (argp->u.iaval != NULL)
					free(argp->u.iaval);
				break;
			case 's':
				if (argp->u.sval != NULL)
					free(argp->u.sval);
				break;
			case 'v':
				if (argp->u.vval != NULL)
					free(argp->u.vval);
				break; */
			case 'G':
				if (argp->u.igval != NULL)
					IntGroupRelease(argp->u.igval);
				break;
			case 'M':
				MoleculeRelease(argp->u.mval);
				break;
			case 'A':
				if (argp->u.aval != NULL) {
					AtomClean(argp->u.aval);
					free(argp->u.aval);
				}
				break;
			case 'r':
			case 'R': {
				int rtype = TYPE(argp->u.vval);
				if (rtype != T_NIL && rtype != T_FIXNUM && rtype != T_SYMBOL && rtype != T_TRUE && rtype != T_FALSE) {
					/*  Look for the same value in sMolActionArgValues and remove it  */
					if (sMolActionArgValues != Qfalse) {
						VALUE *ptr = RARRAY_PTR(sMolActionArgValues);
						int n = RARRAY_LEN(sMolActionArgValues);
						while (--n >= 0) {
							if (ptr[n] == argp->u.vval) {
								rb_ary_delete_at(sMolActionArgValues, n);
								break;
							}
						}
					}
				}
				break;
			}
			default:
				break;
		}
	}
	if (action->args != NULL)
		free(action->args);
	free(action);

	/*  For debug  */
/*	fprintf(stderr, "MolAction %p disposed, count = %d\n", action, --sMolActionCount);	*/
}

void
MolActionSetFrame(MolAction *action, int frame)
{
	if (action != NULL)
		action->frame = frame;
}

typedef struct MolRubyActionInfo {
	Molecule *mol;
	MolAction *action;
	VALUE receiver;
	ID method_id;
	VALUE ary;
	char enable_interrupt;
	char return_type;
	void *return_ptr;
} MolRubyActionInfo;

static VALUE
s_MolActionToRubyArguments(VALUE vinfo)
{
	MolRubyActionInfo *info = (MolRubyActionInfo *)vinfo;
	VALUE val;
	int i, argc;
	char *s, *p;

	argc = info->action->nargs - 1;
	info->return_type = 0;
	info->return_ptr = NULL;
	info->ary = rb_ary_new2(argc);
	for (i = 1; i <= argc; i++) {
		MolActionArg *argp = &(info->action->args[i]);
		if (argp->type & 0x80) {
			/*  Return value specified  */
			info->return_type = (argp->type & 0x7f);
			info->return_ptr = argp->u.pval;
			break;  /*  This should be the last argument  */
		}
		switch (argp->type) {
			case 'i':
				val = INT2NUM(argp->u.ival);
				break;
			case 'd':
				val = rb_float_new(argp->u.dval);
				break;
			case 's':
			case 'C':
				val = rb_str_new((char *)argp->u.arval.ptr, argp->u.arval.nitems);
				break;
			case 'v':
				val = ValueFromVector((Vector *)argp->u.arval.ptr);
				break;
			case 't':
				val = ValueFromTransform((Transform *)argp->u.arval.ptr);
				break;
			case 'I':
			case 'D':
			case 'V':
			case 'T': {
				int j;
				val = rb_ary_new();
				for (j = 0; j < argp->u.arval.nitems; j++) {
					VALUE val1;
					void *ptr = argp->u.arval.ptr;
					switch (argp->type) {
						case 'I': val1 = INT2NUM(((Int *)ptr)[j]); break;
						case 'D': val1 = rb_float_new(((double *)ptr)[j]); break;
						case 'V': val1 = ValueFromVector(&((Vector *)ptr)[j]); break;
						case 'T': val1 = ValueFromTransform(&((Transform *)ptr)[j]); break;
					}
					rb_ary_push(val, val1);
				}
				break;
			}
			case 'G':
				val = ValueFromIntGroup(info->action->args[i].u.igval);
				break;
			case 'M':
				val = ValueFromMolecule(info->action->args[i].u.mval);
				break;
			case 'r':
			case 'b':
				val = (VALUE)info->action->args[i].u.vval;
				break;
			case 'R':
				rb_ary_concat(info->ary, (VALUE)info->action->args[i].u.vval);
				continue;  /*  Skip rb_ary_push()  */
			default:
				rb_raise(rb_eArgError, "Unsupported argument type '%c'", info->action->args[i].type);
		}
		rb_ary_push(info->ary, val);
	}

	if (info->mol != NULL) {
		info->receiver = ValueFromMolecule(info->mol);
	} else {
		info->receiver = rb_cMolecule;  /*  Assume class method  */
	}
	
	s = info->action->args[0].u.arval.ptr;
	for (p = s; *p != 0; p++) {
		if (isalpha(*p) || *p == '_' || (p > s && isdigit(*p)) || (p[1] == 0 && (*p == '?' || *p == '!')))
			continue;
		break;
	}
	
	if (*p == 0) {
		/*  Can be a method name  */
		info->method_id = rb_intern(s);
	} else {
		/*  Cannot be a method name: try to create a proc object and call it  */
		VALUE procstr = rb_str_new2(s);
		VALUE proc = rb_obj_instance_eval(1, &procstr, info->receiver);
		info->receiver = proc;
		info->method_id = rb_intern("call");
	}

	return Qnil;
}

static void
s_MolActionStoreReturnValue(MolRubyActionInfo *info, VALUE val)
{
	if (info->return_type != 0 && info->return_ptr != NULL) {
		switch (info->return_type) {
			case 'i': *((Int *)(info->return_ptr)) = NUM2INT(rb_Integer(val)); break;
			case 'd': *((Double *)(info->return_ptr)) = NUM2DBL(rb_Float(val)); break;
			case 's': *((char **)(info->return_ptr)) = strdup(StringValuePtr(val)); break;
			case 'v': VectorFromValue(val, (Vector *)(info->return_ptr)); break;
			case 't': TransformFromValue(val, (Transform *)(info->return_ptr)); break;
			case 'r': *((RubyValue *)(info->return_ptr)) = (RubyValue)val; break;
		}
	}
}

static VALUE
s_MolActionPerformRubyScript(VALUE vinfo)
{
/*	VALUE save_interrupt_flag; */
	VALUE val;
	MolRubyActionInfo *info = (MolRubyActionInfo *)vinfo;
/*	if (info->enable_interrupt)
		save_interrupt_flag = Ruby_SetInterruptFlag(Qtrue); */
	val = rb_funcall2(info->receiver, info->method_id, RARRAY_LEN(info->ary), RARRAY_PTR(info->ary));
	s_MolActionStoreReturnValue(info, val);
/*	if (info->enable_interrupt)
		Ruby_SetInterruptFlag(save_interrupt_flag); */
	return val;
}

static void
sMolActionUpdateSelectionAndParameterNumbering(Molecule *mol, const IntGroup *ig, int is_insert)
{
	int i, j, n, old_natoms;
	IntGroup *sel, *orig_atoms, *ig3;
	UnionPar *up;
	MolAction *act;

	if (ig == NULL || (n = IntGroupGetCount(ig)) == 0)
		return;
	
	/*  Insert: orig_atoms = set of original atoms in new indices  */
	/*  Delete: orig_atoms = set of remaining atoms in old indices  */
	orig_atoms = IntGroupNew();
	old_natoms = (is_insert ? mol->natoms - n : mol->natoms + n);
	ig3 = IntGroupNewWithPoints(0, (is_insert ? mol->natoms : old_natoms), -1);
	IntGroupXor(ig, ig3, orig_atoms);
	IntGroupRelease(ig3);

	/*  Update selection  */
	sel = MoleculeGetSelection(mol);
	if (sel != NULL && IntGroupGetCount(sel) > 0) {
		IntGroup *sel2 = IntGroupNew();
		if (is_insert)
			IntGroupConvolute(sel, orig_atoms, sel2);  /*  Renumber  */
		else {
			IntGroup *selx = IntGroupNew();
			IntGroupIntersect(sel, orig_atoms, selx);
			IntGroupDeconvolute(selx, orig_atoms, sel2);
			IntGroupRelease(selx);
		}
		IntGroupRetain(sel);  /*  To avoid being released during MoleculeSetSelection()  */
		MoleculeSetSelection(mol, sel2);
		if (IntGroupGetCount(sel) != IntGroupGetCount(sel2)) {
			/*  Register undo action  */
			act = MolActionNew(gMolActionSetSelection, sel);
			MolActionCallback_registerUndo(mol, act);
			MolActionRelease(act);
		}
		IntGroupRelease(sel2);
		IntGroupRelease(sel);
	}
	
	/*  Update parameters  */
	if (mol->par != NULL) {
		Int *ip = (Int *)malloc(sizeof(Int) * old_natoms);
		if (is_insert) {
			for (i = 0; i < old_natoms; i++)
				ip[i] = IntGroupGetNthPoint(orig_atoms, i);
		} else {
			for (i = 0; i < old_natoms; i++)
				ip[i] = IntGroupLookupPoint(orig_atoms, i);
		}
		for (i = kFirstParType; i <= kLastParType; i++) {
			UnionPar usave;
			UnionPar *upary = NULL;
			Int count_upary = 0;
			if (!is_insert)
				ig3 = IntGroupNew();
			for (j = 0; (up = ParameterGetUnionParFromTypeAndIndex(mol->par, i, j)) != NULL; j++) {
				usave = *up;
				if (ParameterRenumberAtoms(i, up, old_natoms, ip) && !is_insert) {
					IntGroupAdd(ig3, j, 1);  /*  This parameter is to be restored on undo  */
					AssignArray(&upary, &count_upary, sizeof(UnionPar), count_upary, &usave);
				}
			}
			if (count_upary > 0) {
				/*  Register undo for modifying parameters  */
				act = MolActionNew(gMolActionAddParameters, i, ig3, count_upary, upary);
				MolActionCallback_registerUndo(mol, act);
				MolActionRelease(act);
				act = MolActionNew(gMolActionDeleteParameters, i, ig3);
				MolActionCallback_registerUndo(mol, act);
				MolActionRelease(act);
				free(upary);
			}
			if (!is_insert)
				IntGroupRelease(ig3);
		}
		free(ip);
	}
	IntGroupRelease(orig_atoms);
}

int
MolActionPerform(Molecule *mol, MolAction *action)
{
	int n1, result, natoms;
	Molecule *mol2;
	IntGroup *ig;
	MolAction *act2 = NULL;
	int needsSymmetryAmendment = 0;
	int needsRebuildMDArena = 0;
	Int *ip;
	Vector *vp, v;
	n1 = 0;
	
	if (action == NULL)
		return -1;
	
	if (action->frame >= 0 && mol->cframe != action->frame)
		MoleculeSelectFrame(mol, action->frame, 1);

	/*  Ruby script execution  */
	if (strncmp(action->name, kMolActionPerformScript, strlen(kMolActionPerformScript)) == 0) {
		MolRubyActionInfo info;
		memset(&info, 0, sizeof(info));
		if (gMolbyIsCheckingInterrupt) {
			MolActionAlertRubyIsRunning();
			return -1;
		}
	/*	gMolbyBacktrace = Qnil; */
		info.mol = mol; /*  May be NULL  */
		info.action = action;
		gMolbyRunLevel++;
		rb_protect(s_MolActionToRubyArguments, (VALUE)&info, &result);
		gMolbyRunLevel--;
		if (result == 0) {
			VALUE save_interrupt;
			MyAppCallback_beginUndoGrouping();
		/*	info.enable_interrupt = 1; */
			save_interrupt = Ruby_SetInterruptFlag(Qtrue);
			gMolbyRunLevel++;
			rb_protect(s_MolActionPerformRubyScript, (VALUE)&info, &result);
			gMolbyRunLevel--;
			Ruby_SetInterruptFlag(save_interrupt);
			MyAppCallback_endUndoGrouping();
		}
		if (result != 0) {
			Molby_showError(result);
		}
		MyAppCallback_hideProgressPanel();  /*  In case when the progress panel is still onscreen */
		return (result == 0 ? 0 : -1);
	}
	
	if (mol == NULL)
		return -1;
	natoms = mol->natoms;
	
	if (strcmp(action->name, gMolActionAddAnAtom) == 0) {
		n1 = MoleculeCreateAnAtom(mol, action->args[0].u.aval, action->args[1].u.ival);
		if ((ip = action->args[2].u.pval) != NULL)
			*ip = n1;
		if (n1 < 0)
			return -1;
		ig = IntGroupNewWithPoints(n1, 1, -1);

		/*  Update selection  */
		sMolActionUpdateSelectionAndParameterNumbering(mol, ig, 1);

		act2 = MolActionNew(gMolActionDeleteAnAtom, n1);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionMergeMolecule) == 0) {
		IntGroup *ig2;
		Molecule *mol2;
		ig = action->args[1].u.igval;
		mol2 = action->args[0].u.mval;
		if (ig != NULL) {
			ig2 = ig;
			IntGroupRetain(ig2);
		} else {
			ig2 = IntGroupNew();
			IntGroupAdd(ig2, mol->natoms, mol2->natoms);
		}
		/*  Calculate the offset for the residue number  */
		{
			int minreg, maxreg;
			Atom *ap;
			maxreg = -1;
			for (n1 = 0, ap = mol->atoms; n1 < mol->natoms; n1++, ap = ATOM_NEXT(ap)) {
				if (ap->resSeq > maxreg)
					maxreg = ap->resSeq;
			}
			minreg = ATOMS_MAX_NUMBER;
			for (n1 = 0, ap = mol2->atoms; n1 < mol2->natoms; n1++, ap = ATOM_NEXT(ap)) {
				if (ap->resSeq < minreg)
					minreg = ap->resSeq;
			}
			if (maxreg < 0)
				maxreg = 0;
			if (minreg == ATOMS_MAX_NUMBER)
				minreg = 0;
			if (maxreg >= minreg)
				n1 = maxreg - minreg + 1;
			else n1 = 0;
		}
		if ((result = MoleculeMerge(mol, mol2, ig, n1)) != 0)
			return result;

		sMolActionUpdateSelectionAndParameterNumbering(mol, ig, 1);

		act2 = MolActionNew(gMolActionUnmergeMolecule, ig2);
		IntGroupRelease(ig2);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteAnAtom) == 0 || (strcmp(action->name, gMolActionUnmergeMolecule) == 0 && (n1 = 1))) {
		IntGroup *bg;
		if (n1 == 0) {
			ig = IntGroupNew();
			if (ig == NULL || IntGroupAdd(ig, action->args[0].u.ival, 1) != 0)
				return -1;
		} else {
			ig = action->args[0].u.igval;
			IntGroupRetain(ig);
		}
		/*  Search bonds crossing the molecule border  */
		bg = MoleculeSearchBondsAcrossAtomGroup(mol, ig);
		ip = NULL;
		if (bg != NULL) {
			n1 = IntGroupGetCount(bg);
			if (n1 > 0) {
				IntGroupIterator iter;
				Int i, idx;
				ip = (Int *)calloc(sizeof(Int), n1 * 2 + 1);
				if (ip == NULL)
					return -1;
				IntGroupIteratorInit(bg, &iter);
				idx = 0;
				while ((i = IntGroupIteratorNext(&iter)) >= 0) {
					/*  The atoms at the border  */
					ip[idx++] = mol->bonds[i * 2];
					ip[idx++] = mol->bonds[i * 2 + 1];
				}
				IntGroupIteratorRelease(&iter);
			}
			IntGroupRelease(bg);
		}
		/*  Unmerge molecule  */
		if ((result = MoleculeUnmerge(mol, &mol2, ig, 0)) != 0) {
			if (ip != NULL)
				free(ip);
			return result;
		}
		
		sMolActionUpdateSelectionAndParameterNumbering(mol, ig, 0);

#if 0
		/*  Update selection  */
		{
			IntGroup *sel = MoleculeGetSelection(mol);
			if (sel != NULL && IntGroupGetCount(sel) > 0) {
				IntGroup *remain_atoms = IntGroupNew();  /*  Remaining atoms (with original indices)  */
				IntGroup *sel2 = IntGroupNew();
				IntGroup *sel3 = IntGroupNew();
				IntGroupXor(ig, IntGroupNewWithPoints(0, mol->natoms + IntGroupGetCount(ig), -1), remain_atoms);
				IntGroupIntersect(sel, remain_atoms, sel2);  /*  Atoms to select after deletion  */
				IntGroupDeconvolute(sel2, remain_atoms, sel3);  /*  Renumber  */
				if (!IntGroupIsEqual(sel, sel3)) {
					MoleculeSetSelection(mol, sel3);
					act2 = MolActionNew(gMolActionSetSelection, sel);
					MolActionCallback_registerUndo(mol, act2);
					act2 = NULL;
				}
				IntGroupRelease(sel3);
				IntGroupRelease(sel2);
				IntGroupRelease(remain_atoms);
			}
		}
#endif
		if (mol2 == NULL)
			act2 = NULL;
		else {
			/*  If there exist bonds crossing the molecule border, then register
			    an action to restore them  */
			if (ip != NULL) {
				act2 = MolActionNew(gMolActionAddBonds, n1 * 2, ip);
				MolActionCallback_registerUndo(mol, act2);
				free(ip);
			}
			act2 = MolActionNew(gMolActionMergeMolecule, mol2, ig);
			MoleculeRelease(mol2);
		}
		IntGroupRelease(ig);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddBonds) == 0) {
		ip = (Int *)action->args[0].u.arval.ptr;
		n1 = action->args[0].u.arval.nitems / 2;
		if ((result = MoleculeAddBonds(mol, n1, ip)) <= 0)
			return result;
		ip = (Int *)malloc(sizeof(Int) * 2 * result);
		if (ip == NULL)
			return -4;
		memmove(ip, mol->bonds - result * 2, sizeof(Int) * 2 * result);
		act2 = MolActionNew(gMolActionDeleteBonds, result * 2, ip);
		free(ip);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteBonds) == 0) {
		ip = (Int *)action->args[0].u.arval.ptr;
		n1 = action->args[0].u.arval.nitems / 2;
		if ((result = MoleculeDeleteBonds(mol, n1, ip)) < 0)
			return result;
		act2 = MolActionNew(gMolActionAddBonds, n1 * 2, ip);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddAngles) == 0) {
		ip = (Int *)action->args[0].u.arval.ptr;
		n1 = action->args[0].u.arval.nitems / 3;
		ig = action->args[1].u.igval;
		if (ig == NULL)
			ig = IntGroupNewWithPoints(mol->nangles, n1, -1);
		else
			IntGroupRetain(ig);
		if ((result = MoleculeAddAngles(mol, ip, ig)) < 0) {
			IntGroupRelease(ig);
			return result;
		}
		act2 = MolActionNew(gMolActionDeleteAngles, ig);
		IntGroupRelease(ig);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteAngles) == 0) {
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig) * 3;
		ip = (Int *)malloc(sizeof(Int) * n1);
		if ((result = MoleculeDeleteAngles(mol, ip, ig)) < 0) {
			free(ip);
			return result;
		}
		act2 = MolActionNew(gMolActionAddAngles, n1, ip, ig);
		free(ip);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddDihedrals) == 0) {
		ip = (Int *)action->args[0].u.arval.ptr;
		n1 = action->args[0].u.arval.nitems / 4;
		ig = action->args[1].u.igval;
		if (ig == NULL)
			ig = IntGroupNewWithPoints(mol->ndihedrals, n1, -1);
		else
			IntGroupRetain(ig);
		if ((result = MoleculeAddDihedrals(mol, ip, ig)) < 0) {
			IntGroupRelease(ig);
			return result;
		}
		act2 = MolActionNew(gMolActionDeleteDihedrals, ig);
		IntGroupRelease(ig);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteDihedrals) == 0) {
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig) * 4;
		ip = (Int *)malloc(sizeof(Int) * n1);
		if ((result = MoleculeDeleteDihedrals(mol, ip, ig)) < 0) {
			free(ip);
			return result;
		}
		act2 = MolActionNew(gMolActionAddDihedrals, n1, ip, ig);
		free(ip);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddImpropers) == 0) {
		ip = (Int *)action->args[0].u.arval.ptr;
		n1 = action->args[0].u.arval.nitems / 4;
		ig = action->args[1].u.igval;
		if (ig == NULL)
			ig = IntGroupNewWithPoints(mol->nimpropers, n1, -1);
		else
			IntGroupRetain(ig);
		if ((result = MoleculeAddImpropers(mol, ip, ig)) < 0) {
			IntGroupRelease(ig);
			return result;
		}
		act2 = MolActionNew(gMolActionDeleteImpropers, ig);
		IntGroupRelease(ig);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteImpropers) == 0) {
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig) * 4;
		ip = (Int *)malloc(sizeof(Int) * n1);
		if ((result = MoleculeDeleteImpropers(mol, ip, ig)) < 0) {
			free(ip);
			return result;
		}
		act2 = MolActionNew(gMolActionAddImpropers, n1, ip, ig);
		free(ip);
		needsRebuildMDArena = 1;
/*	} else if (strcmp(action->name, gMolActionReplaceTables) == 0) {
		Int i1, i2, i3;
		Int *ip1, *ip2, *ip3;
		if ((i1 = MoleculeReplaceAllAngles(mol, action->args[0].u.arval.nitems / 3, (Int *)action->args[0].u.arval.ptr, &ip1)) < 0)
			return -1;
		if ((i2 = MoleculeReplaceAllDihedrals(mol, action->args[1].u.arval.nitems / 4, (Int *)action->args[1].u.arval.ptr, &ip2)) < 0)
			return -1;
		if ((i3 = MoleculeReplaceAllImpropers(mol, action->args[2].u.arval.nitems / 4, (Int *)action->args[2].u.arval.ptr, &ip3)) < 0)
			return -1;
		act2 = MolActionNew(gMolActionReplaceTables, i1, ip1, i2, ip2, i3, ip3);
		free(ip1);
		free(ip2);
		free(ip3);
		needsRebuildMDArena = 1; */
	} else if (strcmp(action->name, gMolActionTranslateAtoms) == 0) {
		vp = (Vector *)action->args[0].u.arval.ptr;
		ig = action->args[1].u.igval;
		if (vp == NULL)
			return -1;
		MoleculeTranslate(mol, vp, ig);
		VecScale(v, *vp, -1);
		act2 = MolActionNew(gMolActionTranslateAtoms, &v, ig);
		act2->frame = mol->cframe;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionRotateAtoms) == 0) {
		Vector *vp2;
		Double ang;
		vp = (Vector *)action->args[0].u.arval.ptr;
		ang = action->args[1].u.dval;
		vp2 = (Vector *)action->args[2].u.arval.ptr;
		ig = action->args[3].u.igval;
		MoleculeRotate(mol, vp, ang, vp2, ig);
		act2 = MolActionNew(gMolActionRotateAtoms, vp, -ang, vp2, ig);
		act2->frame = mol->cframe;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionTransformAtoms) == 0) {
		Transform *trp, tr_inv;
		int atomPositions = 0; /* Save atom positions for undo? */
		trp = (Transform *)action->args[0].u.arval.ptr;
		ig = action->args[1].u.igval;
		if (TransformInvert(tr_inv, *trp) != 0)
			atomPositions = 1;
		else {
			/*  Check if inverse transform is reliable enough  */
			Transform temp;
			int j;
			TransformMul(temp, *trp, tr_inv);
			temp[0] -= 1.0;
			temp[4] -= 1.0;
			temp[8] -= 1.0;
			for (j = 0; j < 12; j++) {
				if (temp[j] < -1e-4 || temp[j] > 1e-4)
					break;
			}
			if (j != 12)
				atomPositions = 1;
		}
		if (atomPositions) {
			IntGroupIterator iter;
			int j, k;
			n1 = IntGroupGetCount(ig);
			vp = (Vector *)malloc(sizeof(Vector) * n1);
			if (vp == NULL)
				return -1;
			IntGroupIteratorInit(ig, &iter);
			k = 0;
			while ((j = IntGroupIteratorNext(&iter)) >= 0) {
				vp[k++] = (ATOM_AT_INDEX(mol->atoms, j))->r;
			}
			act2 = MolActionNew(gMolActionSetAtomPositions, ig, n1, vp);
			free(vp);
			IntGroupIteratorRelease(&iter);
		} else {
			act2 = MolActionNew(gMolActionTransformAtoms, &tr_inv, ig);
		}
		MoleculeTransform(mol, *trp, ig);
		act2->frame = mol->cframe;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetAtomPositions) == 0) {
		IntGroupIterator iter;
		int j, k;
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig);
		vp = (Vector *)malloc(sizeof(Vector) * n1);
		if (vp == NULL)
			return -1;
		IntGroupIteratorInit(ig, &iter);
		k = 0;
		while ((j = IntGroupIteratorNext(&iter)) >= 0) {
			vp[k++] = (ATOM_AT_INDEX(mol->atoms, j))->r;
		}
		act2 = MolActionNew(gMolActionSetAtomPositions, ig, n1, vp);
		free(vp);
		vp = (Vector *)action->args[1].u.arval.ptr;
		IntGroupIteratorReset(&iter);
		k = 0;
		while ((j = IntGroupIteratorNext(&iter)) >= 0) {
			(ATOM_AT_INDEX(mol->atoms, j))->r = vp[k++];
		}
		IntGroupIteratorRelease(&iter);
		act2->frame = mol->cframe;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionInsertFrames) == 0) {
		int old_nframes, new_nframes;
		ig = action->args[0].u.igval;
		vp = (Vector *)action->args[1].u.arval.ptr;
		n1 = IntGroupGetCount(ig);
		if (n1 == 0)
			return 0;  /*  Do nothing  */
		if (vp != NULL && action->args[1].u.arval.nitems != n1 * mol->natoms)
			return -1;  /*  Internal inconsistency  */
		old_nframes = MoleculeGetNumberOfFrames(mol);
		if (MoleculeInsertFrames(mol, ig, vp) < 0)
			return -1;  /*  Error  */
		new_nframes = MoleculeGetNumberOfFrames(mol);
		if (old_nframes + n1 < new_nframes) {
			/*  "Extra" frames were automatically inserted because large frame indices were specified  */
			/*  Register undo operation to remove these extra frames  */
			IntGroup *ig2 = IntGroupNewWithPoints(old_nframes, new_nframes - (old_nframes + n1), -1);
			act2 = MolActionNew(gMolActionRemoveFrames, ig2);
			IntGroupRelease(ig2);
			MolActionCallback_registerUndo(mol, act2);
			MolActionRelease(act2);
		}			
		act2 = MolActionNew(gMolActionRemoveFrames, ig);
		act2->frame = mol->cframe;
	} else if (strcmp(action->name, gMolActionRemoveFrames) == 0) {
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig);
		if (n1 == 0)
			return 0;  /*  Do nothing  */
		vp = (Vector *)calloc(sizeof(Vector), n1 * mol->natoms);
		if (MoleculeRemoveFrames(mol, ig, vp) < 0)
			return -1;  /*  Error  */
		act2 = MolActionNew(gMolActionInsertFrames, ig, n1 * mol->natoms, vp);
		act2->frame = mol->cframe;
		free(vp);
	} else if (strcmp(action->name, gMolActionSetSelection) == 0) {
		IntGroup *ig2;
		ig2 = MoleculeGetSelection(mol);
		if (ig2 != NULL)
			IntGroupRetain(ig2);  /*  To avoid releasing during MoleculeSetSelection() */
		ig = action->args[0].u.igval;
		MoleculeSetSelection(mol, ig);
		act2 = MolActionNew(gMolActionSetSelection, ig2);
		if (ig2 != NULL)
			IntGroupRelease(ig2);
	} else if (strcmp(action->name, gMolActionReorderAtoms) == 0) {
		Int *ip2 = (Int *)malloc(sizeof(Int) * (mol->natoms + 1));
		if (ip2 == NULL)
			return -1;
		ip = (Int *)action->args[0].u.arval.ptr;
		n1 = action->args[0].u.arval.nitems;
	/*	for (n1 = 0; n1 < mol->natoms; n1++) {
			if (ip[n1] < 0)
				break;
		} */
		result = MoleculeReorderAtoms(mol, ip, ip2, n1);
		if (result != 0) {
			free(ip2);
			return result;
		}
		act2 = MolActionNew(gMolActionReorderAtoms, mol->natoms, ip2);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionChangeResidueNumber) == 0 || (strcmp(action->name, gMolActionChangeResidueNumberForUndo) == 0 && (n1 = 1))) {
		IntGroupIterator iter;
		int i, nresidues;
		int forUndo;
		forUndo = n1;
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig);
		ip = (Int *)calloc(sizeof(Int), n1 + 1);
		IntGroupIteratorInit(ig, &iter);
		i = 0;
		while ((n1 = IntGroupIteratorNext(&iter)) >= 0) {
			ip[i++] = ATOM_AT_INDEX(mol->atoms, n1)->resSeq;
		}
		n1 = i;
		ip[i] = kInvalidIndex;
		nresidues = mol->nresidues;
		if (forUndo) {
			MoleculeChangeResidueNumberWithArray(mol, ig, (Int *)action->args[1].u.arval.ptr);
			MoleculeChangeNumberOfResidues(mol, action->args[2].u.ival);
		} else {
			MoleculeChangeResidueNumber(mol, ig, action->args[1].u.ival);
		}
		act2 = MolActionNew(gMolActionChangeResidueNumberForUndo, ig, n1, ip, nresidues);
	} else if (strcmp(action->name, gMolActionChangeResidueNames) == 0) {
		char *new_names, *old_names;
		int i, argc;
		ip = (Int *)action->args[0].u.arval.ptr;
		argc = action->args[0].u.arval.nitems;
		new_names = action->args[1].u.arval.ptr;
		old_names = (char *)malloc(argc * 4 + 1);
		for (i = 0; i < argc; i++) {
			if (ip[i] >= mol->nresidues)
				old_names[i * 4] = 0;
			else
				strncpy(old_names + i * 4, mol->residues[ip[i]], 4);
		}
		MoleculeChangeResidueNames(mol, argc, ip, new_names);
		act2 = MolActionNew(gMolActionChangeResidueNames, argc, ip, argc * 4, old_names);
		free(old_names);
	} else if (strcmp(action->name, gMolActionOffsetResidueNumbers) == 0) {
		ig = action->args[0].u.igval;
		n1 = mol->nresidues;
		result = MoleculeOffsetResidueNumbers(mol, ig, action->args[1].u.ival, action->args[2].u.ival);
		if (result != 0)
			return result;  /*  The molecule is not modified  */
		act2 = MolActionNew(gMolActionOffsetResidueNumbers, ig, -(action->args[1].u.ival), n1);
	} else if (strcmp(action->name, gMolActionChangeNumberOfResidues) == 0) {
		int nresidues = mol->nresidues;
		n1 = action->args[0].u.ival;
		if (n1 < nresidues) {
			/*  The residue names will be lost, so undo must be registered to restore the names  */
			int argc = nresidues - n1;
			char *names = (char *)malloc(4 * argc + 1);
			Int *ip = (Int *)malloc(sizeof(Int) * (argc + 1));
			int i;
			for (i = 0; i < argc; i++) {
				strncpy(names + i * 4, mol->residues[i + n1], 4);
				ip[i] = i + n1;
			}
			ip[i] = kInvalidIndex;
			act2 = MolActionNew(gMolActionChangeResidueNames, argc, ip, argc * 4, names);
			MolActionCallback_registerUndo(mol, act2);
			free(ip);
			free(names);
		}
		MoleculeChangeNumberOfResidues(mol, n1);
		act2 = MolActionNew(gMolActionChangeNumberOfResidues, nresidues);
	} else if (strcmp(action->name, gMolActionExpandBySymmetry) == 0) {
		Symop symop;
		ig = action->args[0].u.igval;
		symop.dx = action->args[1].u.ival;
		symop.dy = action->args[2].u.ival;
		symop.dz = action->args[3].u.ival;
		symop.sym = action->args[4].u.ival;
		symop.alive = (symop.dx != 0 || symop.dy != 0 || symop.dz != 0 || symop.sym != 0);
		n1 = MoleculeAddExpandedAtoms(mol, symop, ig);
		if (n1 > 0) {
			ig = IntGroupNew();
			IntGroupAdd(ig, mol->natoms - n1, n1);
			act2 = MolActionNew(gMolActionUnmergeMolecule, ig);
			IntGroupRelease(ig);
		} else return n1;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddSymmetryOperation) == 0) {
		Transform *trp;
		Transform itr = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
		trp = (Transform *)action->args[0].u.arval.ptr;
		if (mol->nsyms == 0) {
			for (n1 = 0; n1 < 12; n1++) {
				if (fabs((*trp)[n1] - itr[n1]) > 1e-8)
					break;
			}
			if (n1 < 12) {
				if (AssignArray(&mol->syms, &mol->nsyms, sizeof(Transform), mol->nsyms, &itr) == 0)
					return -1;
				act2 = MolActionNew(gMolActionDeleteSymmetryOperation);
				MolActionCallback_registerUndo(mol, act2);
				MolActionRelease(act2);
			}
		}
		if (AssignArray(&mol->syms, &mol->nsyms, sizeof(Transform), mol->nsyms, trp) == 0)
				return -1;
		act2 = MolActionNew(gMolActionDeleteSymmetryOperation);
	} else if (strcmp(action->name, gMolActionDeleteSymmetryOperation) == 0) {
		if (mol->nsyms == 0)
			return -1;
		act2 = MolActionNew(gMolActionAddSymmetryOperation, &(mol->syms[mol->nsyms - 1]));
		mol->nsyms--;
	/*	if (mol->nsyms == 1)
			mol->nsyms--;  *//*  Remove the identity operation  */
		if (mol->nsyms == 0) {
			free(mol->syms);
			mol->syms = NULL;
		}
	} else if (strcmp(action->name, gMolActionSetCell) == 0) {
		double *dp, d[6];
		if (mol->cell == NULL)
			d[0] = 0.0;
		else {
			for (n1 = 0; n1 < 6; n1++)
				d[n1] = mol->cell->cell[n1];
		}
		n1 = action->args[1].u.ival;
		dp = action->args[0].u.arval.ptr;
		if (action->args[0].u.arval.nitems == 0)
			MoleculeSetCell(mol, 0, 0, 0, 0, 0, 0, n1);
		else
			MoleculeSetCell(mol, dp[0], dp[1], dp[2], dp[3], dp[4], dp[5], n1);
		act2 = MolActionNew(gMolActionSetCell, (d[0] == 0.0 ? 0 : 6), d, n1);
	} else if (strcmp(action->name, gMolActionSetBox) == 0) {
		Vector v[4];
		char flags[3];
		if (mol->cell == NULL)
			act2 = MolActionNew(gMolActionClearBox);
		else {
			n1 = ((mol->cell->flags[0] != 0) * 4 + (mol->cell->flags[1] != 0) * 2 + (mol->cell->flags[2] != 0));
			act2 = MolActionNew(gMolActionSetBox, &(mol->cell->axes[0]), &(mol->cell->axes[1]), &(mol->cell->axes[2]), &(mol->cell->origin), n1);
		}
		for (n1 = 0; n1 < 4; n1++)
			v[n1] = *((Vector *)(action->args[n1].u.arval.ptr));
		for (n1 = 0; n1 < 3; n1++)
			flags[n1] = ((action->args[4].u.ival >> (2 - n1)) & 1);
		MoleculeSetPeriodicBox(mol, &v[0], &v[1], &v[2], &v[3], flags);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionClearBox) == 0) {
		if (mol->cell == NULL)
			return 0;  /*  Do nothing  */
		n1 = ((mol->cell->flags[0] != 0) * 4 + (mol->cell->flags[1] != 0) * 2 + (mol->cell->flags[2] != 0));
		act2 = MolActionNew(gMolActionSetBox, &(mol->cell->axes[0]), &(mol->cell->axes[1]), &(mol->cell->axes[2]), &(mol->cell->origin), n1);
		MoleculeSetPeriodicBox(mol, NULL, NULL, NULL, NULL, NULL);
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddParameters) == 0) {
		UnionPar *up;
		Int parType;
		if (mol->par == NULL)
			return -1;
		parType = action->args[0].u.ival;
		ig = action->args[1].u.igval;
		up = action->args[2].u.arval.ptr;
		ParameterInsert(mol->par, parType, up, ig);
		act2 = MolActionNew(gMolActionDeleteParameters, parType, ig);
	} else if (strcmp(action->name, gMolActionDeleteParameters) == 0) {
		UnionPar *up;
		Int parType;
		if (mol->par == NULL)
			return -1;
		parType = action->args[0].u.ival;
		ig = action->args[1].u.igval;
		n1 = IntGroupGetCount(ig);
		up = (UnionPar *)calloc(sizeof(UnionPar), n1);
		ParameterDelete(mol->par, parType, up, ig);
		act2 = MolActionNew(gMolActionAddParameters, parType, ig, n1, up);
		free(up);
/*	} else if (strcmp(action->name, gMolActionCartesianToXtal) == 0 || (strcmp(action->name, gMolActionXtalToCartesian) == 0 && (n1 = 1))) {
		Int i, j;
		Atom *ap;
		if (mol->cell == NULL || (n1 == 0 && mol->is_xtal_coord) || (n1 == 1 && !mol->is_xtal_coord))
			return 0;
		if (n1 == 0) {
			for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
				TransformVec(&(ap->r), mol->cell->tr, &(ap->r));
				if (ap->nframes > 0) {
					for (j = 0; j < ap->nframes; j++)
						TransformVec(&(ap->frames[j]), mol->cell->tr, &(ap->frames[j]));
				}
			}
			mol->is_xtal_coord = 1;
			act2 = MolActionNew(gMolActionXtalToCartesian);
		} else {
			for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap))
				TransformVec(&(ap->r), mol->cell->rtr, &(ap->r));
			if (ap->nframes > 0) {
				for (j = 0; j < ap->nframes; j++)
					TransformVec(&(ap->frames[j]), mol->cell->rtr, &(ap->frames[j]));
			}
			mol->is_xtal_coord = 0;
			act2 = MolActionNew(gMolActionCartesianToXtal);
		} */
	} else {
		fprintf(stderr, "Internal error: unknown action name %s\n", action->name);
		return -1;
	}
	if (act2 != NULL) {
		MolActionCallback_registerUndo(mol, act2);
		MolActionRelease(act2);
	}
	if (needsSymmetryAmendment) {
		//  ig should be valid
		IntGroup *ig2;
		Vector *vp2;
		n1 = MoleculeAmendBySymmetry(mol, ig, &ig2, &vp2);
		if (n1 > 0) {
			act2 = MolActionNew(gMolActionSetAtomPositions, ig2, n1, vp2);
			act2->frame = mol->cframe;
			MolActionCallback_registerUndo(mol, act2);
			MolActionRelease(act2);
			free(vp2);
			IntGroupRelease(ig2);
		}
	}
	
	if (needsRebuildMDArena) {
		mol->needsMDRebuild = 1;
	}
	
	MoleculeCallback_notifyModification(mol, 0);
	return 0;
}

int
MolActionCreateAndPerform(Molecule *mol, const char *name, ...)
{
	va_list ap;
	MolAction *action;
	va_start(ap, name);
	action = MolActionNewArgv(name, ap);
	va_end(ap);
	if (action != NULL) {
		int result = MolActionPerform(mol, action);
		MolActionRelease(action);
		return result;
	} else return -1;
}

void
MolActionAlertRubyIsRunning(void)
{
	MyAppCallback_errorMessageBox("Cannot perform operation (Ruby script is running background)");
}
