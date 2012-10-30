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
const char *gMolActionMergeMoleculeForUndo = "mergeMolForUndo:MG";
const char *gMolActionUnmergeMolecule = "unmergeMol:G";
const char *gMolActionUnmergeMoleculeForUndo = "unmergeMolForUndo:G";
const char *gMolActionAddBonds        = "addBonds:IG";
const char *gMolActionAddBondsForUndo = "addBondsForUndo:IG";
const char *gMolActionDeleteBonds     = "deleteBonds:G";
const char *gMolActionAddAngles       = "addAngles:IG";
const char *gMolActionDeleteAngles    = "deleteAngles:G";
const char *gMolActionAddDihedrals    = "addDihedrals:IG";
const char *gMolActionDeleteDihedrals = "deleteDihedrals:G";
const char *gMolActionAddImpropers    = "addImpropers:IG";
const char *gMolActionDeleteImpropers = "deleteImpropers:G";
const char *gMolActionTranslateAtoms  = "translateAtoms:vG";
const char *gMolActionRotateAtoms     = "rotate:vdvG";
const char *gMolActionTransformAtoms  = "transform:tG";
const char *gMolActionSetAtomPositions = "atomPosition:GV";
const char *gMolActionSetAtomVelocities = "atomVelocity:GV";
const char *gMolActionSetAtomForces   = "atomForce:GV";
const char *gMolActionInsertFrames    = "insertFrames:GVV";
const char *gMolActionRemoveFrames    = "removeFrames:G";
const char *gMolActionSetSelection    = "selection:G";
const char *gMolActionChangeResidueNumber = "changeResSeq:Gi";
const char *gMolActionChangeResidueNumberForUndo = "changeResSeqForUndo:GIi";
const char *gMolActionChangeResidueNames = "changeResNames:IC";
const char *gMolActionOffsetResidueNumbers = "offsetResSeq:Gii";
const char *gMolActionChangeNumberOfResidues = "changeNres:i";
const char *gMolActionRenumberAtoms    = "renumberAtoms:I";
const char *gMolActionExpandBySymmetry = "expandSym:Giiii;I";
const char *gMolActionDeleteSymmetryOperation = "deleteSymop";
const char *gMolActionAddSymmetryOperation = "addSymop:t";
const char *gMolActionSetCell         = "setCell:Di";
const char *gMolActionSetCellPeriodicity = "setCellPeriodicity:i";
const char *gMolActionSetBox          = "setBox:vvvvii";
const char *gMolActionClearBox        = "clearBox";
/*const char *gMolActionSetBoxForFrames = "setBoxForFrames:V"; */
/*const char *gMolActionSetCellFlexibility = "setCellFlexibility:i"; */
const char *gMolActionAddParameters   = "addParameters:iGU";
const char *gMolActionDeleteParameters = "deleteParameters:iG";
const char *gMolActionAmendBySymmetry = "amendBySymmetry:G;G";

/*  A Ruby array to retain objects used in MolActionArg  */
static VALUE sMolActionArgValues = Qfalse;

/*  Action arguments  */
/*  (Simple types)  i: Int, d: double, s: string, v: Vector, t: Transform, u: UnionPar
 (Array types)   I: array of Int, D: array of double, V: array of Vector, C: array of char, T: array of Transform, U: array of UnionPars
 (Complex types) M: Molecule, G: IntGroup, A: Atom
 (Ruby value)    b: Ruby boolean, r: Ruby object, R: an array of Ruby object (a Ruby array)
 (Return value)  i, d, s, v, t, r, G + 0x80 */
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
		struct {
			void *ptr;    /*  Return value pointer  */
			Int *nptr;    /*  If the return value is an array, then *ptr contains a malloc'ed pointer, and *nptr contains the number of items.  */
		} retval;          /*  Store address for the return value; usually used in performing Ruby script */
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
				if (*p == 'i' || *p == 'd' || *p == 's' || *p == 'v' || *p == 't' || *p == 'r' || *p == 'G') {
					arg.type = (*p | 0x80);
					p++;
					arg.u.retval.ptr = va_arg(ap, void *);
				} else if (*p == 'I' || *p == 'D' || *p == 'V') {
					arg.type = (*p | 0x80);
					p++;
					arg.u.retval.nptr = va_arg(ap, Int *);
					arg.u.retval.ptr = va_arg(ap, void *);
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
	Int *return_nptr;
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
			info->return_ptr = argp->u.retval.ptr;
			info->return_nptr = argp->u.retval.nptr;  /*  May be unused  */
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
			case 'G': *((IntGroup **)(info->return_ptr)) = (val == Qnil ? NULL : IntGroupFromValue(val)); break;
			case 'I':
			case 'D':
			case 'V': {
				int i, n, size;
				void *p;
				val = rb_ary_to_ary(val);
				n = RARRAY_LEN(val);
				*(info->return_nptr) = n;
				if (n == 0) {
					*((void **)(info->return_ptr)) = NULL;
					break;
				}
				if (info->return_type == 'I')
					size = sizeof(Int);
				else if (info->return_type == 'D')
					size = sizeof(Double);
				else if (info->return_type == 'V')
					size = sizeof(Vector);
				else break;
				*((void **)(info->return_ptr)) = p = calloc(size, n);
				for (i = 0; i < n; i++) {
					VALUE rval = (RARRAY_PTR(val))[i];
					if (info->return_type == 'I') {
						((Int *)p)[i] = NUM2INT(rb_Integer(rval));
					} else if (info->return_type == 'D') {
						((Double *)p)[i] = NUM2DBL(rb_Float(rval));
					} else if (info->return_type == 'V') {
						VectorFromValue(rval, (Vector *)p + i);
					}
				}
				break;
			}
			default: break;
		}
	}
}

static VALUE
s_MolActionExecuteRubyScript(VALUE vinfo)
{
	VALUE val;
	MolRubyActionInfo *info = (MolRubyActionInfo *)vinfo;
	val = rb_funcall2(info->receiver, info->method_id, RARRAY_LEN(info->ary), RARRAY_PTR(info->ary));
	s_MolActionStoreReturnValue(info, val);
	return val;
}

static int
s_MolActionPerformRubyScript(Molecule *mol, MolAction *action)
{
	int result;
	MolRubyActionInfo info;
	memset(&info, 0, sizeof(info));
	if (gMolbyIsCheckingInterrupt) {
		MolActionAlertRubyIsRunning();
		return -1;
	}
	info.mol = mol; /*  May be NULL  */
	info.action = action;
	gMolbyRunLevel++;
	rb_protect(s_MolActionToRubyArguments, (VALUE)&info, &result);
	gMolbyRunLevel--;
	if (result == 0) {
		VALUE save_interrupt;
		MyAppCallback_beginUndoGrouping();
		save_interrupt = Ruby_SetInterruptFlag(Qtrue);
		gMolbyRunLevel++;
		rb_protect(s_MolActionExecuteRubyScript, (VALUE)&info, &result);
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

static int
s_MolActionLog(Molecule *mol, MolAction *action, FILE *fp)
{
	/*  (Simple types)  i: Int, d: double, s: string, v: Vector, t: Transform, u: UnionPar
	 (Array types)   I: array of Int, D: array of double, V: array of Vector, C: array of char, T: array of Transform, U: array of UnionPars
	 (Complex types) M: Molecule, G: IntGroup, A: Atom
	 (Ruby value)    b: Ruby boolean, r: Ruby object, R: an array of Ruby object (a Ruby array)
	 (Return value)  i, d, s, v, t, r, G + 0x80 */
	
	int i, j, k, lastIval = 0;
	char buf[80];
	MolActionArg *argp;
	fprintf(fp, "MolAction: %s", action->name);
	for (i = 0; i < action->nargs; i++) {
		argp = action->args + i;
		if (argp->type >= 0x80)
			break;
		fprintf(fp, " ");
		switch (argp->type) {
			case 'i': lastIval = argp->u.ival; fprintf(fp, " %d", lastIval); break;
			case 'd': fprintf(fp, " %g", argp->u.dval); break;
			case 's':
			case 'C': fprintf(fp, " %.*s", argp->u.arval.nitems, argp->u.arval.ptr); break;
			case 'v': case 't': case 'u': case 'I': case 'D': case 'V': case 'T': case 'U': {
				char *pp = (char *)argp->u.arval.ptr;
				int n = argp->u.arval.nitems;
				if (argp->type == 'v' || argp->type == 't' || argp->type == 'u')
					n = 1;
				else fprintf(fp, "[");
				for (j = 0; j < n; j++) {
					switch (argp->type) {
						case 'I':
							fprintf(fp, "%d", *((Int *)pp));
							pp += sizeof(Int);
							break;
						case 'D':
							fprintf(fp, "%g", *((Double *)pp));
							pp += sizeof(Double);
							break;
						case 'v': case 'V':
							fprintf(fp, "Vector3d[%g,%g,%g]", ((Double *)pp)[0], ((Double *)pp)[1], ((Double *)pp)[2]);
							pp += sizeof(Double) * 3;
							break;
						case 't': case 'T':
							fprintf(fp, "Transform[");
							for (k = 0; k < 12; k++)
								fprintf(fp, "%g%c", ((Double *)pp)[k], (k == 11 ? ']' : ','));
							pp += sizeof(Double) * 12;
							break;
						case 'u': case 'U':
							fprintf(fp, "Parameter[");
							for (k = 0; k < 8; k++) {
								ParameterItemToString(lastIval, (UnionPar *)pp, k, buf, sizeof buf);
								if (buf[0] == 0)
									break;
								fprintf(fp, "%s%s", (k == 0 ? "" : ","), buf);
							}
							fprintf(fp, "]");
							pp += sizeof(UnionPar);
							break;
					}
					if (j < n - 1)
						fprintf(fp, ",");
				}
				if (argp->type == 'v' || argp->type == 't' || argp->type == 'u')
					n = 1;
				else fprintf(fp, "]");
				break;
			}
			case 'r': {
				VALUE val = rb_inspect(argp->u.vval);
				fprintf(fp, "%s", StringValuePtr(val));
				break;
			}
			case 'M':
				MoleculeCallback_displayName(argp->u.mval, buf, sizeof buf);
				if (buf[0] == 0) {
					/*  No associated document  */
					snprintf(buf, sizeof buf, "#<Molecule:0x%lx>", argp->u.mval);
				}
				fprintf(fp, "%s", buf);
				break;
			case 'G': {
				char *s = IntGroupInspect(argp->u.igval);
				fprintf(fp, "%s", s);
				free(s);
				break;
			}
			case 'A':
				fprintf(fp, "%d:%.4s", argp->u.aval->resSeq, argp->u.aval->aname);
				break;
			default:
				fprintf(fp, "\?\?\?\?");
				break;
		}
	}
	fprintf(fp, "\n");
	return 0;
}

static void
s_UpdateSelection(Molecule *mol, const IntGroup *ig, int is_insert)
{
	int n, old_natoms;
	IntGroup *sel, *orig_atoms, *ig3;
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
	
	IntGroupRelease(orig_atoms);
}

static int
s_MolActionAddAnAtom(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int *ip, n1;
	IntGroup *ig;
	n1 = MoleculeCreateAnAtom(mol, action->args[0].u.aval, action->args[1].u.ival);
	if ((ip = action->args[2].u.retval.ptr) != NULL)
		*ip = n1;
	if (n1 < 0)
		return -1;
	ig = IntGroupNewWithPoints(n1, 1, -1);
	s_UpdateSelection(mol, ig, 1);
	IntGroupRelease(ig);
	*actp = MolActionNew(gMolActionDeleteAnAtom, n1);
	return 0;
}

static int
s_MolActionMergeMolecule(Molecule *mol, MolAction *action, MolAction **actp)
{
	int n1, result, regOffset;
	Atom *ap;
	IntGroup *ig, *ig2;
	Molecule *mol2;
	Int nUndoActions;
	MolAction **undoActions;
	Int forUndo = 0;

	if (strcmp(action->name, gMolActionMergeMoleculeForUndo) == 0)
		forUndo = 1;

	ig = action->args[1].u.igval;
	mol2 = action->args[0].u.mval;
	if (ig != NULL) {
		ig2 = ig;
		IntGroupRetain(ig2);
	} else {
		ig2 = IntGroupNew();
		IntGroupAdd(ig2, mol->natoms, mol2->natoms);
	}

	if (forUndo == 0) {
		int minreg, maxreg;
		/*  Calculate the offset for the residue number  */
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
		regOffset = n1;
	} else regOffset = 0;

	nUndoActions = 0;
	undoActions = NULL;
	if ((result = MoleculeMerge(mol, mol2, ig, regOffset, &nUndoActions, &undoActions, forUndo)) != 0)
		return result;
	
	s_UpdateSelection(mol, ig, 1);
	
	/*  Register undo actions after registering unmerge action  */
	*actp = MolActionNew(gMolActionUnmergeMoleculeForUndo, ig2);
	IntGroupRelease(ig2);
	for (n1 = 0; n1 < nUndoActions; n1++) {
		MolActionCallback_registerUndo(mol, *actp);
		MolActionRelease(*actp);
		*actp = undoActions[n1];
	}
	/*  The last MolAction entry will be returned in *actp  */
	free(undoActions);
	
	/*  Sanity check (for debug)  */
	MoleculeCheckSanity(mol);
	
	return 0;
}

static int
s_MolActionDeleteAtoms(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int result, *ip, forUndo = 0;
	IntGroup *ig;
	Molecule *mol2;
	Int nUndoActions;
	MolAction **undoActions;

	if (strcmp(action->name, gMolActionDeleteAnAtom) == 0) {
		ig = IntGroupNew();
		if (ig == NULL || IntGroupAdd(ig, action->args[0].u.ival, 1) != 0)
			return -1;
	} else {
		ig = action->args[0].u.igval;
		IntGroupRetain(ig);
		if (strcmp(action->name, gMolActionUnmergeMoleculeForUndo) == 0)
			forUndo = 1;
	}
	
	/*  Unmerge molecule  */
	nUndoActions = 0;
	undoActions = NULL;
	if ((result = MoleculeUnmerge(mol, &mol2, ig, 0, &nUndoActions, &undoActions, forUndo)) != 0) {
		if (ip != NULL)
			free(ip);
		return result;
	}
	
	s_UpdateSelection(mol, ig, 0);
	
	if (mol2 == NULL)
		*actp = NULL;
	else {
		/*  Register undo actions created by MoleculeUnmerge()  */
		int i;
		for (i = 0; i < nUndoActions; i++) {
			MolActionCallback_registerUndo(mol, undoActions[i]);
			MolActionRelease(undoActions[i]);
		}
		free(undoActions);
		*actp = MolActionNew(gMolActionMergeMoleculeForUndo, mol2, ig);
		MoleculeRelease(mol2);
	}
	IntGroupRelease(ig);

	/*  Sanity check (for debug)  */
	MoleculeCheckSanity(mol);
	
	return 0;
}

static int
s_MolActionAddStructuralElements(Molecule *mol, MolAction *action, MolAction **actp, int type)
{
	Int *ip, n1, n2, result;
	IntGroup *ig;
	
	ip = (Int *)action->args[0].u.arval.ptr;

	if (type == 0 || type == 100) {  /*  bond  */
		Int na, nd;
		IntGroup *ig2;
		n1 = action->args[0].u.arval.nitems / 2;
		ig = action->args[1].u.igval;
		n2 = mol->nbonds;
		na = mol->nangles;
		nd = mol->ndihedrals;
		if ((result = MoleculeAddBonds(mol, n1, ip, ig, (type == 0))) <= 0)
			return result;
		if (ig == NULL)
			ig = IntGroupNewWithPoints(n2, mol->nbonds - n2, -1);
		else
			IntGroupRetain(ig);
		*actp = MolActionNew(gMolActionDeleteBonds, ig);
		IntGroupRelease(ig);
		/*  Register undo for creation of angle and dihedral  */
		if (mol->nangles > na) {
			MolActionCallback_registerUndo(mol, *actp);
			MolActionRelease(*actp);
			ig2 = IntGroupNewWithPoints(na, mol->nangles - na, -1);
			*actp = MolActionNew(gMolActionDeleteAngles, ig2);
			IntGroupRelease(ig2);
		}
		if (mol->ndihedrals > nd) {
			MolActionCallback_registerUndo(mol, *actp);
			MolActionRelease(*actp);
			ig2 = IntGroupNewWithPoints(na, mol->ndihedrals - nd, -1);
			*actp = MolActionNew(gMolActionDeleteDihedrals, ig2);
			IntGroupRelease(ig2);
		}
	} else if (type == 1) {  /*  angle  */
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
		*actp = MolActionNew(gMolActionDeleteAngles, ig);
		IntGroupRelease(ig);
	} else if (type == 2) {  /*  dihedral  */
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
		*actp = MolActionNew(gMolActionDeleteDihedrals, ig);
		IntGroupRelease(ig);		
	} else if (type == 3) {  /*  improper  */
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
		*actp = MolActionNew(gMolActionDeleteImpropers, ig);
		IntGroupRelease(ig);
	}
	
	/*  Sanity check (for debug)  */
	MoleculeCheckSanity(mol);
		
	return 0;
}

static int
s_MolActionDeleteStructuralElements(Molecule *mol, MolAction *action, MolAction **actp, int type)
{
	Int *ip, *ip2, n1, n2, result;
	IntGroup *ig, *ig2;
	MolAction *act2;
	if (type == 0) {  /*  bond  */
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig);
		if (n1 == 0)
			return 0;
		ip = (Int *)malloc(sizeof(Int) * n1 * 2);
		if ((n2 = MoleculeDeleteBonds(mol, ip, ig, &ip2, &ig2)) < 0) {
			free(ip);
			return n2;
		}
		if (n2 > 0) {
			/*  Register undo for restoring angles, dihedrals and impropers  */
			IntGroup *ig3, *ig4;
			Int *ip3 = ip2;
			/*  Angles  */
			ig3 = IntGroupNewFromIntGroup(ig2);
			ig4 = IntGroupNewWithPoints(ATOMS_MAX_NUMBER, ATOMS_MAX_NUMBER * 2, -1);
			IntGroupRemoveIntGroup(ig3, ig4);
			n2 = IntGroupGetCount(ig3);
			if (n2 > 0) {
				act2 = MolActionNew(gMolActionAddAngles, n2 * 3, ip3, ig3);
				MolActionCallback_registerUndo(mol, act2);
				MolActionRelease(act2);
				ip3 += n2 * 3;
			}
			IntGroupRelease(ig3);
			/*  Dihedrals  */
			ig3 = IntGroupNewFromIntGroup(ig2);
			IntGroupClear(ig4);
			IntGroupAdd(ig4, 0, ATOMS_MAX_NUMBER);
			IntGroupAdd(ig4, ATOMS_MAX_NUMBER * 2, ATOMS_MAX_NUMBER);
			IntGroupRemoveIntGroup(ig3, ig4);
			IntGroupOffset(ig3, -ATOMS_MAX_NUMBER);
			n2 = IntGroupGetCount(ig3);
			if (n2 > 0) {
				act2 = MolActionNew(gMolActionAddDihedrals, n2 * 4, ip3, ig3);
				MolActionCallback_registerUndo(mol, act2);
				MolActionRelease(act2);
				ip3 += n2 * 4;
			}
			IntGroupRelease(ig3);
			/*  Dihedrals  */
			ig3 = IntGroupNewFromIntGroup(ig2);
			IntGroupClear(ig4);
			IntGroupAdd(ig4, 0, ATOMS_MAX_NUMBER * 2);
			IntGroupRemoveIntGroup(ig3, ig4);
			IntGroupOffset(ig3, -ATOMS_MAX_NUMBER * 2);
			n2 = IntGroupGetCount(ig3);
			if (n2 > 0) {
				act2 = MolActionNew(gMolActionAddImpropers, n2 * 4, ip3, ig3);
				MolActionCallback_registerUndo(mol, act2);
				MolActionRelease(act2);
				ip3 += n2 * 4;
			}
			IntGroupRelease(ig3);
			IntGroupRelease(ig4);
			free(ip2);
		}
		*actp = MolActionNew(gMolActionAddBondsForUndo, n1 * 2, ip, ig);
		free(ip);
	} else if (type == 1) {  /*  angle  */
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig) * 3;
		ip = (Int *)malloc(sizeof(Int) * n1);
		if ((result = MoleculeDeleteAngles(mol, ip, ig)) < 0) {
			free(ip);
			return result;
		}
		*actp = MolActionNew(gMolActionAddAngles, n1, ip, ig);
		free(ip);
	} else if (type == 2) {  /*  dihedral  */
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig) * 4;
		ip = (Int *)malloc(sizeof(Int) * n1);
		if ((result = MoleculeDeleteDihedrals(mol, ip, ig)) < 0) {
			free(ip);
			return result;
		}
		*actp = MolActionNew(gMolActionAddDihedrals, n1, ip, ig);
		free(ip);
	} else if (type == 3) {  /*  improper  */
		ig = action->args[0].u.igval;
		n1 = IntGroupGetCount(ig) * 4;
		ip = (Int *)malloc(sizeof(Int) * n1);
		if ((result = MoleculeDeleteImpropers(mol, ip, ig)) < 0) {
			free(ip);
			return result;
		}
		*actp = MolActionNew(gMolActionAddImpropers, n1, ip, ig);
		free(ip);
	}

	/*  Sanity check (for debug)  */
	MoleculeCheckSanity(mol);
	
	return 0;
}

static int
s_MolActionTransformAtoms(Molecule *mol, MolAction *action, MolAction **actp, IntGroup **igp, int type)
{
	Vector *vp, v;
	if (type == 0) {  /*  translate  */
		vp = (Vector *)action->args[0].u.arval.ptr;
		if (vp == NULL)
			return -1;
		*igp = action->args[1].u.igval;
		MoleculeTranslate(mol, vp, *igp);
		VecScale(v, *vp, -1);
		*actp = MolActionNew(gMolActionTranslateAtoms, &v, *igp);
		(*actp)->frame = mol->cframe;
	} else if (type == 1) {  /*  rotate  */
		Double ang;
		Vector *vp2;
		vp = (Vector *)action->args[0].u.arval.ptr;
		ang = action->args[1].u.dval;
		vp2 = (Vector *)action->args[2].u.arval.ptr;
		*igp = action->args[3].u.igval;
		MoleculeRotate(mol, vp, ang, vp2, *igp);
		*actp = MolActionNew(gMolActionRotateAtoms, vp, -ang, vp2, *igp);
		(*actp)->frame = mol->cframe;
	} else if (type == 2) {  /*  general transform  */
		Transform *trp, tr_inv;
		int atomPositions = 0; /* Save atom positions for undo? */
		trp = (Transform *)action->args[0].u.arval.ptr;
		*igp = action->args[1].u.igval;
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
			int j, k, n1;
			n1 = IntGroupGetCount(*igp);
			vp = (Vector *)malloc(sizeof(Vector) * n1);
			if (vp == NULL) {
				*igp = NULL;
				return -1;
			}
			IntGroupIteratorInit(*igp, &iter);
			k = 0;
			while ((j = IntGroupIteratorNext(&iter)) >= 0) {
				vp[k++] = (ATOM_AT_INDEX(mol->atoms, j))->r;
			}
			*actp = MolActionNew(gMolActionSetAtomPositions, *igp, n1, vp);
			free(vp);
			IntGroupIteratorRelease(&iter);
		} else {
			*actp = MolActionNew(gMolActionTransformAtoms, &tr_inv, *igp);
		}
		MoleculeTransform(mol, *trp, *igp);
		(*actp)->frame = mol->cframe;		
	}
	if (*igp != NULL)
		IntGroupRetain(*igp);
	return 0;
}

static int
s_MolActionSetAtomGeometry(Molecule *mol, MolAction *action, MolAction **actp, IntGroup **igp, int type)
{
	IntGroupIterator iter;
	int j, k, n2;
	Vector *vp;
	Atom *ap;
	*igp = action->args[0].u.igval;
	n2 = IntGroupGetCount(*igp);
	vp = (Vector *)malloc(sizeof(Vector) * n2);
	if (vp == NULL) {
		*igp = NULL;
		return -1;
	}
	IntGroupIteratorInit(*igp, &iter);
	k = 0;
	while ((j = IntGroupIteratorNext(&iter)) >= 0) {
		ap = ATOM_AT_INDEX(mol->atoms, j);
		vp[k++] = (type == 0 ? ap->r : (type == 1 ? ap->v : ap->f));
	}
	*actp = MolActionNew(gMolActionSetAtomPositions, *igp, n2, vp);
	free(vp);
	vp = (Vector *)action->args[1].u.arval.ptr;
	IntGroupIteratorReset(&iter);
	k = 0;
	while ((j = IntGroupIteratorNext(&iter)) >= 0) {
		Vector w = vp[k++];
		ap = ATOM_AT_INDEX(mol->atoms, j);
		if (type == 0)
			ap->r = w;
		else if (type == 1)
			ap->v = w;
		else ap->f = w;
	}
	IntGroupIteratorRelease(&iter);
	(*actp)->frame = mol->cframe;
	if (*igp != NULL)
		IntGroupRetain(*igp);
	return 0;
}

static int
s_MolActionInsertFrames(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int n1, old_nframes, new_nframes, old_cframe;
	Vector *vp, *vp2;
	IntGroup *ig;
	MolAction *act2;

	ig = action->args[0].u.igval;
	vp = (Vector *)action->args[1].u.arval.ptr;
	vp2 = (Vector *)action->args[2].u.arval.ptr;
	old_cframe = mol->cframe;
	n1 = IntGroupGetCount(ig);
	if (n1 == 0)
		return 0;  /*  Do nothing  */
	if (vp != NULL && action->args[1].u.arval.nitems != n1 * mol->natoms)
		return -1;  /*  Internal inconsistency  */
	if (vp2 != NULL && action->args[2].u.arval.nitems != n1 * 4)
		return -1;  /*  Internal inconsistency  */
	old_nframes = MoleculeGetNumberOfFrames(mol);
	if (MoleculeInsertFrames(mol, ig, vp, vp2) < 0)
		return -1;  /*  Error  */
	
	/*  Undo action for restoring old cframe  */
	act2 = MolActionNew(gMolActionNone);
	act2->frame = old_cframe;
	MolActionCallback_registerUndo(mol, act2);
	MolActionRelease(act2);
	
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
	*actp = MolActionNew(gMolActionRemoveFrames, ig);
	(*actp)->frame = mol->cframe;
	return 0;
}

static int
s_MolActionRemoveFrames(Molecule *mol, MolAction *action, MolAction **actp)
{
	Vector *vp, *vp2;
	IntGroup *ig, *ig2;
	Int n1, n2, nframes, old_cframe;
	MolAction *act2;
	
	ig = ig2 = action->args[0].u.igval;
	old_cframe = mol->cframe;
	n1 = IntGroupGetCount(ig);
	if (n1 == 0)
		return 0;  /*  Do nothing  */
	nframes = MoleculeGetNumberOfFrames(mol);
	n2 = IntGroupGetEndPoint(ig, IntGroupGetIntervalCount(ig) - 1);  /*  Max point + 1  */
	if (n2 > nframes) {
		/*  Remove extra points  */
		ig2 = IntGroupNewFromIntGroup(ig);
		IntGroupRemove(ig2, nframes, n2 - nframes);
		n1 = IntGroupGetCount(ig2);
	}
	if (nframes == n1 && nframes >= 2) {
		/*  Remove all frames: keep the current frame  */
		if (ig2 == ig)
			ig2 = IntGroupNewFromIntGroup(ig);
		IntGroupRemove(ig2, mol->cframe, 1);
		n1--;
	}
	if (n1 == 0) {
		if (ig2 != ig)
			IntGroupRelease(ig2);
		return 0;  /*  Do nothing  */
	}
	vp = (Vector *)calloc(sizeof(Vector), n1 * mol->natoms);
	if (mol->cell != NULL && mol->frame_cells != NULL)
		vp2 = (Vector *)calloc(sizeof(Vector) * 4, n1);
	else vp2 = NULL;
	if (MoleculeRemoveFrames(mol, ig2, vp, vp2) < 0) {
		if (ig2 != ig)
			IntGroupRelease(ig2);
		return -1;  /*  Error  */
	}
	/*  Undo action for restoring old cframe  */
	act2 = MolActionNew(gMolActionNone);
	act2->frame = old_cframe;
	MolActionCallback_registerUndo(mol, act2);
	MolActionRelease(act2);
	
	*actp = MolActionNew(gMolActionInsertFrames, ig2, n1 * mol->natoms, vp, (vp2 != NULL ? n1 * 4 : 0), vp2);
	(*actp)->frame = mol->cframe;

	free(vp);
	if (vp2 != NULL)
		free(vp2);
	if (ig2 != ig)
		IntGroupRelease(ig2);
	return 0;
}

static int
s_MolActionSetSelection(Molecule *mol, MolAction *action, MolAction **actp)
{
	IntGroup *ig, *ig2;
	ig2 = MoleculeGetSelection(mol);
	if (ig2 != NULL)
		IntGroupRetain(ig2);  /*  To avoid releasing during MoleculeSetSelection() */
	ig = IntGroupNewFromIntGroup(action->args[0].u.igval);  /*  Duplicate so that selection change does not affect the original value */
	MoleculeSetSelection(mol, ig);
	IntGroupRelease(ig);
	*actp = MolActionNew(gMolActionSetSelection, ig2);
	if (ig2 != NULL)
		IntGroupRelease(ig2);
	return 0;
}

static int
s_MolActionRenumberAtoms(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int *ip, n1, result;
	Int *ip2 = (Int *)malloc(sizeof(Int) * (mol->natoms + 1));
	if (ip2 == NULL)
		return -1;
	ip = (Int *)action->args[0].u.arval.ptr;
	n1 = action->args[0].u.arval.nitems;
	result = MoleculeRenumberAtoms(mol, ip, ip2, n1);
	if (result != 0) {
		free(ip2);
		return result;
	}
	*actp = MolActionNew(gMolActionRenumberAtoms, mol->natoms, ip2);
	return 0;
}

static int
s_MolActionChangeResidueNumber(Molecule *mol, MolAction *action, MolAction **actp, int forUndo)
{
	IntGroup *ig;
	IntGroupIterator iter;
	Int i, n1, *ip, nresidues;

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
	*actp = MolActionNew(gMolActionChangeResidueNumberForUndo, ig, n1, ip, nresidues);
	return 0;
}

static int
s_MolActionOffsetResidueNumbers(Molecule *mol, MolAction *action, MolAction **actp)
{
	IntGroup *ig;
	Int n1, result;
	ig = action->args[0].u.igval;
	n1 = mol->nresidues;
	result = MoleculeOffsetResidueNumbers(mol, ig, action->args[1].u.ival, action->args[2].u.ival);
	if (result != 0)
		return result;  /*  The molecule is not modified  */
	*actp = MolActionNew(gMolActionOffsetResidueNumbers, ig, -(action->args[1].u.ival), n1);
	return 0;
}

static int
s_MolActionChangeResidueNames(Molecule *mol, MolAction *action, MolAction **actp)
{
	char *new_names, *old_names;
	Int *ip, i, argc;
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
	*actp = MolActionNew(gMolActionChangeResidueNames, argc, ip, argc * 4, old_names);
	free(old_names);
	return 0;
}

static int
s_MolActionChangeNumberOfResidues(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int n1, nresidues = mol->nresidues;
	n1 = action->args[0].u.ival;
	if (n1 < nresidues) {
		/*  The residue names will be lost, so undo must be registered to restore the names  */
		int argc = nresidues - n1;
		char *names = (char *)malloc(4 * argc + 1);
		Int *ip = (Int *)malloc(sizeof(Int) * (argc + 1));
		int i;
		MolAction *act2;
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
	*actp = MolActionNew(gMolActionChangeNumberOfResidues, nresidues);
	return 0;
}
	
static int
s_MolActionExpandBySymmetry(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int n1, *ip, count;
	Symop symop;
	IntGroup *ig = action->args[0].u.igval;
	count = IntGroupGetCount(ig);
	if (count == 0)
		return 0;  /*  No operation  */
	symop.dx = action->args[1].u.ival;
	symop.dy = action->args[2].u.ival;
	symop.dz = action->args[3].u.ival;
	symop.sym = action->args[4].u.ival;
	symop.alive = (symop.dx != 0 || symop.dy != 0 || symop.dz != 0 || symop.sym != 0);
	if (action->args[5].u.retval.ptr != NULL) {
		/*  Request the indices of the atoms  */
		ip = (Int *)calloc(sizeof(Int), count);
	} else ip = NULL;
	n1 = MoleculeAddExpandedAtoms(mol, symop, ig, ip, 0);
	if (n1 > 0) {
		ig = IntGroupNew();
		IntGroupAdd(ig, mol->natoms - n1, n1);
		(*actp) = MolActionNew(gMolActionUnmergeMolecule, ig);
		IntGroupRelease(ig);
	}
	if (ip != NULL) {
		if (n1 < 0) {
			free(ip);
			ip = NULL;
			count = 0;
		}
		*((Int **)(action->args[5].u.retval.ptr)) = ip;
		*(action->args[5].u.retval.nptr) = count;
	}
	return (n1 >= 0 ? 0 : n1);
}

static int
s_MolActionAmendBySymmetry(Molecule *mol, MolAction *action, MolAction **actp)
{
	IntGroup *ig1, *ig2;
	Vector *vp2;
	int n1;
	ig1 = action->args[0].u.igval;
	n1 = MoleculeAmendBySymmetry(mol, ig1, &ig2, &vp2);
	if (action->args[1].u.retval.ptr != NULL) {
		*((IntGroup **)(action->args[1].u.retval.ptr)) = ig2;
		IntGroupRetain(ig2);
	}
	if (n1 > 0) {
		*actp = MolActionNew(gMolActionSetAtomPositions, ig2, n1, vp2);
		free(vp2);
		IntGroupRelease(ig2);
	}
	mol->needsMDCopyCoordinates = 1;
	return 0;
}

static int
s_MolActionAddSymmetryOperation(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int n1;
	Transform *trp;
	Transform itr = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
	trp = (Transform *)action->args[0].u.arval.ptr;
	if (mol->nsyms == 0) {
		for (n1 = 0; n1 < 12; n1++) {
			if (fabs((*trp)[n1] - itr[n1]) > 1e-8)
				break;
		}
		if (n1 < 12) {
			MolAction *act2;
			if (AssignArray(&mol->syms, &mol->nsyms, sizeof(Transform), mol->nsyms, &itr) == 0)
				return -1;
			act2 = MolActionNew(gMolActionDeleteSymmetryOperation);
			MolActionCallback_registerUndo(mol, act2);
			MolActionRelease(act2);
		}
	}
	if (AssignArray(&mol->syms, &mol->nsyms, sizeof(Transform), mol->nsyms, trp) == 0)
		return -1;
	*actp = MolActionNew(gMolActionDeleteSymmetryOperation);
	return 0;
}

static int
s_MolActionDeleteSymmetryOperation(Molecule *mol, MolAction *action, MolAction **actp)
{
	if (mol->nsyms == 0)
		return -1;
	*actp = MolActionNew(gMolActionAddSymmetryOperation, &(mol->syms[mol->nsyms - 1]));
	mol->nsyms--;
	if (mol->nsyms == 0) {
		free(mol->syms);
		mol->syms = NULL;
	}
	return 0;
}

static int
s_MolActionSetCell(Molecule *mol, MolAction *action, MolAction **actp)
{
	double *dp, d[12];
	Vector vecs[4];
	Int oflags;
	Int convertCoord, n1, n2;
	if (mol->cell == NULL) {
		d[0] = 0.0;
		n1 = 0;
		oflags = 0;
	} else {
		for (n1 = 0; n1 < 6; n1++)
			d[n1] = mol->cell->cell[n1];
		if (mol->cell->has_sigma) {
			for (n1 = 6; n1 < 12; n1++)
				d[n1] = mol->cell->cellsigma[n1 - 6];
		}
		memmove(vecs, mol->cell->axes, sizeof(Vector) * 3);
		vecs[3] = mol->cell->origin;
		oflags = (mol->cell->flags[0] != 0) * 4 + (mol->cell->flags[1] != 0) * 2 + (mol->cell->flags[2] != 0);
	}
	convertCoord = action->args[1].u.ival;
	dp = action->args[0].u.arval.ptr;
	n2 = action->args[0].u.arval.nitems;
	if (n2 == 0)
		MoleculeSetCell(mol, 0, 0, 0, 0, 0, 0, convertCoord);
	else {
		MoleculeSetCell(mol, dp[0], dp[1], dp[2], dp[3], dp[4], dp[5], convertCoord);
		if (n2 == 12) {
			mol->cell->has_sigma = 1;
			for (n2 = 6; n2 < 12; n2++)
				mol->cell->cellsigma[n2 - 6] = dp[n2];
		} else mol->cell->has_sigma = 0;
	}
	if (n1 == 0)
		*actp = MolActionNew(gMolActionClearBox);
	else {
		*actp = MolActionNew(gMolActionSetBox, &vecs[0], &vecs[1], &vecs[2], &vecs[3], oflags, convertCoord);
		if (n1 > 6) {
			/*  Two undo actions are needed: first is for restore the cell vectors, and second is for restore the sigmas  */
			MolAction *act2;
			act2 = MolActionNew(gMolActionSetCell, n1, d, 0);
			MolActionCallback_registerUndo(mol, act2);
			MolActionRelease(act2);
		}
	}
	return 0;
}

static int
s_MolActionSetCellPeriodicity(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int flags, oflags;
	if (mol->cell == NULL)
		return 0;  /*  Do nothing  */
	oflags = (mol->cell->flags[0] != 0) * 4 + (mol->cell->flags[1] != 0) * 2 + (mol->cell->flags[2] != 0);
	flags = action->args[0].u.ival;
	mol->cell->flags[0] = ((flags & 4) != 0);
	mol->cell->flags[1] = ((flags & 2) != 0);
	mol->cell->flags[2] = ((flags & 1) != 0);
	*actp = MolActionNew(gMolActionSetCellPeriodicity, oflags);
	return 0;
}

static int
s_MolActionSetBox(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int n1, n2;
	if (action == NULL) {
		/*  Clear box  */
		if (mol->cell == NULL)
			return 0;  /*  Do nothing  */
		n1 = ((mol->cell->flags[0] != 0) * 4 + (mol->cell->flags[1] != 0) * 2 + (mol->cell->flags[2] != 0));
		*actp = MolActionNew(gMolActionSetBox, &(mol->cell->axes[0]), &(mol->cell->axes[1]), &(mol->cell->axes[2]), &(mol->cell->origin), n1, 0);
		MoleculeSetPeriodicBox(mol, NULL, NULL, NULL, NULL, NULL, 0);
	} else {
		/*  Set box  */
		Vector v[4];
		char flags[3];
		int convertCoordinates = action->args[5].u.ival;
		if (mol->cell == NULL)
			*actp = MolActionNew(gMolActionClearBox);
		else {
			n1 = ((mol->cell->flags[0] != 0) * 4 + (mol->cell->flags[1] != 0) * 2 + (mol->cell->flags[2] != 0));
			*actp = MolActionNew(gMolActionSetBox, &(mol->cell->axes[0]), &(mol->cell->axes[1]), &(mol->cell->axes[2]), &(mol->cell->origin), n1, convertCoordinates);
		}
		for (n1 = 0; n1 < 4; n1++)
			v[n1] = *((Vector *)(action->args[n1].u.arval.ptr));
		n2 = action->args[4].u.ival;
		if (n2 < 0) {
			/*  Keep existing flags; if not present, set all flags to 1.  */
			if (mol->cell == NULL)
				flags[0] = flags[1] = flags[2] = 1;
			else {
				flags[0] = mol->cell->flags[0];
				flags[1] = mol->cell->flags[1];
				flags[2] = mol->cell->flags[2];
			}
		} else {
			for (n1 = 0; n1 < 3; n1++)
				flags[n1] = ((n2 >> (2 - n1)) & 1);
		}
		MoleculeSetPeriodicBox(mol, &v[0], &v[1], &v[2], &v[3], flags, convertCoordinates);
	}
	return 0;
}

/*  This action is used for undoing "cell_flexibility = false"  */
/*
static int
s_MolActionSetBoxForFrames(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int i, n1, n2;
	Vector *vp1, *vp2;
	n2 = MoleculeGetNumberOfFrames(mol);
	if (n2 == 0 || mol->cell == NULL)
		return 0;  //  Do nothing 
	n1 = action->args[0].u.arval.nitems / 4;
	vp1 = (Vector *)(action->args[0].u.arval.ptr);
	if (mol->nframe_cells < n2) {
		//  Expand the array before processing 
		i = mol->nframe_cells * 4;
		AssignArray(&(mol->frame_cells), &(mol->nframe_cells), sizeof(Vector) * 4, n2 - 1, NULL);
		while (i < n2 * 4) {
			//  Copy the current cell 
			mol->frame_cells[i++] = mol->cell->axes[0];
			mol->frame_cells[i++] = mol->cell->axes[1];
			mol->frame_cells[i++] = mol->cell->axes[2];
			mol->frame_cells[i++] = mol->cell->origin;
		}
	}
	
	vp2 = (Vector *)malloc(sizeof(Vector) * n2 * 4);
	memmove(vp2, mol->frame_cells, sizeof(Vector) * n2 * 4);
	memmove(mol->frame_cells, vp1, sizeof(Vector) * 4 * (n1 < n2 ? n1 : n2));
	*actp = MolActionNew(gMolActionSetBoxForFrames, n2 * 4, vp2);
	free(vp2);

	//  Set the current cell (no change on the periodic flags)
	vp2 = mol->frame_cells + mol->cframe * 4;
	MoleculeSetPeriodicBox(mol, vp2, vp2 + 1, vp2 + 2, vp2 + 3, mol->cell->flags, 0);
	
	return 0;
} */

/*
static int
s_MolActionSetCellFlexibility(Molecule *mol, MolAction *action, MolAction **actp)
{
	Int n1;
	n1 = action->args[0].u.ival;
	if ((n1 != 0) == (mol->useFlexibleCell != 0))
		return 0;  //  Do nothing
	mol->useFlexibleCell = (n1 != 0);
	if (n1 == 0) {
		//  Clear the existing cells, and register undo
		if (mol->nframe_cells > 0) {
			MolAction *act2 = MolActionNew(gMolActionSetBoxForFrames, mol->nframe_cells * 4, mol->frame_cells);
			MolActionSetFrame(act2, mol->cframe);
			MolActionCallback_registerUndo(mol, act2);
			MolActionRelease(act2);
		}
		free(mol->frame_cells);
		mol->frame_cells = NULL;
		mol->nframe_cells = 0;
	} else {
		//  Allocate cells for all frames and copy the current cell
		Int i, nframes = MoleculeGetNumberOfFrames(mol);
		if (nframes != 0 && mol->cell != NULL) {
			if (mol->nframe_cells < nframes) {
				//  Expand the array
				AssignArray(&(mol->frame_cells), &(mol->nframe_cells), sizeof(Vector) * 4, nframes - 1, NULL);
			}
			//  Copy the current cell 
			//  (No undo action is registered; actually, the frame_cells array should be empty) 
			for (i = 0; i < nframes; i++) {
				mol->frame_cells[i * 4] = mol->cell->axes[0];
				mol->frame_cells[i * 4 + 1] = mol->cell->axes[1];
				mol->frame_cells[i * 4 + 2] = mol->cell->axes[2];
				mol->frame_cells[i * 4 + 3] = mol->cell->origin;
			}
		}
	}
	*actp = MolActionNew(gMolActionSetCellFlexibility, (n1 == 0));
	return 0;
}
*/

static int
s_MolActionAddParameters(Molecule *mol, MolAction *action, MolAction **actp)
{
	UnionPar *up;
	IntGroup *ig;
	Int parType;
	if (mol->par == NULL)
		return -1;
	parType = action->args[0].u.ival;
	ig = action->args[1].u.igval;
	up = action->args[2].u.arval.ptr;
	ParameterInsert(mol->par, parType, up, ig);
	*actp = MolActionNew(gMolActionDeleteParameters, parType, ig);
	mol->parameterTableSelectionNeedsClear = 1;
	return 0;
}

static int
s_MolActionDeleteParameters(Molecule *mol, MolAction *action, MolAction **actp)
{
	UnionPar *up;
	IntGroup *ig;
	Int parType, n1;
	if (mol->par == NULL)
		return -1;
	parType = action->args[0].u.ival;
	ig = action->args[1].u.igval;
	n1 = IntGroupGetCount(ig);
	up = (UnionPar *)calloc(sizeof(UnionPar), n1);
	ParameterDelete(mol->par, parType, up, ig);
	*actp = MolActionNew(gMolActionAddParameters, parType, ig, n1, up);
	free(up);
	mol->parameterTableSelectionNeedsClear = 1;
	return 0;
}

int
MolActionPerform(Molecule *mol, MolAction *action)
{
	int result;
	IntGroup *ig = NULL;
	MolAction *act2 = NULL;
	int needsSymmetryAmendment = 0;
	int needsRebuildMDArena = 0;

	if (action == NULL)
		return -1;
	
	{  /*  For debug  */
		int status;
		VALUE val = rb_eval_string_protect("$log_molby_actions ? true : false", &status);
		if (status == 0 && val != Qfalse) {
			s_MolActionLog(mol, action, stderr);
		}
	}
	
	if (action->frame >= 0 && mol->cframe != action->frame)
		MoleculeSelectFrame(mol, action->frame, 1);

	/*  Ruby script execution  */
	if (strncmp(action->name, kMolActionPerformScript, strlen(kMolActionPerformScript)) == 0) {
		return s_MolActionPerformRubyScript(mol, action);
	}
	
	if (mol == NULL)
		return -1;
	
	if (strcmp(action->name, gMolActionAddAnAtom) == 0) {
		if ((result = s_MolActionAddAnAtom(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionMergeMolecule) == 0 || strcmp(action->name, gMolActionMergeMoleculeForUndo) == 0) {
		if ((result = s_MolActionMergeMolecule(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteAnAtom) == 0 || strcmp(action->name, gMolActionUnmergeMolecule) == 0 || strcmp(action->name, gMolActionUnmergeMoleculeForUndo) == 0) {
		if ((result = s_MolActionDeleteAtoms(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddBonds) == 0) {
		if ((result = s_MolActionAddStructuralElements(mol, action, &act2, 0)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddBondsForUndo) == 0) {
		if ((result = s_MolActionAddStructuralElements(mol, action, &act2, 100)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteBonds) == 0) {
		if ((result = s_MolActionDeleteStructuralElements(mol, action, &act2, 0)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddAngles) == 0) {
		if ((result = s_MolActionAddStructuralElements(mol, action, &act2, 1)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteAngles) == 0) {
		if ((result = s_MolActionDeleteStructuralElements(mol, action, &act2, 1)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddDihedrals) == 0) {
		if ((result = s_MolActionAddStructuralElements(mol, action, &act2, 2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteDihedrals) == 0) {
		if ((result = s_MolActionDeleteStructuralElements(mol, action, &act2, 2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionAddImpropers) == 0) {
		if ((result = s_MolActionAddStructuralElements(mol, action, &act2, 3)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteImpropers) == 0) {
		if ((result = s_MolActionDeleteStructuralElements(mol, action, &act2, 3)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionTranslateAtoms) == 0) {
		if ((result = s_MolActionTransformAtoms(mol, action, &act2, &ig, 0)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionRotateAtoms) == 0) {
		if ((result = s_MolActionTransformAtoms(mol, action, &act2, &ig, 1)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionTransformAtoms) == 0) {
		if ((result = s_MolActionTransformAtoms(mol, action, &act2, &ig, 2)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetAtomPositions) == 0) {
		if ((result = s_MolActionSetAtomGeometry(mol, action, &act2, &ig, 0)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetAtomVelocities) == 0) {
		if ((result = s_MolActionSetAtomGeometry(mol, action, &act2, &ig, 1)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetAtomForces) == 0) {
		if ((result = s_MolActionSetAtomGeometry(mol, action, &act2, &ig, 2)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionInsertFrames) == 0) {
		if ((result = s_MolActionInsertFrames(mol, action, &act2)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionRemoveFrames) == 0) {
		if ((result = s_MolActionRemoveFrames(mol, action, &act2)) != 0)
			return result;
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetSelection) == 0) {
		if ((result = s_MolActionSetSelection(mol, action, &act2)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionRenumberAtoms) == 0) {
		if ((result = s_MolActionRenumberAtoms(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionChangeResidueNumber) == 0) {
		if ((result = s_MolActionChangeResidueNumber(mol, action, &act2, 0)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionChangeResidueNumberForUndo) == 0) {
		if ((result = s_MolActionChangeResidueNumber(mol, action, &act2, 1)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionOffsetResidueNumbers) == 0) {
		if ((result = s_MolActionOffsetResidueNumbers(mol, action, &act2)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionChangeResidueNames) == 0) {
		if ((result = s_MolActionChangeResidueNames(mol, action, &act2)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionChangeNumberOfResidues) == 0) {
		if ((result = s_MolActionChangeNumberOfResidues(mol, action, &act2)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionAmendBySymmetry) == 0) {
		if ((result = s_MolActionAmendBySymmetry(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionExpandBySymmetry) == 0) {
		if ((result = s_MolActionExpandBySymmetry(mol, action, &act2)) != 0)
			return result;
	} else if (strcmp(action->name, gMolActionAddSymmetryOperation) == 0) {
		if ((result = s_MolActionAddSymmetryOperation(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteSymmetryOperation) == 0) {
		if ((result = s_MolActionDeleteSymmetryOperation(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionSetCell) == 0) {
		if ((result = s_MolActionSetCell(mol, action, &act2)) != 0)
			return result;
		if (mol->arena != NULL)
			md_set_cell(mol->arena);
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetCellPeriodicity) == 0) {
		if ((result = s_MolActionSetCellPeriodicity(mol, action, &act2)) != 0)
			return result;
		if (mol->arena != NULL)
			md_set_cell(mol->arena);
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionSetBox) == 0) {
		if ((result = s_MolActionSetBox(mol, action, &act2)) != 0)
			return result;
		if (mol->arena != NULL)
			md_set_cell(mol->arena);
		needsSymmetryAmendment = 1;
	} else if (strcmp(action->name, gMolActionClearBox) == 0) {
		if ((result = s_MolActionSetBox(mol, NULL, &act2)) != 0)
			return result;
		if (mol->arena != NULL)
			md_set_cell(mol->arena);
		needsSymmetryAmendment = 1;
/*	} else if (strcmp(action->name, gMolActionSetBoxForFrames) == 0) {
		if ((result = s_MolActionSetBoxForFrames(mol, action, &act2)) != 0)
			return result; */
/*	} else if (strcmp(action->name, gMolActionSetCellFlexibility) == 0) {
		if ((result = s_MolActionSetCellFlexibility(mol, action, &act2)) != 0)
			return result; */
	} else if (strcmp(action->name, gMolActionAddParameters) == 0) {
		if ((result = s_MolActionAddParameters(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionDeleteParameters) == 0) {
		if ((result = s_MolActionDeleteParameters(mol, action, &act2)) != 0)
			return result;
		needsRebuildMDArena = 1;
	} else if (strcmp(action->name, gMolActionNone) == 0) {
		/*  Do nothing  */
	} else {
		fprintf(stderr, "Internal error: unknown action name %s\n", action->name);
		return -1;
	}
	if (act2 != NULL) {
		MolActionCallback_registerUndo(mol, act2);
		MolActionRelease(act2);
	}

	if (needsSymmetryAmendment) {
		MolAction *act3;
		act3 = MolActionNew(gMolActionAmendBySymmetry, ig, NULL);
		if (ig != NULL)
			IntGroupRelease(ig);
		act2 = NULL;
		result = s_MolActionAmendBySymmetry(mol, act3, &act2);
		MolActionRelease(act3);
		if (result != 0)
			return result;
		if (act2 != NULL) {
			MolActionCallback_registerUndo(mol, act2);
			MolActionRelease(act2);
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
