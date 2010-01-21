/*
 *  ruby_md.c
 *  Ruby binding
 *
 *  Created by Toshi Nagata on 09/08/11.
 *  Copyright 2007-2009 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
*/

#include "Molby.h"
#include "../MD/MDCore.h"
#include "../MD/MDSurface.h"
#include "../MD/MDPressure.h"

#include "env.h"  /*  For ruby_frame  */

#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#pragma mark ====== Global Values ======

VALUE rb_cMDArena;

#pragma mark ====== MDArena class ======

struct MDArena *
MDArenaFromValue(VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cMDArena)) {
		MDArena *arena;
		Data_Get_Struct(val, MDArena, arena);
		return arena;
	} else return NULL;
}

VALUE
ValueFromMDArena(struct MDArena *arena)
{
	if (arena != NULL)
		md_arena_retain(arena);
	return Data_Wrap_Struct(rb_cMDArena, 0, (void (*)(void *))md_arena_release, arena);
}

static void
s_MDArena_panic_func(MDArena *arena, const char *msg)
{
	rb_raise(rb_eMolbyError, "MD Error: %s", msg);
}

static VALUE
s_MDArena_Run_or_minimize(VALUE self, VALUE arg, int minimize)
{
	MDArena *arena;
	int nsteps, retval, start_step;
	IntGroup *ig;
	Data_Get_Struct(self, MDArena, arena);
	nsteps = NUM2INT(rb_Integer(arg));
	if (arena->is_running)
		rb_raise(rb_eMolbyError, "the simulation is already running. You cannot do simulation recursively.");
	if (nsteps < 0)
		rb_raise(rb_eMolbyError, "the number of steps should be non-negative integer");
	if (!arena->is_initialized || arena->xmol->needsMDRebuild) {
		const char *msg = md_prepare(arena, 0);
		if (msg != NULL)
			rb_raise(rb_eMolbyError, "cannot initialize for MD: %s", msg);
	}
	arena->end_step = arena->start_step + nsteps;
	md_log(arena, "%s for %d steps\n", (minimize ? "Minimizing" : "Running"), nsteps);
	start_step = arena->start_step;

	/*  Create a new frame with current coordinates  */
/*	if (nsteps > 0 && MoleculeGetNumberOfFrames(arena->xmol) == 0) {
		ig = IntGroupNewWithPoints(0, 1, -1);
		MolActionCreateAndPerform(arena->xmol, gMolActionInsertFrames, ig, 0, NULL);
		IntGroupRelease(ig);
	} */

	/*  Run simulation  */
/*	arena->md_panic_func = s_MDArena_panic_func; */
	retval = md_main(arena, minimize);
	if (retval != 0)
		rb_raise(rb_eMolbyError, "MD Error: %s", arena->errmsg);

/*	arena->md_panic_func = NULL; */
	
	if (arena->step > start_step) {
		/*  Create a new frame and copy new coordinates  */
		ig = IntGroupNewWithPoints(MoleculeGetNumberOfFrames(arena->xmol), 1, -1);
		MolActionCreateAndPerform(arena->xmol, gMolActionInsertFrames, ig, 0, NULL);
		IntGroupRelease(ig);
		md_copy_coordinates_from_internal(arena);
	}
	
/*	if (retval != 0)
		rb_raise(rb_eMolbyError, "Calculation aborted with status %d", retval); */
	
	if (minimize && arena->minimize_complete && rb_block_given_p())
		rb_yield(self);

	return self;
}	

/*
 *  call-seq:
 *     mdarena.run(n)       -> self
 *
 *  Run the simulation for n steps. 
 */
static VALUE
s_MDArena_Run(VALUE self, VALUE arg)
{
	return s_MDArena_Run_or_minimize(self, arg, 0);
}

/*
 *  call-seq:
 *     mdarena.minimize(n)       -> self
 *     mdarena.minimize(n) { ... } -> self
 *
 *  Minimize the energy for n steps. If a block is given, it is executed when minimization is complete.
 */
static VALUE
s_MDArena_Minimize(VALUE self, VALUE arg)
{
	return s_MDArena_Run_or_minimize(self, arg, 1);
}

/*
 *  call-seq:
 *     mdarena.prepare(check_only = false)         -> self or nil
 *
 *  Prepare for the MD calculation; refresh the internal information even if initialization is
 *  already done.
 *  If check_only is true, then only the parameter checking is done.
 *  Returns self when initialization is complete; returns nil if some MM parameters are missing;
 *  otherwise, an exception is raised.
 */
static VALUE
s_MDArena_Prepare(int argc, VALUE *argv, VALUE self)
{
	MDArena *arena;
	Molecule *mol;
	const char *msg;
	Int nangles, *angles, ndihedrals, *dihedrals, nimpropers, *impropers;
	Int missing = 0;
	Int check_only = 0;
	IntGroup *ig1, *ig2, *ig3;
	VALUE fval;

	Data_Get_Struct(self, MDArena, arena);
	rb_scan_args(argc, argv, "01", &fval);
	if (RTEST(fval))
		check_only = 1;

	arena->is_initialized = 0;
	mol = arena->xmol;
	
	/*  Rebuild the tables  */
	ig1 = ig2 = ig3 = NULL;
	nangles = MoleculeFindMissingAngles(mol, &angles);
	ndihedrals = MoleculeFindMissingDihedrals(mol, &dihedrals);
	nimpropers = MoleculeFindMissingImpropers(mol, &impropers);
	if (nangles > 0) {
		ig1 = IntGroupNewWithPoints(mol->nangles, nangles, -1);
		MolActionCreateAndPerform(mol, gMolActionAddAngles, nangles * 3, angles, ig1);
		free(angles);
		IntGroupRelease(ig1);
	}
	if (ndihedrals > 0) {
		ig2 = IntGroupNewWithPoints(mol->ndihedrals, ndihedrals, -1);
		MolActionCreateAndPerform(mol, gMolActionAddDihedrals, ndihedrals * 4, dihedrals, ig2);
		free(dihedrals);
		IntGroupRelease(ig2);
	}
	if (nimpropers > 0) {
		ig3 = IntGroupNewWithPoints(mol->nimpropers, nimpropers, -1);
		MolActionCreateAndPerform(mol, gMolActionAddImpropers, nimpropers * 4, impropers, ig3);
		free(impropers);
		IntGroupRelease(ig3);
	}

	/*  Prepare parameters and internal information  */
	msg = md_prepare(arena, check_only);

	/*  Some parameters are missing?  */
	if (msg != NULL) {
		if (strstr(msg, "parameter") != NULL && strstr(msg, "missing") != NULL)
			missing = 1;
		else
			rb_raise(rb_eMolbyError, "cannot initialize for MD: %s", msg);
	}

	/*  The local parameter list is updated  */
	{
		Int parType, idx;
		if (mol->par == NULL)
			mol->par = ParameterNew();
		for (parType = kFirstParType; parType <= kLastParType; parType++) {
			/*  Delete global and undefined parameters  */
			UnionPar *up, *upbuf;
			Int nparams, count;
			ig1 = IntGroupNew();
			for (idx = 0; (up = ParameterGetUnionParFromTypeAndIndex(mol->par, parType, idx)) != NULL; idx++) {
				if (up->bond.src != 0)
					IntGroupAdd(ig1, idx, 1);
			}
			if (IntGroupGetCount(ig1) > 0)
				MolActionCreateAndPerform(mol, gMolActionDeleteParameters, parType, ig1);
			IntGroupRelease(ig1);
			/*  Copy global and undefined parameters from arena and insert to mol->par  */
			nparams = ParameterGetCountForType(arena->par, parType);
			if (nparams == 0)
				continue;
			upbuf = (UnionPar *)calloc(sizeof(UnionPar), nparams);
			ig1 = IntGroupNew();
			ig2 = IntGroupNew();
			for (idx = 0; (up = ParameterGetUnionParFromTypeAndIndex(arena->par, parType, idx)) != NULL; idx++) {
				if (up->bond.src > 0)
					IntGroupAdd(ig1, idx, 1); /* Global parameter */
				else if (up->bond.src < 0)
					IntGroupAdd(ig2, idx, 1); /* Undefined parameter */
			}
			if ((count = IntGroupGetCount(ig1)) > 0) {
				/*  Insert global parameters (at the top)  */
				ParameterCopy(arena->par, parType, upbuf, ig1);
				ig3 = IntGroupNewWithPoints(0, count, -1);
				MolActionCreateAndPerform(mol, gMolActionAddParameters, parType, ig3, count, upbuf);
				IntGroupRelease(ig3);
			}
			if ((count = IntGroupGetCount(ig2)) > 0) {
				/*  Insert undefined parameters (at the bottom)  */
				ParameterCopy(arena->par, parType, upbuf, ig2);
				idx = ParameterGetCountForType(mol->par, parType);
				ig3 = IntGroupNewWithPoints(idx, count, -1);
				MolActionCreateAndPerform(mol, gMolActionAddParameters, parType, ig3, count, upbuf);
				IntGroupRelease(ig3);
			}
			IntGroupRelease(ig2);
			IntGroupRelease(ig1);
			free(upbuf);
		}
	}

	if (missing)
		return Qnil;
	else return self;
}

/*
 *  call-seq:
 *     arena.energies -> [total, bond, angle, dihedral, improper, vdw, electrostatic, auxiliary, surface, kinetic, net]
 *
 *  Get the current energies.
 */
static VALUE
s_MDArena_Energies(VALUE self)
{
	MDArena *arena;
	VALUE val;
	int i;
	Data_Get_Struct(self, MDArena, arena);
	val = rb_ary_new();
	rb_ary_push(val, rb_float_new(arena->total_energy * INTERNAL2KCAL));
	for (i = 0; i < kEndIndex; i++)
		rb_ary_push(val, rb_float_new(arena->energies[i] * INTERNAL2KCAL));
	rb_ary_push(val, rb_float_new((arena->energies[kKineticIndex] + arena->total_energy) * INTERNAL2KCAL));
	return val;
}

static VALUE s_LogFileSym, s_CoordFileSym, s_VelFileSym, s_ForceFileSym, 
s_DebugFileSym, s_DebugOutputLevelSym, s_StepSym, s_CoordOutputFreqSym, 
s_EnergyOutputFreqSym, s_CoordFrameSym, s_TimestepSym, s_CutoffSym, 
s_ElectroCutoffSym, s_PairlistDistanceSym, s_TemperatureSym, s_AndersenFreqSym, 
s_AndersenCouplingSym, s_RandomSeedSym, s_DielectricSym, s_GradientConvergenceSym, 
s_CoordinateConvergenceSym, s_UseXplorShiftSym, s_Scale14VdwSym, s_Scale14ElectSym, 
s_RelocateCenterSym, s_SurfaceProbeRadiusSym, s_SurfaceTensionSym, s_SurfacePotentialFreqSym,
s_UseGraphiteSym;

struct s_MDArenaAttrDef {
	char *name;
	VALUE *symref;  /*  Address of s_LogFileSym etc. */
	ID id;			/*  Ruby ID of the symbol; will be set within Init_MolbyMDTypes()  */
	ID sid;         /*  Ruby ID of the symbol plus '='; will be set within Init_MolbyMDTypes()  */
	char type;      /*  s: string (const char *), i: Int, f: Double. Uppercase: read-only.  */
	int  offset;    /*  Offset in the MDArena structure.  */
};
static struct s_MDArenaAttrDef s_MDArenaAttrDefTable[] = {
	{"log_file",          &s_LogFileSym,          0, 0, 's', offsetof(MDArena, log_result_name)},
	{"coord_file",        &s_CoordFileSym,        0, 0, 's', offsetof(MDArena, coord_result_name)},
	{"vel_file",          &s_VelFileSym,          0, 0, 's', offsetof(MDArena, vel_result_name)}, 
	{"force_file",        &s_ForceFileSym,        0, 0, 's', offsetof(MDArena, force_result_name)},
	{"debug_file",        &s_DebugFileSym,        0, 0, 's', offsetof(MDArena, debug_result_name)},
	{"debug_output_level", &s_DebugOutputLevelSym, 0, 0, 'i', offsetof(MDArena, debug_output_level)},
	{"step",              &s_StepSym,             0, 0, 'I', offsetof(MDArena, step)},
	{"coord_output_freq", &s_CoordOutputFreqSym,  0, 0, 'i', offsetof(MDArena, coord_output_freq)},
	{"energy_output_freq", &s_EnergyOutputFreqSym, 0, 0, 'i', offsetof(MDArena, energy_output_freq)},
	{"coord_frame",       &s_CoordFrameSym,       0, 0, 'I', offsetof(MDArena, coord_result_frame)},
	{"timestep",          &s_TimestepSym,         0, 0, 'f', offsetof(MDArena, timestep)},
	{"cutoff",            &s_CutoffSym,           0, 0, 'f', offsetof(MDArena, cutoff)},
	{"electro_cutoff",    &s_ElectroCutoffSym,    0, 0, 'f', offsetof(MDArena, electro_cutoff)},
	{"pairlist_distance", &s_PairlistDistanceSym, 0, 0, 'f', offsetof(MDArena, pairlist_distance)},
	{"temperature",       &s_TemperatureSym,      0, 0, 'f', offsetof(MDArena, temperature)},
	{"andersen_freq",     &s_AndersenFreqSym,     0, 0, 'i', offsetof(MDArena, andersen_thermo_freq)},
	{"andersen_coupling", &s_AndersenCouplingSym, 0, 0, 'f', offsetof(MDArena, andersen_thermo_coupling)},
	{"random_seed",       &s_RandomSeedSym,       0, 0, 'i', offsetof(MDArena, random_seed)},
	{"dielectric",        &s_DielectricSym,       0, 0, 'f', offsetof(MDArena, dielectric)},
	{"gradient_convergence", &s_GradientConvergenceSym, 0, 0, 'f', offsetof(MDArena, gradient_convergence)},
	{"coordinate_convergence", &s_CoordinateConvergenceSym, 0, 0, 'f', offsetof(MDArena, coordinate_convergence)},
	{"use_xplor_shift",   &s_UseXplorShiftSym,    0, 0, 'i', offsetof(MDArena, use_xplor_shift)},
	{"scale14_vdw",       &s_Scale14VdwSym,       0, 0, 'f', offsetof(MDArena, scale14_vdw)},
	{"scale14_elect",     &s_Scale14ElectSym,     0, 0, 'f', offsetof(MDArena, scale14_elect)},
	{"relocate_center",   &s_RelocateCenterSym,   0, 0, 'i', offsetof(MDArena, relocate_center)},
	{"surface_probe_radius", &s_SurfaceProbeRadiusSym, 0, 0, 'f', offsetof(MDArena, probe_radius)},
	{"surface_tension",   &s_SurfaceTensionSym,   0, 0, 'f', offsetof(MDArena, surface_tension)},
	{"surface_potential_freq", &s_SurfacePotentialFreqSym, 0, 0, 'i', offsetof(MDArena, surface_potential_freq)},
	{"use_graphite",      &s_UseGraphiteSym,      0, 0, 'i', offsetof(MDArena, use_graphite)},
	{NULL} /* Sentinel */
};

static VALUE s_PresFreqSym, s_PresCouplingSym, s_PresSym, s_PresCellFlexSym, s_PresFluctCellOriginSym, s_PresFluctCellOrientSym;
static struct s_MDArenaAttrDef s_MDPressureAttrDefTable[] = {
	{"pressure_freq",     &s_PresFreqSym,         0, 0, 'i', offsetof(MDPressureArena, freq)},
	{"pressure_coupling", &s_PresCouplingSym,     0, 0, 'f', offsetof(MDPressureArena, coupling)},
	{"pressure",          &s_PresSym,             0, 0, 'X', offsetof(MDPressureArena, apply)},
	{"pressure_cell_flexibility", &s_PresCellFlexSym, 0, 0, 'Y', offsetof(MDPressureArena, cell_flexibility)},
	{"pressure_fluctuate_cell_origin", &s_PresFluctCellOriginSym, 0, 0, 'f', offsetof(MDPressureArena, fluctuate_cell_origin)},
	{"pressure_fluctuate_cell_orientation", &s_PresFluctCellOrientSym, 0, 0, 'f', offsetof(MDPressureArena, fluctuate_cell_orientation)},
	{NULL} /* Sentinel */
};

/*
 *  call-seq:
 *     arena[attr]
 *
 *  Get the attribute value.
 */
static VALUE
s_MDArena_Get(VALUE self, VALUE attr)
{
	MDArena *arena;
	int i;
	struct s_MDArenaAttrDef *dp;
	ID aid = rb_to_id(attr);
	Data_Get_Struct(self, MDArena, arena);
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		if (dp->id == aid) {
			char *p = (char *)arena + dp->offset;
			switch (dp->type) {
				case 's':
				case 'S': {
					const char *cp = *((const char **)p);
					if (cp == NULL)
						return Qnil;
					else
						return Ruby_NewFileStringValue(cp);
				}
				case 'i':
				case 'I':
					return INT2NUM(*((Int *)p));
				case 'f':
				case 'F':
					return rb_float_new(*((Double *)p));
				default:
					rb_raise(rb_eMolbyError, "Internal inconsistency: unknown type field");
			}
		}
	}
	for (i = 0, dp = s_MDPressureAttrDefTable; dp->name != NULL; i++, dp++) {
		if (dp->id == aid) {
			char *pp;
			MDPressureArena *pres = arena->pressure;
			if (pres == NULL)
				return Qnil;
			pp = (char *)pres + dp->offset;
			switch (dp->type) {
				case 'i':
				case 'I':
					return INT2NUM(*((Int *)pp));
				case 'f':
				case 'F':
					return rb_float_new(*((Double *)pp));
				case 'X':
					/*  Isotropic pressure only  */
					return rb_ary_new3(3, rb_float_new(pres->apply[0]), rb_float_new(pres->apply[4]), rb_float_new(pres->apply[8]));
				case 'Y': {
					VALUE aval = rb_ary_new();
					int j;
					for (j = 0; j < 8; j++)
						rb_ary_push(aval, rb_float_new(pres->cell_flexibility[j]));
					return aval;
				}
				default:
					rb_raise(rb_eMolbyError, "Internal inconsistency: unknown type field");
			}
		}
	}
	rb_raise(rb_eMolbyError, "unknown attribute name (%s)", rb_id2name(aid));
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     arena.(attrname)
 *
 *  Get the attribute value. The name of the attribute is taken from ruby_frame->last_func.
 */
static VALUE
s_MDArena_GetAttr(VALUE self)
{
	return s_MDArena_Get(self, ID2SYM(ruby_frame->last_func));
}

/*
 *  call-seq:
 *     arena[attr]=
 *
 *  Set the attribute value.
 */
static VALUE
s_MDArena_Set(VALUE self, VALUE attr, VALUE val)
{
	MDArena *arena;
	int i, j;
	struct s_MDArenaAttrDef *dp;
	ID aid = rb_to_id(attr);
	attr = ID2SYM(aid);  /*  May be used later  */
	Data_Get_Struct(self, MDArena, arena);
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		if (dp->id == aid) {
			char *p = (char *)arena + dp->offset;
			switch (dp->type) {
				case 's': {
					const char *cp = (val == Qnil ? NULL : (const char *)strdup(FileStringValuePtr(val)));
					const char **cpp = (const char **)p;
					FILE **fpp;
					if (*cpp == cp || (*cpp != NULL && cp != NULL && strcmp(*cpp, cp) == 0))
						return val;  /*  No need to change  */
					if (*cpp != NULL)
						free((void *)*cpp);
					if (cp != NULL && cp[0] == 0) {
						free((void *)cp);
						cp = NULL;
					}
					/*  Close the corresponding FILE if necessary  */
					if (attr == s_LogFileSym)
						fpp = &(arena->log_result);
					else if (attr == s_CoordFileSym)
						fpp = &(arena->coord_result);
					else if (attr == s_VelFileSym)
						fpp = &(arena->vel_result);
					else if (attr == s_ForceFileSym)
						fpp = &(arena->force_result);
					else if (attr == s_DebugFileSym)
						fpp = &(arena->debug_result);
					else fpp = NULL;
					if (fpp != NULL && *fpp != NULL) {
						fclose(*fpp);
						*fpp = NULL;
					}
					*cpp = cp;
					return val;
				}
				case 'i':
					*((Int *)p) = NUM2INT(rb_Integer(val));
					return val;
				case 'f':
					*((Double *)p) = NUM2DBL(rb_Float(val));
					return val;
				case 'S': case 'I': case 'F':
					rb_raise(rb_eMolbyError, "The attribute is read-only");
				default:
					rb_raise(rb_eMolbyError, "Internal inconsistency: unknown type field");
			}
		}
	}
	for (i = 0, dp = s_MDPressureAttrDefTable; dp->name != NULL; i++, dp++) {
		if (dp->id == aid) {
			char *pp;
			MDPressureArena *pres = arena->pressure;
			if (pres == NULL)
				arena->pressure = pres = pressure_new();
			pp = (char *)pres + dp->offset;
			switch (dp->type) {
				case 'i':
					*((Int *)pp) = NUM2INT(rb_Integer(val));
					return val;
				case 'f':
					*((Double *)pp) = NUM2DBL(rb_Float(val));
					return val;
				case 'X':
					/*  Isotropic pressure only  */
					val = rb_ary_to_ary(val);
					memset(pres->apply, 0, sizeof(Mat33));
					for (j = 0; j < 3 && j < RARRAY_LEN(val); j++)
						pres->apply[j * 4] = NUM2DBL(rb_Float(RARRAY_PTR(val)[j]));
					return val;
				case 'Y':
					val = rb_ary_to_ary(val);
					for (j = 0; j < 8; j++) {
						if (j < RARRAY_LEN(val))
							pres->cell_flexibility[j] = NUM2DBL(rb_Float(RARRAY_PTR(val)[j]));
						else pres->cell_flexibility[j] = 0.0;
					}
					return val;
				case 'S': case 'I': case 'F':
					rb_raise(rb_eMolbyError, "The attribute is read-only");
				default:
					rb_raise(rb_eMolbyError, "Internal inconsistency: unknown type field");
			}
		}
	}
	rb_raise(rb_eMolbyError, "unknown attribute name (%s)", rb_id2name(aid));
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     arena.(attrname)=
 *
 *  Set the attribute value. The name of the attribute is taken from ruby_frame->last_func.
 */
static VALUE
s_MDArena_SetAttr(VALUE self, VALUE val)
{
	int i;
	struct s_MDArenaAttrDef *dp;
	ID aid = ruby_frame->last_func;
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		if (dp->sid == aid)
			return s_MDArena_Set(self, *(dp->symref), val);
	}
	for (i = 0, dp = s_MDPressureAttrDefTable; dp->name != NULL; i++, dp++) {
		if (dp->sid == aid)
			return s_MDArena_Set(self, *(dp->symref), val);
	}
	rb_raise(rb_eMolbyError, "unknown attribute name (%s)", rb_id2name(aid));
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     arena.to_hash
 *
 *  Returns a (frozen) hash that contains the current value for all attribute keys.
 */
static VALUE
s_MDArena_ToHash(VALUE self)
{
	int i;
	VALUE hash;
	struct s_MDArenaAttrDef *dp;
	hash = rb_hash_new();
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		VALUE attr = ID2SYM(dp->id);
		rb_hash_aset(hash, attr, s_MDArena_Get(self, attr));
	}
	rb_obj_freeze(hash);
	return hash;
}

/*
 *  call-seq:
 *     arena.print_surface_area
 *
 *  Print the surface area information to standard output. (for debug)
 */
static VALUE
s_MDArena_PrintSurfaceArea(VALUE self)
{
	MDArena *arena;
	Int i, natoms;
	VALUE retval, outval;
	Data_Get_Struct(self, MDArena, arena);
	if (arena->sp_arena == NULL)
		rb_raise(rb_eMolbyError, "surface potential is not available");
	natoms = arena->mol->natoms;
	retval = rb_str_new2("Atom     area    energy         forcex1000\n");
	for (i = 0; i < natoms; i++) {
		char buf[256];
		Vector f = arena->forces[kSurfaceIndex * natoms + i];
		Double area = arena->sp_arena->atom_area[i];
		Double energy = area * arena->sp_arena->atom_pot[i];
		VecScaleSelf(f, INTERNAL2KCAL * 1000);
		energy *= INTERNAL2KCAL;
		snprintf(buf, sizeof buf, "%5d %11.5f %11.8f %11.5f %11.5f %11.5f\n",
			   i+1, area, energy, f.x, f.y, f.z);
		rb_str_cat(retval, buf, strlen(buf));
	}
	outval = rb_gv_get("$stdout");
	rb_funcall(outval, rb_intern("write"), 1, retval);
	return self;
}

void
Init_MolbyMDTypes(void)
{
	int i;
	struct s_MDArenaAttrDef *dp;
	char name[40];

	/*  class MDArena  */
	rb_cMDArena = rb_define_class("MDArena", rb_cObject);
/*	rb_define_alloc_func(rb_cMDArena, s_MDArena_Alloc);
    rb_define_private_method(rb_cMDArena, "initialize", s_MDArena_Initialize, 1); */
	rb_define_method(rb_cMDArena, "run", s_MDArena_Run, 1);
    rb_define_method(rb_cMDArena, "minimize", s_MDArena_Minimize, 1);
    rb_define_method(rb_cMDArena, "prepare", s_MDArena_Prepare, -1);
    rb_define_method(rb_cMDArena, "energies", s_MDArena_Energies, 0);
	rb_define_method(rb_cMDArena, "[]", s_MDArena_Get, 1);
	rb_define_method(rb_cMDArena, "[]=", s_MDArena_Set, 2);
	rb_define_method(rb_cMDArena, "to_hash", s_MDArena_ToHash, 0);
	rb_define_method(rb_cMDArena, "print_surface_area", s_MDArena_PrintSurfaceArea, 0);

	/*  All setter and getter are handled with the same C function (attribute name is taken
	    from ruby_frame)  */
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		rb_define_method(rb_cMDArena, dp->name, s_MDArena_GetAttr, 0);
		strncpy(name, dp->name, 38);
		strcat(name, "=");
		rb_define_method(rb_cMDArena, name, s_MDArena_SetAttr, 1);
		dp->id = rb_intern(dp->name);
		dp->sid = rb_intern(name);
		*(dp->symref) = ID2SYM(dp->id);
	}
	for (i = 0, dp = s_MDPressureAttrDefTable; dp->name != NULL; i++, dp++) {
		rb_define_method(rb_cMDArena, dp->name, s_MDArena_GetAttr, 0);
		strncpy(name, dp->name, 38);
		strcat(name, "=");
		rb_define_method(rb_cMDArena, name, s_MDArena_SetAttr, 1);
		dp->id = rb_intern(dp->name);
		dp->sid = rb_intern(name);
		*(dp->symref) = ID2SYM(dp->id);
	}
}
