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
	Data_Get_Struct(self, MDArena, arena);
	nsteps = NUM2INT(rb_Integer(arg));
	if (arena->is_running)
		rb_raise(rb_eMolbyError, "the simulation is already running. You cannot do simulation recursively.");
	if (nsteps < 0)
		rb_raise(rb_eMolbyError, "the number of steps should be non-negative integer");
	
	if (arena->is_initialized < 2 || arena->xmol->needsMDRebuild) {
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
	retval = md_main(arena, minimize);

	if (retval == 0 || arena->step > start_step) {
		int i, natoms = arena->xmol->natoms;
		Atom *ap;
		Vector *vp = (Vector *)malloc(sizeof(Vector) * (natoms > 4 ? natoms : 4));
		IntGroup *ig = IntGroupNewWithPoints(0, natoms, -1);
		int copy_cell = 0;
		if (arena->pressure != NULL && arena->pressure->disabled == 0)
			copy_cell = 1;
#if MINIMIZE_CELL
		if (arena->minimize_cell && minimize)
			copy_cell = 1;
#endif
		if (arena->step > start_step) {
			/*  Copy coordinates and velocities from mol to xmol */
			for (i = 0, ap = arena->mol->atoms; i < natoms; i++, ap = ATOM_NEXT(ap))
				vp[i] = ap->r;
			MolActionCreateAndPerform(arena->xmol, gMolActionSetAtomPositions, ig, natoms, vp);
			for (i = 0, ap = arena->mol->atoms; i < natoms; i++, ap = ATOM_NEXT(ap))
				vp[i] = ap->v;
			MolActionCreateAndPerform(arena->xmol, gMolActionSetAtomVelocities, ig, natoms, vp);
			if (copy_cell && arena->mol->cell != NULL) {
				vp[0] = arena->mol->cell->axes[0];
				vp[1] = arena->mol->cell->axes[1];
				vp[2] = arena->mol->cell->axes[2];
				vp[3] = arena->mol->cell->origin;
				MolActionCreateAndPerform(arena->xmol, gMolActionSetBox, vp, vp + 1, vp + 2, vp + 3, -1, 0);
			}
		}
		/*  Copy forces (this is valid even for "zero-step" run)  */
		for (i = 0, ap = arena->mol->atoms; i < natoms; i++, ap = ATOM_NEXT(ap))
			vp[i] = ap->f;
		MolActionCreateAndPerform(arena->xmol, gMolActionSetAtomForces, ig, natoms, vp);
		free(vp);
		IntGroupRelease(ig);
	}

	if (retval != 0)
		rb_raise(rb_eMolbyError, "MD Error: %s", arena->errmsg);

	if (minimize && arena->minimize_complete && rb_block_given_p())
		rb_yield(self);

	return self;
}	

/*
 *  call-seq:
 *     run(n)       -> self
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
 *     minimize(n)       -> self
 *     minimize(n) { ... } -> self
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
 *     prepare(check_only = false)         -> self or nil
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
	Int check_only = 0, status;
	char *msg;
	VALUE fval;
	Data_Get_Struct(self, MDArena, arena);
	rb_scan_args(argc, argv, "01", &fval);
	if (RTEST(fval))
		check_only = 1;
	status = MoleculePrepareMDArena(arena->xmol, check_only, &msg);
	if (status < 0) {
		/*  Exception object is created first to have a chance to do free(msg)  */
		VALUE exval = rb_exc_new2(rb_eMolbyError, msg);
		free(msg);
		rb_exc_raise(exval);
	} else if (status > 0) {
		free(msg);
		return Qnil;
	} else return self;
	
#if 0
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
		mol->needsMDRebuild = 0;  /*  We know the "modified" parameters are consistent with the MDArena  */
	}

	if (missing)
		return Qnil;
	else return self;
#endif
}

/*
 *  call-seq:
 *     energies -> [total, bond, angle, dihedral, improper, vdw, electrostatic, auxiliary, surface, kinetic, net [, ewald]
 *
 *  Get the current energies.
 */
static VALUE
s_MDArena_Energies(VALUE self)
{
	MDArena *arena;
	VALUE val;
	int i;
	static Int s_indices[] = {kBondIndex, kAngleIndex, kDihedralIndex, kImproperIndex, kVDWIndex, kElectrostaticIndex, kAuxiliaryIndex, kSurfaceIndex, kKineticIndex };
	Data_Get_Struct(self, MDArena, arena);
	val = rb_ary_new();
	rb_ary_push(val, rb_float_new(arena->total_energy * INTERNAL2KCAL));
	for (i = 0; i < sizeof(s_indices) / sizeof(s_indices[0]); i++)
		rb_ary_push(val, rb_float_new(arena->energies[s_indices[i]] * INTERNAL2KCAL));
	rb_ary_push(val, rb_float_new((arena->energies[kKineticIndex] + arena->total_energy) * INTERNAL2KCAL));
	if (arena->use_ewald)
		rb_ary_push(val, rb_float_new((arena->energies[kESCorrectionIndex] + arena->energies[kPMEIndex]) * INTERNAL2KCAL));
	return val;
}

static VALUE s_LogFileSym, s_CoordFileSym, s_VelFileSym, s_ForceFileSym, 
s_DebugFileSym, s_DebugOutputLevelSym, s_StepSym, s_CoordOutputFreqSym, 
s_EnergyOutputFreqSym, s_CoordFrameSym, s_TimestepSym, s_CutoffSym, 
s_ElectroCutoffSym, s_PairlistDistanceSym, s_SwitchDistanceSym, s_TemperatureSym, s_TransientTempSym, 
s_AverageTempSym, s_AndersenFreqSym, s_AndersenCouplingSym, s_RandomSeedSym, 
s_DielectricSym, s_GradientConvergenceSym, s_CoordinateConvergenceSym, s_UseXplorShiftSym, 
s_Scale14VdwSym, s_Scale14ElectSym, s_RelocateCenterSym, 
s_SurfaceProbeRadiusSym, s_SurfaceTensionSym, s_SurfacePotentialFreqSym, s_UseGraphiteSym,
s_AlchemicalLambdaSym, s_AlchemicalDeltaLambdaSym, s_AlchemicalEnergySym, s_MinimizeCellSym,
s_UseEwaldSym, s_EwaldBetaSym, s_EwaldGridSym, s_EwaldFreqSym, s_EwaldOrderSym;

struct s_MDArenaAttrDef {
	char *name;
	VALUE *symref;  /*  Address of s_LogFileSym etc. */
	ID id;			/*  Ruby ID of the symbol; will be set within Init_MolbyMDTypes()  */
	ID sid;         /*  Ruby ID of the symbol plus '='; will be set within Init_MolbyMDTypes()  */
	char type;      /*  s: string (const char *), i: Int, f: Double, e: Double in energy dimension (unit conversion is necessary). */
					/*  Uppercase: read-only.  */
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
	{"switch_distance",   &s_SwitchDistanceSym,   0, 0, 'f', offsetof(MDArena, switch_distance)},
	{"temperature",       &s_TemperatureSym,      0, 0, 'f', offsetof(MDArena, temperature)},
	{"transient_temperature", &s_TransientTempSym, 0, 0, 'F', offsetof(MDArena, transient_temperature)},
	{"average_temperature", &s_AverageTempSym,    0, 0, 'F', offsetof(MDArena, average_temperature)},
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
	{"alchemical_lambda", &s_AlchemicalLambdaSym, 0, 0, 'f', offsetof(MDArena, alchem_lambda)},
	{"alchemical_delta_lambda", &s_AlchemicalDeltaLambdaSym, 0, 0, 'f', offsetof(MDArena, alchem_dlambda)},
	{"alchemical_energy", &s_AlchemicalEnergySym, 0, 0, 'E', offsetof(MDArena, alchem_energy)},
	{"minimize_cell",     &s_MinimizeCellSym,     0, 0, 'b', offsetof(MDArena, minimize_cell)},
	{"use_ewald",         &s_UseEwaldSym,         0, 0, 'i', offsetof(MDArena, use_ewald)},
	{"ewald_beta",        &s_EwaldBetaSym,        0, 0, 'f', offsetof(MDArena, ewald_beta)},
	{"ewald_grid",        &s_EwaldGridSym,        0, 0, 0,   offsetof(MDArena, ewald_grid_x)},
	{"ewald_freq",        &s_EwaldFreqSym,        0, 0, 'i', offsetof(MDArena, ewald_freq)},
	{"ewald_order",       &s_EwaldOrderSym,       0, 0, 'i', offsetof(MDArena, ewald_order)},
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
 *     self[attr]
 *
 *  Get the attribute value. The attribute values also can be accessed by self.attribute_name.
 */
static VALUE
s_MDArena_Get(VALUE self, VALUE attr)
{
	MDArena *arena;
	int i;
	struct s_MDArenaAttrDef *dp;
	ID aid = rb_to_id(attr);
	Data_Get_Struct(self, MDArena, arena);
	if (aid == SYM2ID(s_EwaldGridSym)) {
		/*  Array of three grid values  */
		return rb_ary_new3(3, INT2NUM(arena->ewald_grid_x), INT2NUM(arena->ewald_grid_y), INT2NUM(arena->ewald_grid_z));
	}
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
				case 'b':
				case 'B':
					return INT2NUM((Int)(*((Byte *)p)));
				case 'i':
				case 'I':
					return INT2NUM(*((Int *)p));
				case 'f':
				case 'F':
					return rb_float_new(*((Double *)p));
				case 'e':
				case 'E':
					return rb_float_new(*((Double *)p) * INTERNAL2KCAL);
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
				case 'b':
				case 'B':
					return INT2NUM((Int)(*((Byte *)pp)));
				case 'f':
				case 'F':
					return rb_float_new(*((Double *)pp));
				case 'e':
				case 'E':
					return rb_float_new(*((Double *)pp) * INTERNAL2KCAL);
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

static VALUE
s_MDArena_GetAttr(VALUE self)
{
	return s_MDArena_Get(self, ID2SYM(ruby_frame->last_func));
}

/*
 *  call-seq:
 *     self[attr] = value
 *
 *  Set the attribute value. The attributes can also be modified by self.attribute_name = value.
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
	if (aid == SYM2ID(s_EwaldGridSym)) {
		if (rb_obj_is_kind_of(val, rb_cNumeric)) {
			i = NUM2INT(rb_Integer(val));
			if (i <= 0)
				rb_raise(rb_eMolbyError, "The ewald grid must be positive integer");
			arena->ewald_grid_x = arena->ewald_grid_y = arena->ewald_grid_z = i;
		} else {
			int ival[3];
			val = rb_ary_to_ary(val);
			j = RARRAY_LEN(val);
			if (j < 3)
				rb_raise(rb_eMolbyError, "The ewald grid must be an integer or an array of three integers");
			for (i = 0; i < 3; i++) {
				ival[i] = NUM2INT(rb_Integer(RARRAY_PTR(val)[i]));
				if (ival[i] <= 0)
					rb_raise(rb_eMolbyError, "The ewald grid must be positive integer");
			}
			arena->ewald_grid_x = ival[0];
			arena->ewald_grid_y = ival[1];
			arena->ewald_grid_z = ival[2];
		}
		return val;
	}
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
				case 'b':
					*((Byte *)p) = NUM2INT(rb_Integer(val));
					return val;
				case 'f':
					*((Double *)p) = NUM2DBL(rb_Float(val));
					return val;
				case 'e':
					*((Double *)p) = NUM2DBL(rb_Float(val) * KCAL2INTERNAL);
					return val;
				case 'S': case 'I': case 'B': case 'F': case 'E':
					rb_raise(rb_eMolbyError, "The attribute '%s' is read-only", rb_id2name(aid));
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
				case 'b':
					*((Byte *)pp) = NUM2INT(rb_Integer(val));
					return val;
				case 'f':
					*((Double *)pp) = NUM2DBL(rb_Float(val));
					return val;
				case 'e':
					*((Double *)pp) = NUM2DBL(rb_Float(val) * KCAL2INTERNAL);
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
				case 'S': case 'I': case 'B': case 'F': case 'E':
					rb_raise(rb_eMolbyError, "The attribute '%s' is read-only", rb_id2name(aid));
				default:
					rb_raise(rb_eMolbyError, "Internal inconsistency: unknown type field");
			}
		}
	}
	rb_raise(rb_eMolbyError, "unknown attribute name (%s)", rb_id2name(aid));
	return Qnil;  /*  Not reached  */
}

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
 *     to_hash -> Hash
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
 *     set_alchemical_perturbation(group1, group2) -> [group1, group2]
 *
 *  Set vanishing and appearing atom groups for alchemical perturbation.
 */
static VALUE
s_MDArena_SetAlchemicalPerturbation(VALUE self, VALUE gval1, VALUE gval2)
{
	IntGroup *ig1, *ig2;
	MDArena *arena;
	char *flags;
	int i, n;
	Data_Get_Struct(self, MDArena, arena);
	if (gval1 == Qnil && gval2 == Qnil) {
		md_set_alchemical_flags(arena, 0, NULL);
		return Qnil;
	}
	if (arena->mol == NULL)
		rb_raise(rb_eMolbyError, "Molecule is not set");
	n = arena->xmol->natoms;
	flags = (char *)calloc(1, n);
	ig1 = (gval1 == Qnil ? NULL : IntGroupFromValue(gval1));
	ig2 = (gval2 == Qnil ? NULL : IntGroupFromValue(gval2));
	for (i = 0; i < n; i++) {
		if (ig1 != NULL && IntGroupLookupPoint(ig1, i) >= 0)
			flags[i] = 1;
		if (ig2 != NULL && IntGroupLookupPoint(ig2, i) >= 0) {
			if (flags[i] == 1)
				rb_raise(rb_eMolbyError, "duplicate atom (%d) in vanishing and appearing groups", i);
			flags[i] = 2;
		}
	}
	if (md_set_alchemical_flags(arena, n, flags) != 0)
		rb_raise(rb_eMolbyError, "cannot set alchemical flags");
	free(flags);
	if (ig1 != NULL) {
		gval1 = ValueFromIntGroup(ig1);
		IntGroupRelease(ig1);
	}
	if (ig2 != NULL) {
		gval2 = ValueFromIntGroup(ig2);
		IntGroupRelease(ig2);
	}
	return rb_ary_new3(2, gval1, gval2);
}

/*
 *  call-seq:
 *     get_alchemical_perturbation -> [group1, group2] or nil
 *
 *  If alchemical perturbation is enabled, get the current vanishing and appearing atom groups.
 *  Otherwise, return nil.
 */
static VALUE
s_MDArena_GetAlchemicalPerturbation(VALUE self)
{
	IntGroup *ig1, *ig2;
	VALUE gval1, gval2;
	MDArena *arena;
	int i;
	Data_Get_Struct(self, MDArena, arena);
	if (arena->nalchem_flags == 0)
		return Qnil;
	if (arena->mol == NULL)
		rb_raise(rb_eMolbyError, "Molecule is not set");	
	ig1 = IntGroupNew();
	ig2 = IntGroupNew();
	for (i = 0; i < arena->nalchem_flags; i++) {
		if (arena->alchem_flags[i] == 1)
			IntGroupAdd(ig1, i, 1);
		else if (arena->alchem_flags[i] == 2)
			IntGroupAdd(ig2, i, 1);
	}
	gval1 = ValueFromIntGroup(ig1);
	gval2 = ValueFromIntGroup(ig2);
	IntGroupRelease(ig1);
	IntGroupRelease(ig2);
	return rb_ary_new3(2, gval1, gval2);	
}

/*
 *  call-seq:
 *    set_external_forces(ary) -> self
 *
 *  Set external forces. Ary should be an array of objects that can be converted to Vector3D.
 */
static VALUE
s_MDArena_SetExternalForces(VALUE self, VALUE aval)
{
	Vector *vp;
	MDArena *arena;
	int i, n;
	Data_Get_Struct(self, MDArena, arena);
	if (arena->mol == NULL)
		rb_raise(rb_eMolbyError, "Molecule is not set");
	if (aval == Qnil) {
		md_set_external_forces(arena, 0, NULL);
		return self;
	}
	aval = rb_ary_to_ary(aval);
	n = RARRAY_LEN(aval);
	if (n == 0) {
		md_set_external_forces(arena, 0, NULL);
		return self;
	}
	vp = (Vector *)calloc(sizeof(Vector), n);
	for (i = 0; i < n; i++) {
		VectorFromValue(RARRAY_PTR(aval)[i], vp + i);
		VecScaleSelf(vp[i], KCAL2INTERNAL);
	}
	md_set_external_forces(arena, n, vp);
	free(vp);
	return self;
}

/*
 *  call-seq:
 *    get_external_force(index) -> Vector3D or nil
 *
 *  Get the current external force for the atom. If the external force is not set, nil is returned.
 */
static VALUE
s_MDArena_GetExternalForce(VALUE self, VALUE ival)
{
	int i;
	VALUE vval;
	Vector v;
	MDArena *arena;
	Data_Get_Struct(self, MDArena, arena);
	if (arena->mol == NULL)
		rb_raise(rb_eMolbyError, "Molecule is not set");
	i = NUM2INT(rb_Integer(ival));
	if (i < 0 || i >= arena->nexforces)
		return Qnil;
	v = arena->exforces[i];
	VecScaleSelf(v, INTERNAL2KCAL);
	vval = ValueFromVector(&v);
	return vval;
}

/*
 *  call-seq:
 *     init_velocities([temperature]) -> self
 *
 *  Give random (Boltzmann-weighted) velocities to all atoms. If temperature is given,
 *  it is set before giving velocities.
 */
static VALUE
s_MDArena_InitVelocities(int argc, VALUE *argv, VALUE self)
{
	MDArena *arena;
	VALUE tval;
	Data_Get_Struct(self, MDArena, arena);
	rb_scan_args(argc, argv, "01", &tval);
	if (tval != Qnil)
		s_MDArena_Set(self, s_TemperatureSym, tval);
	md_init_velocities(arena);
	return self;
}

/*
 *  call-seq:
 *     scale_velocities([temperature]) -> self
 *
 *  Scale atom velocities to match the current temperature. If temperature is given,
 *  it is set before giving velocities.
 */
static VALUE
s_MDArena_ScaleVelocities(int argc, VALUE *argv, VALUE self)
{
	MDArena *arena;
	VALUE tval;
	Data_Get_Struct(self, MDArena, arena);
	rb_scan_args(argc, argv, "01", &tval);
	if (tval != Qnil)
		s_MDArena_Set(self, s_TemperatureSym, tval);
	md_scale_velocities(arena);
	return self;
}

/*
 *  call-seq:
 *     keys -> Array
 *
 *  Returns an array of valid attributes.
 */
static VALUE
s_MDArena_Keys(VALUE self)
{
	int i;
	VALUE ary;
	struct s_MDArenaAttrDef *dp;
	ary = rb_ary_new();
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		VALUE attr = ID2SYM(dp->id);
		rb_ary_push(ary, attr);
	}
	return ary;
}

/*
 *  call-seq:
 *     print_surface_area
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

/*
 *  call-seq:
 *     bond_par(idx) -> ParameterRef
 *
 *  Returns a parameter that is used for the idx-th bond.
 */
static VALUE
s_MDArena_BondPar(VALUE self, VALUE val)
{
	MDArena *arena;
	Int i, j;
	Data_Get_Struct(self, MDArena, arena);
	i = NUM2INT(rb_Integer(val));
	if (arena->xmol->needsMDRebuild || arena->par == NULL)
		md_prepare(arena, 1);
	if (i < 0 || i >= arena->xmol->nbonds) {
		rb_raise(rb_eMolbyError, "bond index (%d) out of range (0...%d)", i, arena->xmol->nbonds);
	}
	j = arena->bond_par_i[i];
	if (j < 0 || j >= arena->par->nbondPars)
		return Qnil;  /*  No parameter assigned  */
	return ValueFromMoleculeWithParameterTypeAndIndex(arena->xmol, kBondParType, -j - 1);
}

/*
 *  call-seq:
 *     angle_par(idx) -> ParameterRef
 *
 *  Returns a parameter that is used for the idx-th angle.
 */
static VALUE
s_MDArena_AnglePar(VALUE self, VALUE val)
{
	MDArena *arena;
	Int i, j;
	Data_Get_Struct(self, MDArena, arena);
	i = NUM2INT(rb_Integer(val));
	if (arena->xmol->needsMDRebuild || arena->par == NULL)
		md_prepare(arena, 1);
	if (i < 0 || i >= arena->xmol->nangles) {
		rb_raise(rb_eMolbyError, "angle index (%d) out of range (0...%d)", i, arena->xmol->nangles);
	}
	j = arena->angle_par_i[i];
	if (j < 0 || j >= arena->par->nanglePars)
		return Qnil;  /*  No parameter assigned  */
	return ValueFromMoleculeWithParameterTypeAndIndex(arena->xmol, kAngleParType, -j - 1);
}

/*
 *  call-seq:
 *     dihedral_par(idx) -> ParameterRef
 *
 *  Returns a parameter that is used for the idx-th dihedral.
 */
static VALUE
s_MDArena_DihedralPar(VALUE self, VALUE val)
{
	MDArena *arena;
	Int i, j;
	Data_Get_Struct(self, MDArena, arena);
	i = NUM2INT(rb_Integer(val));
	if (arena->xmol->needsMDRebuild || arena->par == NULL)
		md_prepare(arena, 1);
	if (i < 0 || i >= arena->xmol->ndihedrals) {
		rb_raise(rb_eMolbyError, "dihedral index (%d) out of range (0...%d)", i, arena->xmol->ndihedrals);
	}
	j = arena->dihedral_par_i[i];
	if (j < 0 || j >= arena->par->ndihedralPars)
		return Qnil;  /*  No parameter assigned  */
	return ValueFromMoleculeWithParameterTypeAndIndex(arena->xmol, kDihedralParType, -j - 1);
}

/*
 *  call-seq:
 *     improper_par(idx) -> ParameterRef
 *
 *  Returns a parameter that is used for the idx-th improper.
 */
static VALUE
s_MDArena_ImproperPar(VALUE self, VALUE val)
{
	MDArena *arena;
	Int i, j;
	Data_Get_Struct(self, MDArena, arena);
	i = NUM2INT(rb_Integer(val));
	if (arena->xmol->needsMDRebuild || arena->par == NULL)
		md_prepare(arena, 1);
	if (i < 0 || i >= arena->xmol->nimpropers) {
		rb_raise(rb_eMolbyError, "improper index (%d) out of range (0...%d)", i, arena->xmol->nimpropers);
	}
	j = arena->improper_par_i[i];
	if (j < 0 || j >= arena->par->nimproperPars)
		return Qnil;  /*  No parameter assigned  */
	return ValueFromMoleculeWithParameterTypeAndIndex(arena->xmol, kImproperParType, -j - 1);
}

/*
 *  call-seq:
 *     vdw_par(idx) -> ParameterRef
 *
 *  Returns a vdw parameter that is used for the idx-th atom.
 */
static VALUE
s_MDArena_VdwPar(VALUE self, VALUE val)
{
	MDArena *arena;
	Int i, j;
	Data_Get_Struct(self, MDArena, arena);
	i = NUM2INT(rb_Integer(val));
	if (arena->xmol->needsMDRebuild || arena->par == NULL)
		md_prepare(arena, 1);
	if (i < 0 || i >= arena->xmol->natoms) {
		rb_raise(rb_eMolbyError, "atom index (%d) out of range (0...%d)", i, arena->xmol->natoms);
	}
	j = arena->vdw_par_i[i];
	if (j < 0 || j >= arena->par->nvdwPars)
		return Qnil;  /*  No parameter assigned  */
	return ValueFromMoleculeWithParameterTypeAndIndex(arena->xmol, kVdwParType, -j - 1);
}

static VALUE
s_MDArena_testPME(int argc, VALUE *argv, VALUE self)
{
	extern pme_test(MDArena *);
	MDArena *arena;
	if (rb_obj_is_kind_of(self, rb_cMDArena)) {
		Data_Get_Struct(self, MDArena, arena);
		if (argc == 0 || RTEST(argv[0])) {
			arena->use_ewald = 1;
			arena->xmol->needsMDRebuild = 1;
			if (argc > 0 && rb_obj_is_kind_of(argv[0], rb_cNumeric)) {
				Int n = NUM2INT(rb_Integer(argv[0]));
				if (n > 1) {
					arena->ewald_grid_x = n;
					arena->ewald_grid_y = n;
					arena->ewald_grid_z = n;
				} else {
					arena->use_ewald = 2;  /*  Direct Ewald  */
				}
			}
		} else {
			arena->use_ewald = 0;
			arena->xmol->needsMDRebuild = 1;
		}
		s_MDArena_Prepare(0, NULL, self);
	} else {
		rb_raise(rb_eMolbyError, "class method MDArena.test_pme is not supported now");
	}
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
	rb_define_method(rb_cMDArena, "keys", s_MDArena_Keys, 0);
	rb_define_method(rb_cMDArena, "set_alchemical_perturbation", s_MDArena_SetAlchemicalPerturbation, 2);
	rb_define_method(rb_cMDArena, "get_alchemical_perturbation", s_MDArena_GetAlchemicalPerturbation, 0);
	rb_define_method(rb_cMDArena, "set_external_forces", s_MDArena_SetExternalForces, 1);
	rb_define_method(rb_cMDArena, "get_external_force", s_MDArena_GetExternalForce, 1);
	rb_define_method(rb_cMDArena, "init_velocities", s_MDArena_InitVelocities, -1);
	rb_define_method(rb_cMDArena, "scale_velocities", s_MDArena_ScaleVelocities, -1);
	rb_define_method(rb_cMDArena, "print_surface_area", s_MDArena_PrintSurfaceArea, 0);

	rb_define_method(rb_cMDArena, "bond_par", s_MDArena_BondPar, 1);
	rb_define_method(rb_cMDArena, "angle_par", s_MDArena_AnglePar, 1);
	rb_define_method(rb_cMDArena, "dihedral_par", s_MDArena_DihedralPar, 1);
	rb_define_method(rb_cMDArena, "improper_par", s_MDArena_ImproperPar, 1);
	rb_define_method(rb_cMDArena, "vdw_par", s_MDArena_VdwPar, 1);

	rb_define_method(rb_cMDArena, "pme_test", s_MDArena_testPME, -1);
//	rb_define_singleton_method(rb_cMDArena, "pme_test", s_MDArena_testPME, -1);

	/*  All setter and getter are handled with the same C function (attribute name is taken
	    from ruby_frame)  */
	for (i = 0, dp = s_MDArenaAttrDefTable; dp->name != NULL; i++, dp++) {
		rb_define_method(rb_cMDArena, dp->name, s_MDArena_GetAttr, 0);
		strncpy(name, dp->name, 38);
		name[38] = 0;
		strcat(name, "=");
		rb_define_method(rb_cMDArena, name, s_MDArena_SetAttr, 1);
		dp->id = rb_intern(dp->name);
		dp->sid = rb_intern(name);
		*(dp->symref) = ID2SYM(dp->id);
	}
	for (i = 0, dp = s_MDPressureAttrDefTable; dp->name != NULL; i++, dp++) {
		rb_define_method(rb_cMDArena, dp->name, s_MDArena_GetAttr, 0);
		strncpy(name, dp->name, 38);
		name[38] = 0;
		strcat(name, "=");
		rb_define_method(rb_cMDArena, name, s_MDArena_SetAttr, 1);
		dp->id = rb_intern(dp->name);
		dp->sid = rb_intern(name);
		*(dp->symref) = ID2SYM(dp->id);
	}
}
