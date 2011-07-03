/*
 *  MDCore.h
 *
 *  Created by Toshi Nagata on 2005/06/06.
 *  Copyright 2005 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __MDCORE_H__
#define __MDCORE_H__

#include "../Molecule.h"
#include "../Parameter.h"

#include <stdarg.h>
#include <time.h>
#include <setjmp.h>

#ifdef __cplusplus
extern "C" {
#endif
	
#define UNITCHARGE 1.602e-19
/*  Vacuum permittability in internal units  */
#define EPSILON0 0.572798
/* = 8.854e-12*1e-10/(6.023e23*1e-3*KJ2INTERNAL)/(UNITCHARGE*UNITCHARGE)  */
/* = 8.854e-12 (C^2 kg^-1 m^-3 s^2) * (UNITCHARGE C/e)^(-2) * (1e10 ang/m)^(-3) * (6.02e26 am/kg)^(-1) * (1e15 fs/s)^2  */

#define COULOMBIC (1/(4*PI*EPSILON0))

/*  1 bar = 1e5 Pa = 1e2 kJ m^(-3) = (1e2*6.023e23*KJ2INTERNAL)internal_eu * (1e10 ang)^(-3) = 1e-32*6.023e23 internal_eu ang^(-3)  */
#define BAR2INTERNALPRESSURE (6.02e-9)

#define ERROR_out_of_memory "Low memory"

/*  Index enumeration for storing partial energy/force  */
enum {
	kBondIndex = 0,
	kAngleIndex,
	kDihedralIndex,
	kImproperIndex,
	kVDWIndex,
	kElectrostaticIndex,
	kAuxiliaryIndex,
	kSurfaceIndex,
	kKineticIndex,
	kEndIndex,
	kSlowIndex = kSurfaceIndex
};

/*  Type of custom bond parameter  */
enum {
	kTCLProcType = 255,
	kMorseType = 0
};

typedef struct MDCustomPar {
	Byte   type;  /* 0:bond, 1:angle, 2:dihedral, 3:improper, 4:vdw */
	Int    n1, n2, n3, n4;  /*  Atom indices  */
	Int    index; /* index of bond_pars, angle_pars, etc. */
} MDCustomPar;

/*  Verlet list item (for non-bonding force)  */
typedef struct MDVerlet {
	Int    n1, n2;
	Symop symop;         /*  The symmetry operation for the 2nd atom.
	 Only the translational components are used
	 in the force calculations, because the non-
	 translational components should be handled
	 in amend_by_symmetry() operation.  */
	Int    mult;                /*  How many equivalent interactions will appear */
	Byte   vdw_type;            /*  0: vdw, 1: scaled (1-4) vdw  */
	unsigned int index;         /*  The index to arena->vdw_cache  */
	Double  vdw_cutoff;          /*  Specific vdw cutoff  */
	Double  length;
	/*	Double  vcut;   */           /*  Value at the specific cutoff distance  */
} MDVerlet;

/*  MDExclusion list  */
/*  exlist[exinfo[i].index0 .. exinfo[i].index1-1] : special exclusion  */
/*  exlist[exinfo[i].index1 .. exinfo[i].index2-1] : 1-2 (directly bonded)  */
/*  exlist[exinfo[i].index2 .. exinfo[i].index3-1] : 1-3  */
/*  exlist[exinfo[i].index3 .. exinfo[i+1].index0-1] : 1-4  */
/*  exlist[] is a (long) array of Int; exinfo[] is MDExclusion[natoms+1]; +1 needed because exinfo[natoms].index0 must be accessible  */
/*  1-4 comes last because it sometimes needs special treatments  */
typedef struct MDExclusion {
	Int index0, index1, index2, index3;
} MDExclusion;

/*  van der Waals parameter cache  */
typedef struct MDVdwCache {
	VdwPairPar par;
	Double vcut, vcut14;  /*  The values at the cutoff distance  */
} MDVdwCache;

/*  Custom bond parameter  */
typedef struct MDCustomBondPar {
	Byte type;
	union {
		void *proc;  /*  TCL callback proc (stored as a TCL object)  */
		/*  The TCL callback proc takes one parameter (interatomic distance) and returns a list of two floating numbers (the energy and force). Force = -d(energy)/dr; note the minus sign. */
		struct { /*  The Morse potential  */
			/*  V(r) = D(1-s)^2, F(r) = -2aDs(1-s); s = exp(-a(r-r0)) */
			Double D, a, r0;
		} morse;
	} u;
} MDCustomBondPar;

/*  Snapshots  */
typedef struct MDSnapshot {
	Int    step;
	Int    natoms;
	Vector rvf[3];    /*  Variable length array; r[0], v[0], f[0], r[1], ... */
} MDSnapshot;

/*  Group flags  */
typedef struct MDGroupFlags {
	unsigned char *flags;
	int natoms;
} MDGroupFlags;
#define group_flag_mask "\001\002\004\010\020\040\100\200"
#define get_group_flag(gf, n) (n < (gf)->natoms ? (gf)->flags[(n)/8] & group_flag_mask[(n)%8] : 0)
#define set_group_flag(gf, n, bool) (n < (gf)->natoms ? (bool ? ((gf)->flags[(n)/8] |= group_flag_mask[(n)%8]) : ((gf)->flags[(n)/8] &= ~(group_flag_mask[(n)%8]))) : 0)

/*  Ring buffer for sending coordinates to another thread  */
typedef struct MDRing {
	Vector *buf;           /*  Vector buffer  */
	Int    size;           /*  Number of vectors per frame; natoms + (use_cell ? 4 : 0)  */
	Int    nframes;        /*  Number of frames in the ring buffer (2000 / natoms)  */
	Int    next;           /*  Next frame index to store data  */
	Int    count;          /*  Number of frames currently in the ring buffer  */
} MDRing;
	
/*  Everything needed for MD  */
typedef struct MDArena {
	Int refCount;
	Molecule *xmol;         /*  Molecule given from the outside world  */
	Molecule *mol;          /*  Private copy of the Molecule for MD calc  */

	const char *coord_input_name;
	FILE *coord_input;
	Byte  coord_input_type;
	Int   coord_input_lineno;
	Int   coord_input_frame;
	void *coord_input_data;  /*  pointer to dcd_header_record etc. */

	/*  Results output streams  */
	const char *log_result_name;
	const char *coord_result_name;
	const char *vel_result_name;
	const char *force_result_name;
	const char *extend_result_name;
	const char *debug_result_name;
	FILE *log_result;
	FILE *coord_result;
	FILE *vel_result;
	FILE *force_result;
	FILE *extend_result;
	FILE *debug_result;
	Int debug_output_level;
	
	/*  Error and interrupt handling  */
	void (*md_panic_func)(struct MDArena *arena, const char *msg);
	int (*md_callback_func)(struct MDArena *arena);
	Int callback_freq;

	char errmsg[256];
	jmp_buf *setjmp_buf;
	
	/*  Time  */
	time_t start_time;
	
	/*  MD parameters  */
	Int step, start_step, end_step;
	Int coord_output_freq;
	Int energy_output_freq;
	Int coord_result_frame;
	Double timestep;
	Double cutoff;
	Double electro_cutoff;
	Double pairlist_distance;
	Double temperature;
	Int rescale_temp_freq;
	Int reinit_temp_freq;
	Int velocities_read;
	Int random_seed;
	Double dielectric;  /*  Bulk dielectric constant of the medium  */
	Double gradient_convergence;    /*  Convergent criterion for energy minimization  */
	Double coordinate_convergence;  /*  Convergent criterion for energy minimization  */
	Int use_xplor_shift; /* Use X-Plor type shift function? (default = 1)  */
	
	Double scale14_vdw;   /* Scaling factor for 1-4 vdw interactions (default = 1)  */
	Double scale14_elect; /* Scaling factor for 1-4 electrostatic interactions (default = 1) */
	
	Int relocate_center; /* Fix center of mass at the original position (default = 1) */
	Int quench_translational_momentum;
	Int quench_angular_momentum;
	Int output_expanded_atoms;  /*  Include symmetry-expanded atoms in output  */
	
	/*  Andersen thermostat  */
	Int andersen_thermo_freq;
	Double andersen_thermo_coupling;
	
	/*  Surface potential  */
	Double probe_radius;  /*  The probe radius for surface area calculation  */
	Double surface_tension; /*  The microscopic surface tension for SASA potential  */
	Int surface_potential_freq;
//	Int sphere_points;   /*  The number of sphere points  */
	
	/*  Attempt to handle _extremely_ distorted structures. */
	Double anbond_thres;  /* The bond with initial length >= (1+anbond_thres)*r0 is an abnormal bond */
	Double anbond_anneal_rate; /* For the abnormal bonds, r0 is set to the initial length, and decreased by this amount every timestep */
	
	/*  Custom bond potential  */
	Int    ncustom_bond_pars;
	MDCustomBondPar *custom_bond_pars;
	Int    *custom_bond_par_i;  /*  Int[nbonds]  */
	/*  0: no custom, >=1: custom_bond_pars[i-1]  */
	
	/*  Custom parameters  */
	Int    ncustom_pars;
	MDCustomPar *custom_pars;
	
	/*  Fix atoms  */
	/*	Int nfix_atoms;
	 fix_atom_record *fix_atoms; */
	
	/*  Switch the bond potential from harmonic to linear at r >= (1 + bond_switch) * r0  */
	/*	Double bond_switch; */
	/*	Double bond_switch_thres; *//* Bond switch is ON only for those bonds whose initial length >= (1 + bond_switch_thres) * r0  */
	
	/*  Artificial, quadratic potential centered at a certain atom - OBSOLETE  */
	/*	Double  centric_potential_force;
	 Int    centric_potential_center;
	 Double  centric_potential_inner_limit;  *//*  potential = 0 for r < inner_limit  */
	/*	Double  centric_potential_outer_limit;  *//*  Switch the potential to linear at r > outer_limit  */
	
	/*  Spherical boundary condition  */
	Double  spherical_bc_force;
	Vector spherical_bc_center;
	Double  spherical_bc_inner_limit; /*  potential = 0 for r < inner_limit  */
	Double  spherical_bc_outer_limit; /*  Switch the potential to linear at r > outer_limit  */
	
	/*  Artificial box potential centered at the origin  */
	/*  The box is [-xsize, -ysize, -zsize]-[xsize, ysize, zsize]  */
	/*  Inside the box potential = 0. Outside the box, a potential k{(x-xsize)^2+(y-ysize)^2+(z-zsize)^2} is applied.  */
	Double  box_potential_xsize, box_potential_ysize, box_potential_zsize;
	Double  box_potential_force;
	
	/*  Graphite potential  */
	Int use_graphite;

	/*  Velocity limit  */
	Double  velocity_limit;  /*  Default = 100  */
	
	/*  TCL interface  */
//	void *tcl_interp;       /*  TCL interpreter; the type is actually (Tcl_Interp *) */
//	void *tcl_retval;       /*  The return value from the callback procedure; set only by md_abort. */
//	const char *callback;   /*  TCL callback procedure  */
	
	/*  Runtime fields  */
	
	/*  Initialize flag  */
	Byte   is_initialized;  /*  0: not initialized, 1: only the static fields (structure-related) are initialized, 2: the runtime fields are initialized (i.e. MD is ready to go)  */
	Byte   is_running;
	
	Byte   request_abort;	/*  If early return is necesarry, assert this flag.  */
	Byte   minimize_complete;  /*  Becomes non-zero if minimization completed.  */

	Int    natoms_uniq;     /*  Number of symmetry-unique atoms  */
	
	MDRing *ring;

	/*  Parameters are copied from mol->par and gBuiltinParameters for each call to md_prepare()  */
	Parameter *par;
	
	/*  Number of missing parameters and suspicious parameters  */
	Int    nmissing;        /*  This must be zero  */
	Int    nsuspicious;     /*  This can be non-zero, but attention should be paid  */
	
	/*  Indices to the parameter records  */
	Int    *bond_par_i;  /*  The index of the bond parameter record: Int[nbonds]  */
	Int    *angle_par_i; /*  The index of the angle parameter record: Int[nangles]  */
	Int    *dihedral_par_i;  /*  The index of the dihedral parameter record: Int[ndihedrals]  */
	Int    *improper_par_i;  /*  The index of the improper parameter record: Int[nimpropers]  */
	
	/*  The van der Waals parameters are precalculated for each pair of atom types  */
	/*  For each pair, the A/B parameters are calculated from the values for the component types  */
	/*  If the pair parameters are specifically given those are used  */
	/*  The precalculated parameters for the atom pair (i,j) are given as 
	    vdw_cache[vdw_par_i[i] * par->nvdwPars + vdw_par_i[j]]  */
	Int    *vdw_par_i;
	MDVdwCache *vdw_cache;
	
	Double  *anbond_r0;   /*  The r0 parameters for abnormal bonds  */
	
	/*  Energies and partial forces  */
	Double  *energies;    /*  bond, angle, dihedral, improper, */
	/*  vdw, electrostatic, surface, auxiliary  */
	Double   total_energy;
	Vector *forces;      /*  forces[idx*natoms+i], idx=bond, angle, etc. i=atom index */
	
	/*  Handling transient temperature  */
	Int    degree_of_freedom;  /*  3*(natoms - (number of fixed atoms))  */
	Double  sum_temperature;    /*  For rescaling temperature  */
	Int    nsum_temperature;   /*  ditto  */
	Double  transient_temperature;
	Double  average_temperature;

	/*  Temporary storage for pair interaction calculations  */
	/*	char   *group_flags; */
//	void *group_flags_1, *group_flags_2;  /*  hold the TCL objects  */
	Double   pair_energies[2];  /* [0]: vdw, [1]: electrostatic */
	Vector *pair_forces;
	
	/*  Symmetry operation  */
//	Int    nsyms;
//	Transform  *syms;
	Transform  *cellsyms;  /*  celltr * syms[i] * rcelltr  */
	
	/*  Unit cell vectors  */
//	Vector cella, cellb, cellc, cello;
//	Transform celltr;  /*  fractional coordinate -> cartesian  */
//	Transform rcelltr; /*  cartesian -> fractional coordinate  */
	
	/*  Use periodic boundary conditions along a/b/c axes?  */
	Int    periodic_a, periodic_b, periodic_c;
	
	/*  Wrap the coordinates to the unit cell. This affects only the output (not internal coordinates)  */
	Int    wrap_coordinates;
	
	/*  Pressure control  */
	Int    pressure_freq;
	struct MDPressureArena *pressure;
	
	/*  Fragment informations  */
	/*  Fragment is a cluster of symmetry-unique atoms that are connected by bonds */
	Int    nfragments;        /*  Number of fragments  */
	Int    *fragment_indices; /*  Int[natoms_uniq]; the fragment index for each atom */
	struct MDFragmentInfo {
		Vector pos;
		Double  mass;
	} *fragment_info;     /*  array[nfragments]; internally used  */
	
	/*  "Image" atoms for non-bonding calculations under periodic boundary conditions */
	Int    nexatoms;          /*  May fluctuate during simulation  */
	Atom  *exatoms;
	Int    nexatoms_current;  /*  Current number of "image" atoms  */
	
	/*  sym_relate[i*nsyms+j] is the index of atom that coincide with atom i transformed by inverse of syms[j] */
	/*	Int   *sym_relate; */
	
	/*  The tolerance (in angstrom) to find symmetry-equivalent atoms. Default = 5e-3  */
	Double sym_tolerance;
	
	/*  Multiplicity of each atom implied by symmetry operation  */
	/*	Int   *sym_mult; */
	
	/*  The following 4 arrays are allocated as a contiguous chunk. Call free() for the first one only. */
	/*  Meaning of the value (v) in sym_****_uniq (**** = bond, angle, dihedral, improper):
	 v < 0: the **** is unique and its multiplicity is -v.
	 v >= 0: the **** is not unique and equivalent to v-th ****.  */
	/*	Int   *sym_bond_uniq;
	 Int   *sym_angle_uniq;
	 Int   *sym_dihedral_uniq;
	 Int   *sym_improper_uniq; */
	
	/*  Information for special positions  */
	Int    nspecial_positions;
	Int    (*special_positions)[ATOMS_MAX_SYMMETRY];
	
	/*  MDExclusion list  */
	Int    nexlist;
	Int    *exlist;
	MDExclusion *exinfo;
	
	/*  Verlet list (the pairlist for non-bonding interaction)  */
	Int    max_nverlets;
	Int    nverlets;
	MDVerlet *verlets; /* Variable size; the current limit is max_nverlets */
	Vector *verlets_dr;  /*  The movement of atoms since the Verlet list was updated last  */
	Int    last_verlet_step;  /*  The timestep when the Verlet list was updated last */
	Int    *verlet_i;    /*  The verlet entries for atom i begins at verlet_i[i]  */
	
	/*  Save snapshot and restore  */
	Int    nsnapshots;
	MDSnapshot **snapshots;
	
	/*  For relocating center-of-mass  */
	Vector initial_center;
	
	/*  Minimize: conjugate gradient algorithm  */
	Vector *old_forces;		  /*  The forces in the previous step  */
	Vector *old_pos;          /*  = old_forces + natoms, The original positions before trial movements  */
	Double  f_len2;            /*  Sum of square lengths of force  */
	Double  old_f_len2;        /*  Old f_len2 */
	Double  v_len2;            /*  Sum of square lengths of "velocity" (= current searching direction)  */
	Double  max_gradient;      /*  maximum gradient  */
	Int    conv_flag;         /*  0: not converged, 1: converged by coordinate, 2: converged by gradient  */
	
	/*  Surface potential calculation  */
	struct SPArena *sp_arena;
//	Double  surface_potential; /*  The surface potential  */
//	Vector *sp_grad;          /*  The gradient of the surface potential for each atom coordinates */
	
	/*  Graphite calculation  */
	struct MDGraphiteArena *graphite;
	
} MDArena;

Double md_rand(void);
void md_srand(unsigned int seed);

int md_fit_coordinates(MDArena *arena, Int refno, Double *weights, Transform trans);

void md_scale_velocities(MDArena *arena);
void md_init_velocities(MDArena *arena);
/* int md_symmetry_relate(MDArena *arena, int n1, int n2, int sym_op); */
void md_transform_vec_by_symmetry(MDArena *arena, Vector *dst, const Vector *src, Symop rec, int no_transform);

void md_init_for_positions(MDArena *arena);
const char *md_prepare(MDArena *arena, int check_only);
int md_check_abnormal_bond(MDArena *arena, Molecule *mol, int idx);
int md_check_abnormal_angle(MDArena *arena, Molecule *mol, int idx);
int md_check_abnormal_dihedral(MDArena *arena, Molecule *mol, int idx);
int md_check_abnormal_improper(MDArena *arena, Molecule *mol, int idx);

void md_amend_by_symmetry(MDArena *arena);
void md_cell_recalculate(MDArena *arena);
void md_scale_cell(MDArena *arena, const Transform tf, int scale_atoms);
void md_wrap_coordinates(MDArena *arena);

void md_snapshot(MDArena *arena, int idx);
void md_restore(MDArena *arena, int idx);

void md_calc_kinetic_energy(MDArena *arena);

int md_output_results(MDArena *arena);

int md_copy_coordinates_from_internal(MDArena *arena);

int md_is_running(MDArena *arena);
	
int md_main(MDArena *arena, int minimize);
void md_set_default(MDArena *arena);

MDArena *md_arena_new(Molecule *mol);
MDArena *md_arena_set_molecule(MDArena *arena, Molecule *mol);
MDArena *md_arena_retain(MDArena *arena);
void md_arena_release(MDArena *arena);
void md_finish(MDArena *arena);

void md_panic(MDArena *arena, const char *fmt,...);
void md_debug(MDArena *arena, const char *fmt,...);
int md_log(MDArena *arena, const char *fmt, ...);
int md_warning(MDArena *arena, const char *fmt, ...);

#ifdef __cplusplus
}
#endif
		
#endif /* __MDCORE_H__ */
