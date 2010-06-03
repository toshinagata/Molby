/*
 *  MDForce.c
 *
 *  Created by Toshi Nagata on 2005/06/07.
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

#include "MDForce.h"
#include "MDGraphite.h"

#include <stdlib.h>
#include <string.h>

extern int do_custom_bond_callback(MDArena *arena, Double r, void *procObj, Double *energy, Double *force);

static void
s_custom_bond_force(MDArena *arena, Double r, MDCustomBondPar *cp, Double *energy, Double *force)
{
	float s;
	switch (cp->type) {
		case kMorseType:
			s = exp(-cp->u.morse.a * (r - cp->u.morse.r0));
			*energy = cp->u.morse.D * (1.0 + s * (s - 2.0));
			*force = -2.0 * cp->u.morse.a * cp->u.morse.D * s * (1.0 - s);
			break;
		case kTCLProcType:
		/*	retval = do_custom_bond_callback(arena, r, cp->u.proc, energy, force);
			if (retval != 0)
				md_panic(arena, NULL); */
			break;
		default:
			*energy = *force = 0.0;
			break;
	}
}

static void
s_calc_bond_force(MDArena *arena)
{
	Atom *ap = arena->mol->atoms;
	Int nbonds = arena->mol->nbonds;
	Int *bonds = arena->mol->bonds;
	Int *bond_par_i = arena->bond_par_i;
	Int *custom_bond_par_i = arena->custom_bond_par_i;
	BondPar *bond_pars = arena->par->bondPars;
	Double *anbond_r0 = arena->anbond_r0;
	Double *energies = &arena->energies[kBondIndex];
	Vector *forces = &arena->forces[kBondIndex * arena->mol->natoms];
	int i;
	for (i = 0; i < nbonds; i++, bonds += 2, bond_par_i++) {
		Vector r12;
		Double k0, k1, w1, w2, r0;
		Int idx;
		BondPar *bp;
		if (*bond_par_i < 0)
			continue;  /*  Ignore this entry  */
	/*	if (bond_uniq != NULL && bond_uniq[i] >= 0)
			continue;  *//*  Non-unique bond  */
		VecSub(r12, ap[bonds[1]].r, ap[bonds[0]].r);
		w1 = VecLength(r12);
		if (custom_bond_par_i != NULL && (idx = custom_bond_par_i[i]) > 0 && idx <= arena->ncustom_bond_pars) {
			Double energy, force;
			s_custom_bond_force(arena, w1, arena->custom_bond_pars + (idx - 1), &energy, &force);
			k0 = energy;
			k1 = force / w1;
			r0 = 0.0;  /*  To keep complier happy  */
		} else {
			bp = bond_pars + *bond_par_i;
			r0 = bp->r0;
			if (anbond_r0 != NULL) {
				if (anbond_r0[i] > 0.0)
					r0 += anbond_r0[i];
			}
			w2 = w1 - r0;
			k0 = bp->k * w2 * w2;         /*  Energy  */
			k1 = -2.0 * bp->k * w2 / w1;  /*  Force / r  */
		}
		VecScaleSelf(r12, k1);
		*energies += k0;
		VecDec(forces[bonds[0]], r12);
		VecInc(forces[bonds[1]], r12);
		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "bond force %d-%d: r=%f, r0=%f, k1=%f, {%f %f %f}\n", bonds[0]+1, bonds[1]+1, w1, r0, k1, r12.x, r12.y, r12.z);
		}
	}
}

static void
s_calc_angle_force(MDArena *arena)
{
	Atom *ap = arena->mol->atoms;
	Int nangles = arena->mol->nangles;
	Int *angles = arena->mol->angles;
	Int *angle_par_i = arena->angle_par_i;
	AnglePar *angle_pars = arena->par->anglePars;
/*	Int *angle_uniq = (arena->nsyms > 0 ? arena->sym_angle_uniq : NULL); */
	Double *energies = &arena->energies[kAngleIndex];
	Vector *forces = &arena->forces[kAngleIndex * arena->mol->natoms];
	int i;
	for (i = 0; i < nangles; i++, angles += 3, angle_par_i++) {
		Vector r21, r23, f1, f2, v1;
		Double k0, k1, w1, w2, w3;
		Double cost, sint, t;
		AnglePar *anp;
		if (*angle_par_i < 0)
			continue;  /*  Ignore this entry  */
	/*	if (angle_uniq != NULL && angle_uniq[i] >= 0)
			continue; */ /*  Non-unique angle  */
		anp = angle_pars + *angle_par_i;
		VecSub(r21, ap[angles[0]].r, ap[angles[1]].r);
		VecSub(r23, ap[angles[2]].r, ap[angles[1]].r);
		w1 = 1.0 / VecLength(r21);
		w2 = 1.0 / VecLength(r23);
		VecScaleSelf(r21, w1);
		VecScaleSelf(r23, w2);
		cost = VecDot(r21, r23);
		if (cost > 0.999999 || cost < -0.999999) {
		/*	printf("Cannot handle linear angle %d-%d-%d: skipped.\n", angles[0]+1, angles[1]+1, angles[2]+1); */
			continue;  /*  Cannot handle linear angle  */
		}
		sint = sqrt(1.0 - cost * cost);
	/*	if (sint < 1e-5 && sint > -1e-5)
			continue;  *//*  Cannot handle linear angle  */
		t = atan2(sint, cost) - anp->a0;
		k0 = anp->k * t * t;
	
		if (sint < 0.1 && sint > -0.1) {
			/* ---- Method 1 ---- */
			/* This is slower than method 2, but it is safer when sint is close to 0 */
			k1 = -2 * anp->k * t;
			VecCross(v1, r21, r23);
			VecCross(f1, r21, v1);
			w3 = w1 * k1 / VecLength(f1);
			if (!isfinite(w3))
				continue;
			VecScaleSelf(f1, w3);
			VecCross(f2, v1, r23);
			w3 = w2 * k1 / VecLength(f2);
			if (!isfinite(w3))
				continue;
			VecScaleSelf(f2, w3);
		} else {
			/* ---- Method 2 ---- */
			k1 = -2 * anp->k * t / sint;
			VecScale(f1, r21, cost);
			VecDec(f1, r23);
			w3 = w1 * k1;
			VecScaleSelf(f1, w3);
			VecScale(f2, r23, cost);
			VecDec(f2, r21);
			w3 = w2 * k1;
			VecScaleSelf(f2, w3);
		}

		*energies += k0;
		VecInc(forces[angles[0]], f1);
		VecInc(forces[angles[2]], f2);
		VecInc(f1, f2);
		VecDec(forces[angles[1]], f1);

		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "angle force %d-%d-%d: a=%f, a0=%f, k1=%f, {%f %f %f}, {%f %f %f}\n", angles[0]+1, angles[1]+1, angles[2]+1, (t+anp->a0)*180/PI, anp->a0*180/PI, k1, f1.x, f1.y, f1.z, f2.x, f2.y, f2.z);
		}

	}
}

static void
s_calc_dihedral_force_sub(MDArena *arena, Atom *ap, Int ndihedrals, Int *dihedrals, Int *dihedral_par_i, TorsionPar *dihedral_pars, Int *dihedral_uniq, Double *energies, Vector *forces)
{
	int i;
	if (arena->mol->nsyms == 0)
		dihedral_uniq = NULL;  /*  Ignore the symmetry info  */
	for (i = 0; i < ndihedrals; i++, dihedrals += 4, dihedral_par_i++) {
		Vector r21, r32, r43;
		Vector v1, v2, v3;
		Double w1, w2, w3, k0, k1;
		Double cosphi, sinphi, phi;
		Vector f1, f2, f3;
		TorsionPar *tp;
		int n;

		if (*dihedral_par_i < 0)
			continue;  /*  Ignore this entry  */
	/*	if (dihedral_uniq != NULL && dihedral_uniq[i] >= 0)
			continue;  *//*  Non-unique dihedral  */
		else tp = dihedral_pars + *dihedral_par_i;
		VecSub(r21, ap[dihedrals[0]].r, ap[dihedrals[1]].r);
		VecSub(r32, ap[dihedrals[1]].r, ap[dihedrals[2]].r);
		VecSub(r43, ap[dihedrals[2]].r, ap[dihedrals[3]].r);
		VecCross(v1, r21, r32);
		VecCross(v2, r32, r43);
		VecCross(v3, r32, v1);
		w1 = VecLength(v1);
		w2 = VecLength(v2);
		w3 = VecLength(v3);
		if (w1 < 1e-5 || w2 < 1e-5 || w3 < 1e-5)
			continue;  /*  The dihedral cannot be defined  */
		w1 = 1.0/w1;
		w2 = 1.0/w2;
		w3 = 1.0/w3;
		VecScaleSelf(v1, w1);
		VecScaleSelf(v2, w2);
		VecScaleSelf(v3, w3);
		
		/*  The dihedral angle value  */
		cosphi = VecDot(v1, v2);
		sinphi = VecDot(v3, v2);
		phi = -atan2(sinphi, cosphi);
		
		/*  Repeat for multiple dihedral terms  */
		k0 = k1 = 0.0;
		for (n = 0; n < tp->mult; n++) {
			Double k = tp->k[n];
			Double phi0 = tp->phi0[n];
			Int period = tp->period[n];
			if (period > 0) {
				k0 += k * (1 + cos(period * phi + phi0));
				k1 -= period * k * sin(period * phi + phi0);
			} else {
				/*  A simple quadratic term  */
				phi -= phi0;
				if (phi < -PI)
					phi += 2 * PI;
				else if (phi > PI)
					phi -= 2 * PI;
				k0 += k * phi * phi;
				k1 += 2 * k * phi;
			}
		}
		
		if (sinphi < -0.1 || sinphi > 0.1) {
			/*  The sin form  */
			Vector v4, v5, vw0;
			VecScale(v4, v1, cosphi);
			VecDec(v4, v2);
			VecScaleSelf(v4, w1);
			VecScale(v5, v2, cosphi);
			VecDec(v5, v1);
			VecScaleSelf(v5, w2);
			k1 /= sinphi;
			VecCross(f1, r32, v4);
			VecCross(f3, v5, r32);
			VecCross(f2, v4, r21);
			VecCross(vw0, r43, v5);
			VecInc(f2, vw0);
			VecScaleSelf(f1, k1);
			VecScaleSelf(f2, k1);
			VecScaleSelf(f3, k1);
		} else {
			/*  The cos form  */
			Vector v6, v7, vw1, vw2;
			VecScale(v6, v3, sinphi);
			VecDec(v6, v2);
			VecScaleSelf(v6, w3);
			VecScale(v7, v2, sinphi);
			VecDec(v7, v3);
			VecScaleSelf(v7, w2);
			k1 = -k1 / cosphi;
			VecCross(vw1, v6, r32);
			VecCross(f1, r32, vw1);
			VecCross(f3, v7, r32);
			VecCross(f2, r43, v7);
			VecCross(vw1, r32, r21);
			VecCross(vw2, v6, vw1);
			VecInc(f2, vw2);
			VecCross(vw1, r32, v6);
			VecCross(vw2, r21, vw1);
			VecInc(f2, vw2);
			VecScaleSelf(f1, k1);
			VecScaleSelf(f2, k1);
			VecScaleSelf(f3, k1);
		}
	/*	if (dihedral_uniq != NULL)
			k0 *= -dihedral_uniq[i]; */
		*energies += k0;
		VecInc(forces[dihedrals[0]], f1);
		VecDec(f1, f2);
		VecDec(forces[dihedrals[1]], f1);
		VecDec(f2, f3);
		VecDec(forces[dihedrals[2]], f2);
		VecDec(forces[dihedrals[3]], f3);

		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "dihedral(improper) force %d-%d-%d-%d: phi=%f, k1=%f, {%f %f %f}, {%f %f %f}, {%f %f %f}\n", dihedrals[0]+1, dihedrals[1]+1, dihedrals[2]+1, dihedrals[3]+1, phi*180/PI, k1, f1.x, f1.y, f1.z, f2.x, f2.y, f2.z, f3.x, f3.y, f3.z);
		}
	}
}

static void
s_calc_dihedral_force(MDArena *arena)
{
	Molecule *mol = arena->mol;
	Parameter *par = arena->par;
	Double *energies = &arena->energies[kDihedralIndex];
	Vector *forces = &arena->forces[kDihedralIndex * arena->mol->natoms];
	s_calc_dihedral_force_sub(arena, mol->atoms, mol->ndihedrals, mol->dihedrals, arena->dihedral_par_i, par->dihedralPars, NULL/*arena->sym_dihedral_uniq*/, energies, forces);
}

static void
s_calc_improper_force(MDArena *arena)
{
	Molecule *mol = arena->mol;
	Parameter *par = arena->par;
	Double *energies = &arena->energies[kImproperIndex];
	Vector *forces = &arena->forces[kImproperIndex * arena->mol->natoms];
	s_calc_dihedral_force_sub(arena, mol->atoms, mol->nimpropers, mol->impropers, arena->improper_par_i, par->improperPars, NULL/*arena->sym_improper_uniq*/, energies, forces);
}

/*  ==================================================================== */
/*  md_check_verlet_list: only for debugging the verlet list generation  */
/*  ==================================================================== */

typedef struct md_check_record {
	Int i, j, dx, dy, n;
	Double len;
} md_check_record;

static int
s_md_check_verlet_comparator(const void *ap, const void *bp)
{
	Double d = ((const md_check_record *)ap)->len - ((const md_check_record *)bp)->len;
	return (d < 0 ? -1 : (d > 0 ? 1 : 0));
}

void
md_check_verlet_list(MDArena *arena)
{
	int i, j, k;
	int dx, dy, dz, nn;
	int ndx, ndy, ndz;
	Atom *api, *apj;
	Vector cell_offsets[27];
	XtalCell *cell = arena->mol->cell;
	md_check_record *cr;
	int ncr;

	ndx = (arena->periodic_a ? 1 : 0);
	ndy = (arena->periodic_b ? 1 : 0);
	ndz = (arena->periodic_c ? 1 : 0);
	if (cell != NULL && (ndx != 0 || ndy != 0 || ndz != 0)) {
		nn = 0;
		for (i = -1; i <= 1; i++) {
			for (j = -1; j <= 1; j++) {
				for (k = -1; k <= 1; k++) {
					VecZero(cell_offsets[nn]);
					VecScaleInc(cell_offsets[nn], cell->axes[0], i);
					VecScaleInc(cell_offsets[nn], cell->axes[1], j);
					VecScaleInc(cell_offsets[nn], cell->axes[2], k);
					nn++;
				}
			}
		}
	}
	/*  Debugging is done with x-y 2-dimensional system  */
	ncr = arena->mol->natoms * (arena->mol->natoms - 1) / 2 * 9;
	cr = (md_check_record *)malloc(sizeof(md_check_record) * ncr);
	k = 0;
	for (i = 0, api = arena->mol->atoms; i < arena->mol->natoms; i++, api = ATOM_NEXT(api)) {
		for (j = i + 1, apj = ATOM_AT_INDEX(arena->mol->atoms, j); j < arena->mol->natoms; j++, apj = ATOM_NEXT(apj)) {
			Vector dr;
			VecSub(dr, apj->r, api->r);
			for (dx = -1; dx <= 1; dx++) {
				for (dy = -1; dy <= 1; dy++) {
					Vector dr2 = dr;
					nn = dx * 9 + dy * 3 + 13;
					VecInc(dr2, cell_offsets[nn]);
					cr[k].i = i;
					cr[k].j = j;
					cr[k].dx = dx;
					cr[k].dy = dy;
					cr[k].n = -1;
					cr[k].len = VecLength(dr2);
					k++;
				}
			}
		}
	}
	for (k = 0; k < arena->nverlets; k++) {
		i = arena->verlets[k].n1;
		j = arena->verlets[k].n2;
		if (i > j) {
			j = i;
			i = arena->verlets[k].n2;
		}
		dx = arena->verlets[k].symop.dx;
		dy = arena->verlets[k].symop.dy;
		nn = (i * arena->mol->natoms - i * (i + 1) / 2 + (j - i) - 1) * 9 + (dx + 1) * 3 + dy + 1;
		if (cr[nn].i != i || cr[nn].j != j || cr[nn].dx != dx || cr[nn].dy != dy)
			fprintf(arena->log_result, "Debug: internal inconsistency\n");
		cr[nn].n = k;
	}
		
	mergesort(cr, ncr, sizeof(md_check_record), s_md_check_verlet_comparator);
	for (i = 0; i < ncr; i++) {
		Double len2;
		len2 = -1;
		if (cr[i].n >= 0)
			len2 = arena->verlets[cr[i].n].length;
		fprintf(arena->log_result, "Debug: %5d (i=%4d,j=%4d) [dx=%4d, dy=%4d] n=%4d %f %f\n", i, cr[i].i, cr[i].j, cr[i].dx, cr[i].dy, cr[i].n, cr[i].len, len2);
		if (cr[i].len > arena->pairlist_distance)
			break;
	}
	while (i < ncr) {
		if (cr[i].n != -1) {
			fprintf(arena->log_result, "Debug: %5d (i=%4d,j=%4d) [dx=%4d, dy=%4d] n=%4d %f %f\n", i, cr[i].i, cr[i].j, cr[i].dx, cr[i].dy, cr[i].n, cr[i].len, arena->verlets[cr[i].n].length);
		}
		i++;
	}
	fflush(arena->log_result);
	free(cr);
}

/*  ==================================================================== */


/*  Update the Verlet list (for non-bonding interaction)  */
static void
s_make_verlet_list(MDArena *arena)
{
	int i, j, k, n, nn, natoms;
	int dx, dy, dz, ndx, ndy, ndz;
	Molecule *mol = arena->mol;
	Atom *atoms = mol->atoms;
	Atom *api, *apj;
	MDVerlet *vl = arena->verlets;
	Vector *vdr = arena->verlets_dr;
	Double limit;
	Byte use_sym;
	XtalCell *cell = mol->cell;
	Vector cell_offsets[27];
	Parameter *par = arena->par;
	Double vdw_cutoff;
	
	natoms = mol->natoms;
	for (i = 0; i < natoms; i++)
		VecZero(vdr[i]);

	ndx = (arena->periodic_a ? 1 : 0);
	ndy = (arena->periodic_b ? 1 : 0);
	ndz = (arena->periodic_c ? 1 : 0);

	if (cell != NULL && (ndx != 0 || ndy != 0 || ndz != 0)) {
		nn = 0;
		for (i = -1; i <= 1; i++) {
			for (j = -1; j <= 1; j++) {
				for (k = -1; k <= 1; k++) {
					VecZero(cell_offsets[nn]);
					VecScaleInc(cell_offsets[nn], cell->axes[0], i);
					VecScaleInc(cell_offsets[nn], cell->axes[1], j);
					VecScaleInc(cell_offsets[nn], cell->axes[2], k);
					nn++;
				}
			}
		}
	}
	
	limit = arena->pairlist_distance * arena->pairlist_distance;
	n = 0;
	use_sym = (arena->mol->nsyms > 0);
	for (i = 0, api = atoms; i < natoms; i++, api++) {
		MDExclusion *exinfo = arena->exinfo + i;
		Int index4 = (exinfo + 1)->index0;
		Int vdw_idx1 = arena->vdw_par_i[i];
		Int vdw_idx, vdw_idx2;
		arena->verlet_i[i] = n;

		for (j = i + 1, apj = atoms + j; j < natoms; j++, apj++) {
			Vector rij;
			Double lenij2;
			Int *ip, index0;
			int exflag = 0;
			int mult = 1;
			int dxbase, dybase, dzbase;

			/*  Fixed atoms  */
			if (api->fix_force < 0 && apj->fix_force < 0)
				continue;

			VecSub(rij, apj->r, api->r);

			/*  Calculate the cell offset for the nearest neighbor  */
			/*  NOTE: the offset is calculated independently for each axis. This may result
			    in unexpected choice when the angles between the axes are far from 90 deg */
			if (apj->periodic_exclude == 0 && cell != NULL) {
				TransformPtr rtp = cell->rtr;
				TransformPtr tp = cell->tr;
				Double w;
				if (arena->periodic_a) {
					w = rtp[0] * rij.x + rtp[1] * rij.y + rtp[2] * rij.z;
					dxbase = floor(-w + 0.5);
				} else dxbase = 0;
				if (arena->periodic_b) {
					w = rtp[3] * rij.x + rtp[4] * rij.y + rtp[5] * rij.z;
					dybase = floor(-w + 0.5);
				} else dybase = 0;
				if (arena->periodic_c) {
					w = rtp[6] * rij.x + rtp[7] * rij.y + rtp[8] * rij.z;
					dzbase = floor(-w + 0.5);
				} else dzbase = 0;
				rij.x += tp[0] * dxbase + tp[1] * dybase + tp[2] * dzbase;
				rij.y += tp[3] * dxbase + tp[4] * dybase + tp[5] * dzbase;
				rij.z += tp[6] * dxbase + tp[7] * dybase + tp[8] * dzbase;
			} else dxbase = dybase = dzbase = 0;

			/*  Non unique atom pair  */
		/*	if (use_sym && api->symop.alive && apj->symop.alive)
				continue; */

			/*  Check the specific cutoff table for (i, j) pair  */
			vdw_idx2 = arena->vdw_par_i[j];
			if (vdw_idx1 < 0 || vdw_idx2 < 0)
				vdw_idx = par->nvdwPars * par->nvdwPars;  /*  A null record  */
			else if (vdw_idx1 < vdw_idx2)
				vdw_idx = vdw_idx1 * par->nvdwPars + vdw_idx2;
			else
				vdw_idx = vdw_idx2 * par->nvdwPars + vdw_idx1;
			vdw_cutoff = arena->cutoff;
			for (k = par->nvdwCutoffPars - 1; k >= 0; k--) {
				VdwCutoffPar *cr = par->vdwCutoffPars + k;
				if (cr->type == 0) {
					if (((cr->n1 < 0 || cr->n1 == api->type) && (cr->n2 < 0 || cr->n2 == apj->type)) ||
						((cr->n1 < 0 || cr->n1 == apj->type) && (cr->n2 < 0 || cr->n2 == api->type))) {
						vdw_cutoff = cr->cutoff;
						break;
					}
				} else {
					if ((cr->n1 <= i && i <= cr->n2 && cr->n3 <= j && j <= cr->n4)||
						(cr->n3 <= i && i <= cr->n4 && cr->n1 <= j && j <= cr->n2)) {
						vdw_cutoff = cr->cutoff;
						break;
					}
				}
			}
			if (vdw_cutoff < 0) {
				/*  A negative value of specific cutoff means "cutoff at r_eq * (-specific_cutoff)  */
				MDVdwCache *cp = &(arena->vdw_cache[vdw_idx]);
				Double r_eq = pow(cp->par.A / cp->par.B * 2.0, 1.0/6.0);
				vdw_cutoff = r_eq * (-vdw_cutoff);
			}

			/*  Search for pairs, taking periodicity into account  */
			for (dx = -ndx; dx <= ndx; dx++) {
				for (dy = -ndy; dy <= ndy; dy++) {
					for (dz = -ndz; dz <= ndz; dz++) {
						Vector rij0 = rij;
						nn = dx * 9 + dy * 3 + dz + 13; 
						if (nn == 13) {
							/*  Pair within the unit cell  */
							/*  Is this pair to be excluded?  */
							for (index0 = exinfo->index0, ip = arena->exlist + index0; index0 < index4; index0++, ip++) {
								if (*ip == j)
									break;
							}
							if (index0 < exinfo->index3)
								continue;  /*  Special exclusion, 1-2, 1-3  */
							if (index0 < index4)
								exflag = 1;  /*  1-4 interaction  */
						} else if (apj->periodic_exclude) {
							continue;
						} else {
							VecInc(rij0, cell_offsets[nn]);
							exflag = 0;
						}

						lenij2 = VecLength2(rij0);
						if (lenij2 <= limit) {
							MDVerlet *vlp;
							if (n >= arena->max_nverlets) {
								arena->max_nverlets += 32;
								vl = (MDVerlet *)realloc(vl, sizeof(MDVerlet) * arena->max_nverlets);
								if (vl == NULL)
									md_panic(arena, "Low memory");
							}
							vlp = &vl[n];
							vlp->vdw_type = (exflag ? 1 : 0);
							vlp->mult = mult;
							vlp->symop.dx = dx + dxbase;
							vlp->symop.dy = dy + dybase;
							vlp->symop.dz = dz + dzbase;
							vlp->symop.sym = 0;
							vlp->symop.alive = (vlp->symop.dx != 0 || vlp->symop.dy != 0 || vlp->symop.dz != 0);
							vlp->index = vdw_idx;
							vlp->n1 = i;
							vlp->n2 = j;
							vlp->vdw_cutoff = vdw_cutoff;
							vlp->length = sqrt(lenij2);
							n++;
						} /* end if lenij2 <= limit */
					} /* end loop dz */
				} /* end loop dy */
			} /* end loop dx */
		} /* end loop j */

	} /* end loop i */
	arena->verlet_i[natoms] = n;
	arena->nverlets = n;
	arena->verlets = vl;
	arena->last_verlet_step = arena->step;
		
	if (arena->debug_result && arena->debug_output_level > 2) {
		fprintf(arena->debug_result, "\n  Verlet list at step %d\n", arena->step);
		fprintf(arena->debug_result, "  {atom1 atom2 vdw_type (=1 if 1-4 bonded) vdw_index}\n");
		for (i = 0; i < arena->nverlets; i++) {
			fprintf(arena->debug_result, "{%d %d %d %d}%c", vl[i].n1+1, vl[i].n2+1, vl[i].vdw_type, vl[i].index, (i % 4 == 3 ? '\n' : ' '));
		}
		if (i % 4 != 0)
			fprintf(arena->debug_result, "\n");
	}
}

/*  Calculate the nonbonded force  */
/*  group_flags[] is used for calculation of pair-specific interaction between
    two groups of atoms  */
static void
s_calc_nonbonded_force_sub(MDArena *arena, Double *energies, Double *eenergies, Vector *forces, Vector *eforces, const MDGroupFlags *group_flags_1, const MDGroupFlags *group_flags_2)
{
	int i;
	Vector *vdr;
	Double limit, elimit, dielec_r;
	MDVerlet *vl;
/*	MDVdwCache *vdw_cache = arena->vdw_cache; */
	Atom *atoms = arena->mol->atoms;
	XtalCell *cell = arena->mol->cell;

	/*  Check whether the Verlet list needs update  */
	if (arena->last_verlet_step >= 0) {
		i = arena->mol->natoms - 1;
		vdr = arena->verlets_dr + i;
		if (arena->cutoff > arena->electro_cutoff)
			limit = (arena->pairlist_distance - arena->cutoff) * 0.5;
		else
			limit = (arena->pairlist_distance - arena->electro_cutoff) * 0.5;
		limit = limit * limit;
		for ( ; i >= 0; i--, vdr--) {
			if (VecLength2(*vdr) >= limit)
				break;
		}
	} else i = 0;
	if (i >= 0)
		s_make_verlet_list(arena);
	
	/*  Calculate the non-bonded interaction for each pair in the Verlet list  */
/*	limit = arena->cutoff * arena->cutoff; */
	elimit = arena->electro_cutoff * arena->electro_cutoff;
	dielec_r = COULOMBIC / arena->dielectric;
	for (i = arena->nverlets - 1, vl = arena->verlets + i; i >= 0; i--, vl--) {
		Double A, B, vofs, k0, k1;
		MDVdwCache *vp = arena->vdw_cache + vl->index;
		Atom *ap1, *ap2;
		Vector rij, fij;
		Double r2, w2, w6, w12;
		if (group_flags_1 != NULL && group_flags_2 != NULL) {
			if (!(get_group_flag(group_flags_1, vl->n1) && get_group_flag(group_flags_2, vl->n2))
			&& !(get_group_flag(group_flags_1, vl->n2) && get_group_flag(group_flags_2, vl->n1)))
				continue;
		}
		if (vl->vdw_type == 1) {
			A = vp->par.A14;
			B = vp->par.B14;
		/*	vofs = vp->vcut14; */
		} else {
			A = vp->par.A;
			B = vp->par.B;
		/*	vofs = vp->vcut; */
		}
		ap1 = &atoms[vl->n1];
		ap2 = &atoms[vl->n2];
		rij = ap2->r;
		if (vl->symop.alive) {
			VecScaleInc(rij, cell->axes[0], vl->symop.dx);
			VecScaleInc(rij, cell->axes[1], vl->symop.dy);
			VecScaleInc(rij, cell->axes[2], vl->symop.dz);
		}
		VecDec(rij, ap1->r);
		r2 = VecLength2(rij);
		if (r2 >= elimit)
			continue;
		limit = vl->vdw_cutoff * vl->vdw_cutoff;
		if (r2 >= limit)
			continue;
		
		fij.x = fij.y = fij.z = 0.0;
		
		/*  Coulombic force  */
		w12 = ap1->charge * ap2->charge * dielec_r;
		if (w12 != 0.0) {
			if (arena->use_xplor_shift) {
				w2 = 1.0 / sqrt(r2);
				k0 = r2 / elimit - 1.0;
				k0 = k0 * k0;
				k0 = w12 * k0 * w2;
				k1 = (3.0 * r2 / elimit - 2.0) / elimit - 1.0 / r2;
				k1 = -w12 * k1 * w2;
			} else {
				w2 = 1.0 / sqrt(r2);
				w6 = w2 / r2;
				k0 = w12 * w2;
				k1 = w12 * w6;
				vofs = -w12 / arena->electro_cutoff;
				k0 += vofs;
			}
			if (vl->vdw_type == 1) {
				k0 *= arena->scale14_elect;
				k1 *= arena->scale14_elect;
			}
			VecScale(fij, rij, k1);
			*eenergies += k0 / vl->mult;
			if (eforces != NULL) {
				VecDec(eforces[vl->n1], fij);
				VecInc(eforces[vl->n2], fij);
			}
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "nonbonded(electrostatic) force %d-%d: r=%f, k0=%f, k1=%f, {%f %f %f}\n", vl->n1+1, vl->n2+1, sqrt(r2), k0/KCAL2INTERNAL, k1*sqrt(r2)/KCAL2INTERNAL, fij.x/KCAL2INTERNAL, fij.y/KCAL2INTERNAL, fij.z/KCAL2INTERNAL);
			}
		}

		/*  van der Waals force  */
		w2 = 1.0 / r2;
		w6 = w2 * w2 * w2;
		w12 = w6 * w6;
		k0 = A * w12 - B * w6;
		k1 = 6 * w2 * (2.0 * A * w12 - B * w6);
		w2 = 1.0 / limit;
		w6 = w2 * w2 * w2;
		w12 = w6 * w6;
		vofs = -A * w12 + B * w6;
		k0 += vofs;
		if (vl->vdw_type == 1) {
			k0 *= arena->scale14_vdw;
			k1 *= arena->scale14_vdw;
		}
		VecScale(fij, rij, k1);
		*energies += k0 / vl->mult;
		if (forces != NULL) {
			VecDec(forces[vl->n1], fij);
			VecInc(forces[vl->n2], fij);
		}
		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "nonbonded(vdw) force %d-%d: r=%f, k0=%f, k1=%f, {%f %f %f}\n", vl->n1+1, vl->n2+1, sqrt(r2), k0/KCAL2INTERNAL, k1*sqrt(r2)/KCAL2INTERNAL, fij.x/KCAL2INTERNAL, fij.y/KCAL2INTERNAL, fij.z/KCAL2INTERNAL);
		}
	}
}

static void
s_calc_nonbonded_force(MDArena *arena)
{
	Double *energies = &arena->energies[kVDWIndex];
	Double *eenergies = &arena->energies[kElectrostaticIndex];
	Vector *forces = &arena->forces[kVDWIndex * arena->mol->natoms];
	Vector *eforces = &arena->forces[kElectrostaticIndex * arena->mol->natoms];
	s_calc_nonbonded_force_sub(arena, energies, eenergies, forces, eforces, NULL, NULL);
}

static void
s_calc_auxiliary_force(MDArena *arena)
{
	Double *energies = &arena->energies[kAuxiliaryIndex];
	Vector *forces = &arena->forces[kAuxiliaryIndex * arena->mol->natoms];
	Atom *atoms = arena->mol->atoms;
	int natoms = arena->mol->natoms;
	int i, j;

	/*  Centric force - OBSOLETE */
/*	if (arena->centric_potential_force > 0) {
		Vector center = arena->mol->atoms[arena->centric_potential_center].r;
		Double rin = arena->centric_potential_inner_limit;
		Double rout = arena->centric_potential_outer_limit;
		Double k = arena->centric_potential_force;
		for (i = 0; i < natoms; i++) {
			Vector r21;
			Double w1, w2, k0, k1;
			VecSub(r21, atoms[i].r, center);
			w2 = VecLength2(r21);
			w1 = sqrt(w2);
			if (w1 <= rin) {
				k0 = k1 = 0.0;
			} else if (rout > 0 && w1 > rout) {
				k0 = k * (rout - rin) * (w1 + w1 - rout + rin);
				k1 = -2.0 * k * (rout - rin) / w1;
			} else {
				k0 = k * (w1 - rin) * (w1 - rin);
				k1 = -2.0 * k * (1 - rin / w1);
			}
			*energies += k0;
			VecScaleSelf(r21, k1);
			VecInc(forces[i], r21);
#if DEBUG
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "auxiliary(centric) force %d: r=%f, k1=%f, {%f %f %f}\n", i, w1, k1*w1, r21.x, r21.y, r21.z);
			}
#endif
		}
	}
	*/

	/*  Spherical boundary conditions  */
	if (arena->spherical_bc_force > 0) {
		Vector center = arena->spherical_bc_center;
		Double rin = arena->spherical_bc_inner_limit;
		Double rout = arena->spherical_bc_outer_limit;
		Double k = arena->spherical_bc_force;
		for (i = 0; i < natoms; i++) {
			Vector r21;
			Double w1, w2, k0, k1;
			VecSub(r21, atoms[i].r, center);
			w2 = VecLength2(r21);
			w1 = sqrt(w2);
			if (w1 <= rin) {
				k0 = k1 = 0.0;
			} else if (rout > 0 && w1 > rout) {
				k0 = k * (rout - rin) * (w1 + w1 - rout + rin);
				k1 = -2.0 * k * (rout - rin) / w1;
			} else {
				k0 = k * (w1 - rin) * (w1 - rin);
				k1 = -2.0 * k * (1 - rin / w1);
			}
			*energies += k0;
			VecScaleSelf(r21, k1);
			VecInc(forces[i], r21);
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "auxiliary(spherical BC) force %d: r=%f, k1=%f, {%f %f %f}\n", i, w1, k1*w1, r21.x, r21.y, r21.z);
			}
		}
	}

	/*  Box force  */
	if (arena->box_potential_force > 0) {
		Double xsize = arena->box_potential_xsize;
		Double ysize = arena->box_potential_ysize;
		Double zsize = arena->box_potential_zsize;
		Double k = arena->box_potential_force;
		for (i = 0; i < natoms; i++) {
			Vector r = atoms[i].r;
			if (r.x > xsize)
				r.x -= xsize;
			else if (r.x < -xsize)
				r.x += xsize;
			else r.x = 0.0;
			if (r.y > ysize)
				r.y -= ysize;
			else if (r.y < -ysize)
				r.y += ysize;
			else r.y = 0.0;
			if (r.z > zsize)
				r.z -= zsize;
			else if (r.z < -zsize)
				r.z += zsize;
			else r.z = 0.0;
			*energies += k * (r.x * r.x + r.y * r.y + r.z * r.z);
			VecScaleInc(forces[i], r, -2);
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "auxiliary(box) force %d: {%f %f %f}\n", i, -2*r.x, -2*r.y, -2*r.z);
			}
		}
	}

	/*  Fix atoms  */
	for (i = j = 0; i < natoms; i++) {
		Atom *ap = &atoms[i];
		Vector r21;
		Double w1, w2, k0, k1;
		if (ap->fix_force <= 0.0)
			continue;
		VecSub(r21, ap->r, ap->fix_pos);
		w2 = VecLength2(r21);
		w1 = sqrt(w2);
		k0 = ap->fix_force * w2;
		k1 = -2.0 * ap->fix_force * w1;
		*energies += k0;
		VecScaleSelf(r21, k1);
		VecInc(forces[i], r21);
		j++;
		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "auxiliary(fix) force %d: r=%f, k1=%f, {%f %f %f}\n", i, w1, k1*w1, r21.x, r21.y, r21.z);
		}
	}
	
	/*  Graphite  */
	if (arena->graphite != NULL && arena->use_graphite) {
		graphite_force(arena->graphite, arena, energies, forces);
	}
}

static void
s_md_debug_output_positions(MDArena *arena)
{
	int i;
	if (arena->debug_result == NULL || arena->debug_output_level == 0)
		return;
	fprintf(arena->debug_result, "\n  Atom positions, velocities, and forces at step %d\n", arena->step);
	fprintf(arena->debug_result, "%5s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n", "No.", "x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz");
	for (i = 0; i < arena->mol->natoms; i++) {
		Atom *ap = &arena->mol->atoms[i];
		fprintf(arena->debug_result, "%5d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", i+1, ap->r.x, ap->r.y, ap->r.z, ap->v.x, ap->v.y, ap->v.z, ap->f.x, ap->f.y, ap->f.z);
	}
}

void
calc_force(MDArena *arena)
{
	Int i, j, natoms;
	Vector *ff, *fa;
	Molecule *mol = arena->mol;
	Int doSurface = 0;

	natoms = mol->natoms;

	/*  Clear the energy/force storage  */
	for (i = 0; i < natoms; i++)
		VecZero(mol->atoms[i].f);
	
	memset(arena->energies, 0, sizeof(Double) * kSlowIndex);
	memset(arena->forces, 0, sizeof(Vector) * kSlowIndex * natoms);
	arena->total_energy = 0.0;
	
	if (arena->step == arena->start_step || arena->step % arena->surface_potential_freq == 0) {
		doSurface = 1;
		arena->energies[kSurfaceIndex] = 0.0;
		memset(arena->forces + kSurfaceIndex * natoms, 0, sizeof(Vector) * natoms);
	}

	s_calc_bond_force(arena);
	s_calc_angle_force(arena);
	s_calc_dihedral_force(arena);
	s_calc_improper_force(arena);
	s_calc_nonbonded_force(arena);
	s_calc_auxiliary_force(arena);

	if (doSurface && arena->probe_radius > 0.0) {
	/*	if (arena->sphere_points == 0)
			calc_surface_force(arena);
		else 
			calc_surface_force_2(arena); */
		calc_surface_force(arena);
	}
	
	/*  Sum up all partial forces and energies  */
	arena->total_energy = 0.0;
	for (i = 0; i < kKineticIndex; i++)
		arena->total_energy += arena->energies[i];

	for (i = 0; i < natoms; i++) {
		fa = &mol->atoms[i].f;
		ff = &arena->forces[i];
		for (j = 0; j < kKineticIndex; j++) {
			VecInc(*fa, *ff);
			ff += natoms;
		}
	}
	
	/*  The total (internal) force must be zero  */
/*	{
		Vector f;
		Double w;
		VecZero(f);
		for (i = 0; i < natoms; i++) {
			Vector *fp = &(mol->atoms[i].f);
			VecInc(f, *fp);
		}
		w = (natoms > 0 ? 1.0 / natoms : 0.0);
		VecScaleSelf(f, w);
		for (i = 0; i < natoms; i++) {
			Vector *fp = &(mol->atoms[i].f);
			VecDec(*fp, f);
		}
	} */
	
	s_md_debug_output_positions(arena);
}

void
calc_pair_interaction(MDArena *arena, const MDGroupFlags *group_flags_1, const MDGroupFlags *group_flags_2)
{
	Vector *forces, *eforces;
	if (arena->pair_forces != NULL) {
		forces = arena->pair_forces;
		eforces = forces + arena->mol->natoms;
		memset(forces, 0, sizeof(Vector) * arena->mol->natoms * 2);
	} else forces = eforces = NULL;
	arena->pair_energies[0] = arena->pair_energies[1] = 0.0;
	s_calc_nonbonded_force_sub(arena, &arena->pair_energies[0], &arena->pair_energies[1], forces, eforces, group_flags_1, group_flags_2);
}
