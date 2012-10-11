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
		if (ap[bonds[0]].mm_exclude || ap[bonds[1]].mm_exclude)
			continue;  /*  Skip non-occupied atoms  */
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

static inline int
s_calc_angle_force_one(const Vector *r1, const Vector *r2, const Vector *r3, const AnglePar *anp, Double *angle, Double *energy, Double *k1out, Vector *force1, Vector *force2, Vector *force3)
{
	Vector r21, r23, f1, f2, v1;
	Double w1, w2, w3, cost, sint, t, k0, k1;
	VecSub(r21, *r1, *r2);
	VecSub(r23, *r3, *r2);
	w1 = 1.0 / VecLength(r21);
	w2 = 1.0 / VecLength(r23);
	VecScaleSelf(r21, w1);
	VecScaleSelf(r23, w2);
	cost = VecDot(r21, r23);
	if (cost > 0.999999 || cost < -0.999999) {
		/*	printf("Cannot handle linear angle %d-%d-%d: skipped.\n", angles[0]+1, angles[1]+1, angles[2]+1); */
		return -1;  /*  Cannot handle linear angle  */
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
			return -1;
		VecScaleSelf(f1, w3);
		VecCross(f2, v1, r23);
		w3 = w2 * k1 / VecLength(f2);
		if (!isfinite(w3))
			return -1;
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
	
	*energy = k0;
	*force1 = f1;
	*force3 = f2;
	*angle = t + anp->a0;
	*k1out = k1;
	VecInc(f1, f2);
	VecScale(*force2, f1, -1);
	return 0;
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
		Atom *ap1, *ap2, *ap3;
		AnglePar *anp;
		Double en, ang, k1;
		Vector f1, f2, f3;
		if (*angle_par_i < 0)
			continue;  /*  Ignore this entry  */
		ap1 = ATOM_AT_INDEX(ap, angles[0]);
		ap2 = ATOM_AT_INDEX(ap, angles[1]);
		ap3 = ATOM_AT_INDEX(ap, angles[2]);
	/*	if (angle_uniq != NULL && angle_uniq[i] >= 0)
			continue; */ /*  Non-unique angle  */
		if (ap1->mm_exclude || ap2->mm_exclude || ap3->mm_exclude)
			continue;  /*  Skip non-occupied atoms  */
		if (arena->nalchem_flags > 0) {
			if (angles[0] < arena->nalchem_flags && angles[1] < arena->nalchem_flags
				&& angles[2] < arena->nalchem_flags
				&& (arena->alchem_flags[angles[0]] | arena->alchem_flags[angles[1]] | arena->alchem_flags[angles[2]]) == 3)
				continue;  /*  Interaction between vanishing and appearing groups is ignored  */
		}
		anp = angle_pars + *angle_par_i;
		if (s_calc_angle_force_one(&(ap1->r), &(ap2->r), &(ap3->r), angle_pars + *angle_par_i, &ang, &en, &k1, &f1, &f2, &f3) != 0)
			continue;
		*energies += en;
		VecInc(forces[angles[0]], f1);
		VecInc(forces[angles[1]], f2);
		VecInc(forces[angles[2]], f3);
		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "angle force %d-%d-%d: a=%f, a0=%f, k1=%f, {%f %f %f}, {%f %f %f}\n", angles[0]+1, angles[1]+1, angles[2]+1, ang*180/PI, anp->a0*180/PI, k1, -f2.x, -f2.y, -f2.z, f3.x, f3.y, f3.z);
		}
	}
}

static inline int
s_calc_dihedral_force_one(const Vector *r1, const Vector *r2, const Vector *r3, const Vector *r4, const TorsionPar *tp, Double *phi, Double *energy, Double *k1out, Vector *force1, Vector *force2, Vector *force3, Vector *force4)
{
	Vector r21, r32, r43;
	Vector v1, v2, v3;
	Double w1, w2, w3, k0, k1;
	Double cosphi, sinphi;
	Vector f1, f2, f3;
	Int n;
	VecSub(r21, *r1, *r2);
	VecSub(r32, *r2, *r3);
	VecSub(r43, *r3, *r4);
	VecCross(v1, r21, r32);
	VecCross(v2, r32, r43);
	VecCross(v3, r32, v1);
	w1 = VecLength(v1);
	w2 = VecLength(v2);
	w3 = VecLength(v3);
	if (w1 < 1e-5 || w2 < 1e-5 || w3 < 1e-5)
		return -1;  /*  The dihedral cannot be defined  */
	w1 = 1.0/w1;
	w2 = 1.0/w2;
	w3 = 1.0/w3;
	VecScaleSelf(v1, w1);
	VecScaleSelf(v2, w2);
	VecScaleSelf(v3, w3);
	
	/*  The dihedral angle value  */
	cosphi = VecDot(v1, v2);
	sinphi = VecDot(v3, v2);
	*phi = -atan2(sinphi, cosphi);
	
	/*  Repeat for multiple dihedral terms  */
	k0 = k1 = 0.0;
	for (n = 0; n < tp->mult; n++) {
		Double k = tp->k[n];
		Double phi0 = tp->phi0[n];
		Int period = tp->period[n];
		if (period > 0) {
			k0 += k * (1 + cos(period * *phi + phi0));
			k1 -= period * k * sin(period * *phi + phi0);
		} else {
			/*  A simple quadratic term  */
			phi0 = *phi - phi0;
			if (phi0 < -PI)
				phi0 += 2 * PI;
			else if (phi0 > PI)
				phi0 -= 2 * PI;
			k0 += k * phi0 * phi0;
			k1 += 2 * k * phi0;
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
	*energy = k0;
	*k1out = k1;
	*force1 = f1;
	VecDec(f1, f2);
	VecScale(*force2, f1, -1);
	VecDec(f2, f3);
	VecScale(*force3, f2, -1);
	VecScale(*force4, f3, -1);
	return 0;
}

static void
s_calc_dihedral_force_sub(MDArena *arena, Atom *ap, Int ndihedrals, Int *dihedrals, Int *dihedral_par_i, TorsionPar *dihedral_pars, Int *dihedral_uniq, Double *energies, Vector *forces)
{
	int i;
	if (arena->mol->nsyms == 0)
		dihedral_uniq = NULL;  /*  Ignore the symmetry info  */
	for (i = 0; i < ndihedrals; i++, dihedrals += 4, dihedral_par_i++) {
		Atom *ap1, *ap2, *ap3, *ap4;
		Double phi, en, k1;
		Vector f1, f2, f3, f4;
		TorsionPar *tp;

		if (*dihedral_par_i < 0)
			continue;  /*  Ignore this entry  */
	/*	if (dihedral_uniq != NULL && dihedral_uniq[i] >= 0)
			continue;  *//*  Non-unique dihedral  */
		ap1 = ATOM_AT_INDEX(ap, dihedrals[0]);
		ap2 = ATOM_AT_INDEX(ap, dihedrals[1]);
		ap3 = ATOM_AT_INDEX(ap, dihedrals[2]);
		ap4 = ATOM_AT_INDEX(ap, dihedrals[3]);
		if (ap1->mm_exclude || ap2->mm_exclude || ap3->mm_exclude || ap4->mm_exclude)
			continue;  /*  Skip non-occupied atoms  */
		if (arena->nalchem_flags > 0) {
			if (dihedrals[0] < arena->nalchem_flags && dihedrals[1] < arena->nalchem_flags
				&& dihedrals[2] < arena->nalchem_flags && dihedrals[3] < arena->nalchem_flags
				&& (arena->alchem_flags[dihedrals[0]] | arena->alchem_flags[dihedrals[1]]
					| arena->alchem_flags[dihedrals[2]] | arena->alchem_flags[dihedrals[3]]) == 3)
				continue;  /*  Interaction between vanishing and appearing groups is ignored  */
		}
		tp = dihedral_pars + *dihedral_par_i;
		if (s_calc_dihedral_force_one(&(ap1->r), &(ap2->r), &(ap3->r), &(ap4->r), tp, &phi, &en, &k1, &f1, &f2, &f3, &f4) != 0)
			continue;
		
		/*
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
			continue;  //  The dihedral cannot be defined 
		w1 = 1.0/w1;
		w2 = 1.0/w2;
		w3 = 1.0/w3;
		VecScaleSelf(v1, w1);
		VecScaleSelf(v2, w2);
		VecScaleSelf(v3, w3);
		
		//  The dihedral angle value
		cosphi = VecDot(v1, v2);
		sinphi = VecDot(v3, v2);
		phi = -atan2(sinphi, cosphi);
		
		//  Repeat for multiple dihedral terms
		k0 = k1 = 0.0;
		for (n = 0; n < tp->mult; n++) {
			Double k = tp->k[n];
			Double phi0 = tp->phi0[n];
			Int period = tp->period[n];
			if (period > 0) {
				k0 += k * (1 + cos(period * phi + phi0));
				k1 -= period * k * sin(period * phi + phi0);
			} else {
				//  A simple quadratic term
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
			//  The sin form
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
			//  The cos form
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
		*energies += k0;
		VecInc(forces[dihedrals[0]], f1);
		VecDec(f1, f2);
		VecDec(forces[dihedrals[1]], f1);
		VecDec(f2, f3);
		VecDec(forces[dihedrals[2]], f2);
		VecDec(forces[dihedrals[3]], f3);
*/
		*energies += en;
		VecInc(forces[dihedrals[0]], f1);
		VecInc(forces[dihedrals[1]], f2);
		VecInc(forces[dihedrals[2]], f3);
		VecInc(forces[dihedrals[3]], f4);
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

static void
s_calc_pibond_force(MDArena *arena)
{
	Atom *ap = arena->mol->atoms;
	Int npibonds = arena->mol->npibonds;
	Int *pibonds = arena->mol->pibonds;
	UnionPar *pars = arena->pi_pars;
	Double energy, val, valk1;
//	Double *energies = &arena->energies[kAngleIndex];
//	Vector *forces = &arena->forces[kAngleIndex * arena->mol->natoms];
	Int i;
	for (i = 0; i < npibonds; i++, pibonds += 4, pars++) {
		Vector r[4];
		Vector f[4];
		PiAtom *pp[4];
		Int j, n;
		for (j = 0; j < 4; j++) {
			n = pibonds[j];
			if (n < 0)
				break;
			else if (n < ATOMS_MAX_NUMBER) {
				r[j] = ATOM_AT_INDEX(ap, n)->r;
				pp[j] = NULL;
			} else {
				n -= ATOMS_MAX_NUMBER;
				MoleculeCalculatePiAtomPosition(arena->mol, n, &r[j]);
				pp[j] = arena->mol->piatoms + n;
			}
		}
		if (j == 2) {  /*  Bond  */
			Vector r12;
			Double w1, w2, k0, k1;
			VecSub(r12, r[1], r[0]);
			w1 = VecLength(r12);
			w2 = w1 - pars->bond.r0;
			k0 = pars->bond.k * w2 * w2;         /*  Energy  */
			k1 = -2.0 * pars->bond.k * w2 / w1;  /*  Force / r  */
			VecScaleSelf(r12, k1);
			energy = k0;
			VecScale(f[0], r12, -1);
			f[1] = r12;
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "pi bond force %d-%d: r=%f, r0=%f, k1=%f, {%f %f %f}\n", pibonds[0]+1, pibonds[1]+1, w1, pars->bond.r0, k1, r12.x, r12.y, r12.z);
			}
		} else if (j == 3) {  /*  Angle  */
			if (s_calc_angle_force_one(&r[0], &r[1], &r[2], &pars->angle, &val, &energy, &valk1, &f[0], &f[1], &f[2]) != 0)
				continue;
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "pi angle force %d-%d-%d: a=%f, a0=%f, k1=%f, {%f %f %f}, {%f %f %f}\n", pibonds[0]+1, pibonds[1]+1, pibonds[2]+1, val*180/PI, pars->angle.a0*180/PI, energy, f[0].x, f[0].y, f[0].z, f[1].x, f[1].y, f[1].z);
			}
		} else if (j == 4) {  /*  Dihedral  */
			if (s_calc_dihedral_force_one(&r[0], &r[1], &r[2], &r[3], &pars->torsion, &val, &energy, &valk1, &f[0], &f[1], &f[2], &f[3]) != 0)
				continue;
			if (arena->debug_result && arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "pi dihedral force %d-%d-%d-%d: phi=%f, k1=%f, {%f %f %f}, {%f %f %f}, {%f %f %f}\n", pibonds[0]+1, pibonds[1]+1, pibonds[2]+1, pibonds[3]+1, val*180/PI, energy, f[0].x, f[0].y, f[0].z, f[1].x, f[1].y, f[1].z, f[2].x, f[2].y, f[2].z);
			}
		} else continue;
	}
	
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
	int dx, dy, nn;
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

/*  Count all lattice points within distance d from point v  */
/*  Ref. "Algorithms for the Shortest and Closest Lattice Vector Problems"
    Guillaume Hanrot, Xavier Pujol, Damien Stehle
    Coding and Cryptology; Lecture Notes in Computer Science, 2011 vol. 6639/2011 pp. 159-190
    DOI: 10.1007/978-3-642-20901-7_10  */
static int
s_enum_neighbors(MDArena *arena, Vector v, Double d)
{
	Vector bn[3];  /*  Base vectors for the cell axes  */
	Double bl[3];  /*  Square of the length of base vectors  */
	Double mu[3];  /*  Non-diagonal terms for Gram-Schmidt orthogonalization  */
	XtalCell *cell = arena->mol->cell;
	Int dim;       /*  Number of periodic axes  */
	Int ijk[3];    /*  Renumber the indices. For 1-dimensional case, axes[ijk[0]] is the periodic axis.
					   For 2-dimensional case, axes[ijk[0]] and axes[ijk[1]] are the periodic axes.
					   Only uses exchange of two elements, so that ijk[ijk[i]] == i for i=0,1,2  */
	Vector axes[3]; /*  Renumbered cell axes  */
	Double t[3];   /*  Projected vector v on bn[]  */
	Int count;     /*  Number of results  */
	Int x[3];      /*  Lattice point candidate  */
	Double r[3];
	Int i;
	Double w;
	
	dim = (arena->periodic_a != 0) + (arena->periodic_b != 0) + (arena->periodic_c != 0);
	ijk[0] = 0; ijk[1] = 1; ijk[2] = 2;
	if (dim == 0) {
		/*  Non-periodic case: check whether (0,0,0) is within d or not  */
		if (VecLength(v) < d) {
			x[0] = x[1] = x[2] = 0;
			AssignArray(&arena->lattice_offsets, &arena->nlattice_offsets, sizeof(Int) * 3, 0, x);
			return 1;
		} else return 0;
	}
	
	if (dim == 1) {
		if (arena->periodic_b) {
			ijk[0] = 1; ijk[1] = 0;
		} else if (arena->periodic_c) {
			ijk[0] = 2; ijk[2] = 0;
		}
	} else if (dim == 2) {
		if (!arena->periodic_a) {
			ijk[2] = 0; ijk[0] = 2;
		} else if (!arena->periodic_b) {
			ijk[1] = 0; ijk[0] = 1;
		}
	}
	axes[0] = cell->axes[ijk[0]];
	axes[1] = cell->axes[ijk[1]];
	axes[2] = cell->axes[ijk[2]];

	/*  Gram-Schmidt orthogonalization  */
	bn[0] = axes[0];
	bl[0] = VecLength2(bn[0]);
	mu[0] = VecDot(axes[1], bn[0]) / bl[0];  /*  bl[0] should be non-zero  */
	if (dim == 1) {
		VecZero(bn[1]);
		VecZero(bn[2]);
		bl[1] = bl[2] = 0.0;
		mu[1] = mu[2] = 0.0;
	} else {
		bn[1] = axes[1];
		VecScaleInc(bn[1], bn[0], -mu[0]);
		bl[1] = VecLength2(bn[1]);
		mu[1] = VecDot(axes[2], bn[0]) / bl[0];
		if (dim == 2 || bl[1] < 1e-10) {
			VecZero(bn[2]);
			bl[2] = mu[2] = 0.0;
			if (bl[1] < 1e-10)
				dim = 1;
			else dim = 2;
		} else {
			mu[2] = VecDot(axes[3], bn[1]) / bl[1];
			bn[2] = axes[2];
			VecScaleInc(bn[2], bn[0], -mu[1]);
			VecScaleInc(bn[2], bn[1], -mu[2]);
			bl[2] = VecLength2(bn[2]);
			if (bl[2] < 1e-10)
				dim = 2;
			else dim = 3;
		}
	}
	
	/*  Project the target vector  */
	t[0] = t[1] = t[2] = 0.0;
	t[0] = VecDot(v, bn[0]) / bl[0];
	if (dim >= 2) {
		t[1] = VecDot(v, bn[1]) / bl[1];
		if (dim >= 3) {
			t[2] = VecDot(v, bn[2]) / bl[2];
		}
	}
	
	/*  Enumerate  */
	count = 0;
	x[0] = x[1] = x[2] = 0;
	r[0] = r[1] = r[2] = 0.0;
	i = dim - 1;
	w = -10000.0;
	x[i] = ceil(t[i] - d / sqrt(bl[i]));
	while (1) {
		Int j;
		Double w2;
		if (w == -10000.0) {
			w = 0.0;
			for (j = i + 1; j < dim; j++) {
				w += x[j] * mu[j + i - 1];
			}
		}
		w2 = x[i] - t[i] + w;
		r[i] = w2 * w2 * bl[i];
		w2 = 0.0;
		for (j = i; j < dim; j++) {
			w2 += r[j];
		}
		w2 = d * d - w2;
		if (w2 >= 0.0) {
			if (i == 0) {
				/*  Found  */
				Int xx[3];
				xx[0] = x[ijk[0]];
				xx[1] = x[ijk[1]];
				xx[2] = x[ijk[2]];
				AssignArray(&arena->lattice_offsets, &arena->nlattice_offsets, sizeof(Int) * 3, count, xx);
				count++;
				x[0]++;
				w = -10000.0;
			} else {
				/*  Step down  */
				i--;
				w = 0.0;
				for (j = i + 1; j < dim; j++) {
					w += x[j] * mu[j + i - 1];
				}
				x[i] = ceil(t[i] - w - sqrt(w2 / bl[i]));
			}
		} else {
			i++;
			if (i >= dim)
				break;
			x[i]++;
			w = -10000.0;
		}
	}
	return count;
}

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
	Parameter *par = arena->par;
	Double vdw_cutoff;
	
	natoms = mol->natoms;
	for (i = 0; i < natoms; i++)
		VecZero(vdr[i]);

	ndx = (arena->periodic_a ? 1 : 0);
	ndy = (arena->periodic_b ? 1 : 0);
	ndz = (arena->periodic_c ? 1 : 0);
	limit = arena->pairlist_distance * arena->pairlist_distance;
	n = 0;
	use_sym = (arena->mol->nsyms > 0);
	for (i = 0, api = atoms; i < natoms; i++, api++) {
		MDExclusion *exinfo = arena->exinfo + i;
		Int index4 = (exinfo + 1)->index0;
		Int vdw_idx1 = arena->vdw_par_i[i];
		Int vdw_idx, vdw_idx2;
		arena->verlet_i[i] = n;

		for (j = i, apj = atoms + j; j < natoms; j++, apj++) {
			Vector rij;
			Double lenij2;
			Int *ip, index0;
			int exflag = 0;
			int mult = 1;
		/*	int dxbase, dybase, dzbase; */
			int count;

			/*  Fixed atoms  */
			if (api->fix_force < 0 && apj->fix_force < 0)
				continue;

			/*  Non-occupied atoms  */
			if (api->mm_exclude || apj->mm_exclude)
				continue;

			/*  If no symmetry expansion, then no interaction with self  */
			if (ndx + ndy + ndz == 0 && i == j)
				continue;

			VecSub(rij, apj->r, api->r);

			/*  Calculate the cell offset for the nearest neighbor  */
			if (ndx + ndy + ndz != 0 && apj->periodic_exclude == 0 && cell != NULL) {
				count = s_enum_neighbors(arena, rij, limit);
			} else {
				static Int sZeros[3] = {0, 0, 0};
				AssignArray(&arena->lattice_offsets, &arena->nlattice_offsets, sizeof(Int) * 3, 0, sZeros);
				count = 1;
			}
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
				if (((cr->type1 == kAtomTypeWildcard || cr->type1 == api->type) && (cr->type2 == kAtomTypeWildcard || cr->type2 == apj->type)) ||
					((cr->type1 == kAtomTypeWildcard || cr->type1 == apj->type) && (cr->type2 == kAtomTypeWildcard || cr->type2 == api->type))) {
					vdw_cutoff = cr->cutoff;
					break;
				}
				if ((cr->type1 == i && cr->type2 == j) || (cr->type1 == j && cr->type2 == i)) {
					vdw_cutoff = cr->cutoff;
					break;
				}
			}
			if (vdw_cutoff < 0) {
				/*  A negative value of specific cutoff means "cutoff at r_eq * (-specific_cutoff)  */
				MDVdwCache *cp = &(arena->vdw_cache[vdw_idx]);
				Double r_eq = pow(cp->par.A / cp->par.B * 2.0, 1.0/6.0);
				vdw_cutoff = r_eq * (-vdw_cutoff);
			}

			/*  Search for pairs  */
			for (nn = 0; nn < count; nn++) {
				Vector rij0 = rij;
				dx = arena->lattice_offsets[nn * 3];
				dy = arena->lattice_offsets[nn * 3 + 1];
				dz = arena->lattice_offsets[nn * 3 + 2];
				if (dx == 0 && dy == 0 && dz == 0) {
					/*  Pair within the unit cell  */
					if (i == j)
						continue;  /*  No interaction with self  */
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
					VecScaleInc(rij0, cell->axes[0], dx);
					VecScaleInc(rij0, cell->axes[1], dy);
					VecScaleInc(rij0, cell->axes[2], dz);
				/*	VecInc(rij0, cell_offsets[nn]); */
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
					vlp->symop.dx = dx;
					vlp->symop.dy = dy;
					vlp->symop.dz = dz;
					vlp->symop.sym = 0;
					vlp->symop.alive = (vlp->symop.dx != 0 || vlp->symop.dy != 0 || vlp->symop.dz != 0);
					vlp->index = vdw_idx;
					vlp->n1 = i;
					vlp->n2 = j;
					vlp->vdw_cutoff = vdw_cutoff;
					vlp->length = sqrt(lenij2);
					n++;
				} /* end if lenij2 <= limit */
			} /* end loop nn */
		} /* end loop j */

	} /* end loop i */
	arena->verlet_i[natoms] = n;
	arena->nverlets = n;
	arena->verlets = vl;
	arena->last_verlet_step = arena->step;
		
	if (arena->debug_result && arena->debug_output_level > 2) {
		fprintf(arena->debug_result, "\n  Verlet list at step %d\n", arena->step);
		fprintf(arena->debug_result, "  {atom1 atom2 vdw_type (=1 if 1-4 bonded) vdw_index cell_translation}\n");
		for (i = 0; i < arena->nverlets; i++) {
			fprintf(arena->debug_result, "{%d %d %d %d %06d}%c", vl[i].n1+1, vl[i].n2+1, vl[i].vdw_type, vl[i].index, (vl[i].symop.dx + 50) * 10000 + (vl[i].symop.dy + 50) * 100 + (vl[i].symop.dz + 50), (i % 4 == 3 ? '\n' : ' '));
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
	Double lambda, dlambda;

	MDVerlet *vl;
/*	MDVdwCache *vdw_cache = arena->vdw_cache; */
	Atom *atoms = arena->mol->atoms;
	XtalCell *cell = arena->mol->cell;

#if MINIMIZE_CELL
	memset(arena->cell_forces, 0, sizeof(Double) * 12);
#endif
	
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
		Vector rij, fij, rij2, fij2;
		Double r2, w2, w6, w12;
		if (group_flags_1 != NULL && group_flags_2 != NULL) {
			if (!(get_group_flag(group_flags_1, vl->n1) && get_group_flag(group_flags_2, vl->n2))
			&& !(get_group_flag(group_flags_1, vl->n2) && get_group_flag(group_flags_2, vl->n1)))
				continue;
		}

		if (arena->nalchem_flags > 0) {
			char c1, c2;
			if (vl->n1 < arena->nalchem_flags)
				c1 = arena->alchem_flags[vl->n1];
			else c1 = 0;
			if (vl->n2 < arena->nalchem_flags)
				c2 = arena->alchem_flags[vl->n2];
			else c2 = 0;
			if ((c1 == 1 && c2 == 2) || (c1 == 2 && c2 == 1))
				continue;
			if (c1 == 1 || c2 == 1) {
				lambda = (1.0 - arena->alchem_lambda);
				dlambda = -arena->alchem_dlambda;
			} else if (c1 == 2 || c2 == 2) {
				lambda = arena->alchem_lambda;
				dlambda = arena->alchem_dlambda;
			} else {
				lambda = 1.0;
				dlambda = 0.0;
			}
		} else {
			lambda = 1.0;
			dlambda = 0.0;
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
		fij2 = fij;

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
			fij2 = fij;
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
		k0 /= vl->mult;
		k1 *= lambda;
		VecScale(fij, rij, k1);
		*energies += k0 * lambda;
		if (forces != NULL) {
			VecDec(forces[vl->n1], fij);
			VecInc(forces[vl->n2], fij);
		}
		VecInc(fij2, fij);
		if (dlambda != 0.0)
			arena->alchem_energy += k0 * dlambda;

		if (arena->debug_result && arena->debug_output_level > 1) {
			fprintf(arena->debug_result, "nonbonded(vdw) force %d-%d: r=%f, k0=%f, k1=%f, {%f %f %f}\n", vl->n1+1, vl->n2+1, sqrt(r2), k0/KCAL2INTERNAL, k1*sqrt(r2)/KCAL2INTERNAL, fij.x/KCAL2INTERNAL, fij.y/KCAL2INTERNAL, fij.z/KCAL2INTERNAL);
		}
		
#if MINIMIZE_CELL
		if (arena->is_minimizing && arena->minimize_cell && (ap2->symop.alive || vl->symop.alive)) {
			/*  Force for changing cell parameters  */
			Vector rr;
			Int j, k;
			Transform da, tr, symtr, symtr2;
			if (ap2->symop.alive)
				rij2 = ATOM_AT_INDEX(arena->mol->atoms, ap2->symbase)->r;
			else rij2 = ap2->r;
			TransformVec(&rr, cell->rtr, &rij2);
			memmove(symtr, SYMMETRY_AT_INDEX(arena->mol->syms, ap2->symop.sym), sizeof(Transform));
			symtr[9] += ap2->symop.dx + vl->symop.dx;
			symtr[10] += ap2->symop.dy + vl->symop.dy;
			symtr[11] += ap2->symop.dz + vl->symop.dz;
			TransformMul(symtr2, symtr, cell->rtr);
			TransformMul(symtr2, cell->tr, symtr2);
			for (j = 0; j < 12; j++) {
				Vector rrr;
				if (arena->periodic_a == 0 && (j == 0 || j == 1 || j == 2 || j == 3 || j == 6 || j == 9))
					continue;
				if (arena->periodic_b == 0 && (j == 1 || j == 3 || j == 4 || j == 5 || j == 7 || j == 10))
					continue;
				if (arena->periodic_c == 0 && (j == 2 || j == 5 || j == 6 || j == 7 || j == 8 || j == 11))
					continue;
				memset(da, 0, sizeof(Transform));
				da[j] = 1.0;
				TransformMul(tr, symtr2, da);
				/*  The translation part should be adjusted, because the derivative matrix
					has the (3,3) element as 0.0 (instead of 1.0). Thus,
				    | R v | | R' v' | = | R*R' R*v' |   (instead of R*v'+v)
				    | 0 1 | | 0  0  |   | 0    0    |
				    da * symtr does not have this problem, because the derivative matrix is
				    multipled from the left.  */
				tr[9] -= symtr2[9];
				tr[10] -= symtr2[10];
				tr[11] -= symtr2[11];
				TransformMul(da, da, symtr);
				for (k = 0; k < 12; k++)
					da[k] -= tr[k];
				TransformVec(&rrr, da, &rr);
				arena->cell_forces[j] += VecDot(fij2, rrr);
			}
		}
#endif
		
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
	
	/*  External forces (used for artificial deformation of the structure)  */
	if (arena->nexforces > 0) {
		for (i = 0; i < arena->nexforces; i++) {
			VecInc(forces[i], arena->exforces[i]);
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
	arena->alchem_energy = 0.0;

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
		if (mol->atoms[i].mm_exclude)
			continue;
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
