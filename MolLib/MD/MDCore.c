/*
 *  MDCore.c
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

#include "MDCore.h"
#include "MDGraphite.h"
#include "MDPressure.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#if __WXMSW__
#define ftello(x) ftell(x)
#endif

static char sErrBuf[256];

void
md_panic(MDArena *arena, const char *fmt,...)
{
	va_list ap;
	jmp_buf *envp;

	/*  Clean up the running simulation  */
	md_finish(arena);
	arena->is_running = 0;
	arena->mol->needsMDRebuild = 1;

	va_start(ap, fmt);
	vsnprintf(arena->errmsg, sizeof(arena->errmsg), fmt, ap);
	va_end(ap);

	envp = arena->setjmp_buf;
	arena->setjmp_buf = NULL;
	if (arena->md_panic_func != NULL)
		(*(arena->md_panic_func))(arena, arena->errmsg);
	if (envp != NULL)
		longjmp(*envp, 1);
	else {
		fprintf(stderr, "%s\n", sErrBuf);
		exit(1);
	}
}

/*  Message output  */
extern int MyAppCallback_showScriptMessage(const char *fmt, ...);
	
int
md_warning(MDArena *arena, const char *fmt, ...)
{
	va_list ap;
	char buf[1024];
	va_start(ap, fmt);
	vsnprintf(buf, sizeof buf, fmt, ap);
	va_end(ap);
	MyAppCallback_showScriptMessage("%s", buf);
	if (arena->log_result != NULL)
		fputs(buf, arena->log_result);
	return strlen(buf);
}

int
md_log(MDArena *arena, const char *fmt, ...)
{
	va_list ap;
	char buf[1024];
	if (arena->log_result == NULL)
		return 0;
	if (fmt == NULL) {
		fflush(arena->log_result);
		return 0;
	}
	va_start(ap, fmt);
	vsnprintf(buf, sizeof buf, fmt, ap);
	va_end(ap);
	fputs(buf, arena->log_result);
	return strlen(buf);
}

void
md_debug(MDArena *arena, const char *fmt,...)
{
	va_list ap;
	if (arena->debug_result == NULL || arena->debug_output_level == 0)
		return;
	va_start(ap, fmt);
	vfprintf(arena->debug_result, fmt, ap);
	va_end(ap);
	fflush(arena->debug_result);
}

#if DEBUG
#define DEBUG_FIT_COORDINATES 0  /*  set to 1 while debugging md_fit_coordinates */
#endif

/*  Calculate the transform that moves the current coordinates to the reference
  coordinates with least displacements. The reference coordinates are given
  as one of the snapshots (refno = zero-based)  */
int
md_fit_coordinates(MDArena *arena, Int refno, Double *weights, Transform trans)
{
	Atom *ap;
	Vector *rvf;
	Int natoms, nn;
	Vector org1, org2;
	Int i, j, k;
	Double w, w1;
	Mat33 r, q, u;
	Double eigen_val[3];
	Vector eigen_vec[3];
	Vector s[3];
	if (arena == NULL || arena->mol == NULL || arena->mol->natoms == 0)
		return 1;  /*  Molecule is empty  */
	if (refno < 0 || refno >= arena->nsnapshots)
		return 2;  /*  No such snapshot  */
	natoms = arena->mol->natoms;
	ap = arena->mol->atoms;
	rvf = arena->snapshots[refno]->rvf;
	
	/*  Calculate the weighted center  */
	for (j = 0; j < 2; j++) {
		Vector *op = (j == 0 ? &org1 : &org2);
		VecZero(*op);
		w = 0.0;
		for (i = 0; i < natoms; i++) {
			w1 = (weights != NULL ? weights[i] : ap[i].weight);
			if (w1 == 0.0)
				continue;
			if (j == 0)
				VecScaleInc(*op, ap[i].r, w1);
			else
				VecScaleInc(*op, rvf[i * 3], w1);
			w += w1;
		}
		w = 1.0 / w;
		VecScaleSelf(*op, w);
	}
	
    /*  R = sum(weight[n] * x[n] * t(y[n]));  */
    /*  Matrix to diagonalize = R * tR    */
	memset(r, 0, sizeof(Mat33));
	memset(q, 0, sizeof(Mat33));
	memset(u, 0, sizeof(Mat33));
	nn = 0;
	for (i = 0; i < natoms; i++) {
		Vector v1, v2;
		w1 = (weights != NULL ? weights[i] : ap[i].weight);
		if (w1 == 0.0)
			continue;
		VecSub(v1, ap[i].r, org1);
		VecSub(v2, rvf[3 * i], org2);
		r[0] += w1 * v2.x * v1.x;
		r[1] += w1 * v2.x * v1.y;
		r[2] += w1 * v2.x * v1.z;
		r[3] += w1 * v2.y * v1.x;
		r[4] += w1 * v2.y * v1.y;
		r[5] += w1 * v2.y * v1.z;
		r[6] += w1 * v2.z * v1.x;
		r[7] += w1 * v2.z * v1.y;
		r[8] += w1 * v2.z * v1.z;
/*		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				r[3*j+k] += w1 * VecIndex(&v2, j) * VecIndex(&v1, k);
#if DEBUG_FIT_COORDINATES
				printf("r(%d,%d) += %.6g * %.6g * %.6g (%.6g)\n", j, k, w1, VecIndex(&v2, j), VecIndex(&v1, k), r[3*j+k]);
#endif
			}
		}
*/
		nn++;
	}
	for (i = 0; i < 9; i++)
		r[i] /= (nn * nn);
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				q[i*3+j] += r[i*3+k] * r[j*3+k];
			}
		}
	}

#if DEBUG_FIT_COORDINATES
	printf("Matrix to diagonalize:\n");
	for (i = 0; i < 3; i++) {
		printf("%10.6g %10.6g %10.6g\n", q[i*3], q[i*3+1], q[i*3+2]);
	}
#endif

	if (MatrixSymDiagonalize(q, eigen_val, eigen_vec) != 0)
		return 3;  /*  Cannot determine the eigenvector  */

#if DEBUG_FIT_COORDINATES
	for (i = 0; i < 3; i++) {
		printf("Eigenvalue %d = %.6g\n", i+1, eigen_val[i]);
		printf("Eigenvector %d: %.6g %.6g %.6g\n", i+1, eigen_vec[i].x, eigen_vec[i].y, eigen_vec[i].z);
	}
#endif

    /*  s[i] = normalize(tR * v[i])  */
    /*  U = s0*t(v0) + s1*t(v1) + s2*t(v2)  */
	MatrixTranspose(r, r);
	for (i = 0; i < 3; i++) {
		MatrixVec(&s[i], r, &eigen_vec[i]);
		w1 = 1.0 / VecLength(s[i]);
		VecScaleSelf(s[i], w1);
#if DEBUG_FIT_COORDINATES
		printf("s[%d] = %.6g %.6g %.6g\n", i+1, s[i].x, s[i].y, s[i].z);
#endif
	}
/*	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				u[3*i+j] += VecIndex(&s[k], j) * VecIndex(&eigen_vec[k], i);
			}
		}
	}
*/
	for (k = 0; k < 3; k++) {
		u[0] += s[k].x * eigen_vec[k].x;
		u[1] += s[k].y * eigen_vec[k].x;
		u[2] += s[k].z * eigen_vec[k].x;
		u[3] += s[k].x * eigen_vec[k].y;
		u[4] += s[k].y * eigen_vec[k].y;
		u[5] += s[k].z * eigen_vec[k].y;
		u[6] += s[k].x * eigen_vec[k].z;
		u[7] += s[k].y * eigen_vec[k].z;
		u[8] += s[k].z * eigen_vec[k].z;
	}
	
	/*  y = U*(x - org1) + org2 = U*x + (org2 - U*org1)  */
	MatrixVec(&org1, u, &org1);
	VecDec(org2, org1);
	for (i = 0; i < 9; i++)
		trans[i] = u[i];
	trans[9] = org2.x;
	trans[10] = org2.y;
	trans[11] = org2.z;
	
	return 0;
}

static void
s_register_missing_parameters(Int **missing, Int *nmissing, Int type, Int t1, Int t2, Int t3, Int t4)
{
	Int *mp;
	Int i;
	for (i = 0, mp = *missing; i < *nmissing; i++, mp += 5) {
		if (mp[0] == type && mp[1] == t1 && mp[2] == t2 && mp[3] == t3 && mp[4] == t4)
			return; /* Already registered */
	}
	mp = (Int *)AssignArray(missing, nmissing, sizeof(Int)*5, *nmissing, NULL);
	if (mp != NULL) {
		mp[0] = type;
		mp[1] = t1;
		mp[2] = t2;
		mp[3] = t3;
		mp[4] = t4;
	}
}

/*  Check the bonded atoms and append to results if not already present */
/*  results[] is terminated by -1, hence must be at least (natom+1) size  */
static int
s_check_bonded(Atom *ap, Int *results)
{
	int i, n, *ip;
	for (i = 0; i < ap->nconnects; i++) {
		n = ap->connects[i];
		for (ip = results; *ip >= 0; ip++) {
			if (n == *ip)
				break;
		}
		if (*ip < 0) {
			*ip++ = n;
			*ip = -1;
		}
	}
	for (ip = results; *ip >= 0; ip++)
		;
	return ip - results;
}

static void
s_make_exclusion_list(MDArena *arena)
{
	Int *results;
	Int natoms = arena->mol->natoms;
	Atom *atoms = arena->mol->atoms;
	MDExclusion *exinfo;
	int next_index, i, j;

	results = (Int *)calloc(sizeof(Int), natoms + 1);
	if (results == NULL)
		md_panic(arena, ERROR_out_of_memory);

	if (arena->exlist != NULL) {
		free(arena->exlist);
		arena->exlist = NULL;
	}
	arena->nexlist = 0;

	if (arena->exinfo != NULL)
		free(arena->exinfo);
	arena->exinfo = (MDExclusion *)calloc(sizeof(MDExclusion), natoms + 1);
	if (arena->exinfo == NULL)
		md_panic(arena, ERROR_out_of_memory);
	exinfo = arena->exinfo;

	next_index = 0;
	for (i = 0; i < natoms; i++) {
		int n;
		exinfo[i].index0 = 0;
		/* special exclusion: only self */
		results[0] = i;
		results[1] = -1;
		exinfo[i].index1 = 1;
		/*  1-2 exclusion (directly bonded)  */
		exinfo[i].index2 = s_check_bonded(&atoms[i], results);
		n = exinfo[i].index2;
		/*  1-3 exclusion: atoms bonded to 1-2 exclusions  */
		for (j = exinfo[i].index1; j < exinfo[i].index2; j++)
			n = s_check_bonded(&atoms[results[j]], results);
		exinfo[i].index3 = n;
		/*  1-4 exclusion: atoms bonded to 1-3 exclusions  */
		for (j = exinfo[i].index2; j < exinfo[i].index3; j++)
			n = s_check_bonded(&atoms[results[j]], results);
		AssignArray(&arena->exlist, &arena->nexlist, sizeof(Int), next_index + n, NULL);
		memcpy(arena->exlist + next_index, results, n * sizeof(Int));
		exinfo[i].index0 += next_index;
		exinfo[i].index1 += next_index;
		exinfo[i].index2 += next_index;
		exinfo[i].index3 += next_index;
		next_index += n;
	}
	exinfo[natoms].index0 = next_index;  /*  End of exlist  */
	
	free(results);
}

static int
s_lookup_improper_pars(Parameter *par, Int type1, Int type2, Int type3, Int type4)
{
	Int idx;
	Int t1, t2, t3, t4;
	for (idx = par->nimproperPars - 1; idx >= 0; idx--) {
		t1 = par->improperPars[idx].type1;
		t2 = par->improperPars[idx].type2;
		t3 = par->improperPars[idx].type3;
		t4 = par->improperPars[idx].type4;
		if (t1 == -3)
			continue;  /*  Custom parameter  */
		if (t3 >= 0 && t3 != type3)
			continue;
		if ((t1 < 0 || t1 == type1) &&
			(t2 < 0 || t2 == type2) &&
			(t4 < 0 || t4 == type4))
			break;
		if ((t1 < 0 || t1 == type1) &&
			(t2 < 0 || t2 == type4) &&
			(t4 < 0 || t4 == type2))
			break;
		if ((t1 < 0 || t1 == type2) &&
			(t2 < 0 || t2 == type1) &&
			(t4 < 0 || t4 == type4))
			break;
		if ((t1 < 0 || t1 == type2) &&
			(t2 < 0 || t2 == type4) &&
			(t4 < 0 || t4 == type1))
			break;
		if ((t1 < 0 || t1 == type4) &&
			(t2 < 0 || t2 == type1) &&
			(t4 < 0 || t4 == type2))
			break;
		if ((t1 < 0 || t1 == type4) &&
			(t2 < 0 || t2 == type2) &&
			(t4 < 0 || t4 == type1))
			break;
	}
	return idx;
}

int
md_check_abnormal_bond(MDArena *arena, Molecule *mol, int idx)
{
	BondPar *bp;
	Atom *ap1, *ap2;
	Int idx2;
	Double d;
	if (mol == NULL)
		mol = arena->mol;
	if (arena->par == NULL || idx < 0 || idx >= mol->nbonds || (idx2 = arena->bond_par_i[idx]) < 0 || idx2 >= arena->par->nbondPars)
		return -1;
	ap1 = ATOM_AT_INDEX(mol->atoms, mol->bonds[idx * 2]);
	ap2 = ATOM_AT_INDEX(mol->atoms, mol->bonds[idx * 2 + 1]);
	bp = arena->par->bondPars + idx2;
	if (bp->k == 0.0 || bp->r0 == 0.0)
		return 0;
	d = MoleculeMeasureBond(mol, &(ap1->r), &(ap2->r));
	return (fabs(d / bp->r0 - 1.0) >= 0.2);
}

int
md_check_abnormal_angle(MDArena *arena, Molecule *mol, int idx)
{
	AnglePar *anp;
	Atom *ap1, *ap2, *ap3;
	Int idx2;
	Double d;
	if (mol == NULL)
		mol = arena->mol;
	if (arena->par == NULL || idx < 0 || idx >= mol->nangles || (idx2 = arena->angle_par_i[idx]) < 0 || idx2 >= arena->par->nanglePars)
		return -1;
	ap1 = ATOM_AT_INDEX(mol->atoms, mol->angles[idx * 3]);
	ap2 = ATOM_AT_INDEX(mol->atoms, mol->angles[idx * 3 + 1]);
	ap3 = ATOM_AT_INDEX(mol->atoms, mol->angles[idx * 3 + 2]);
	anp = arena->par->anglePars + idx2;
	if (anp->k == 0.0 || anp->a0 == 0.0)
		return 0;
	d = MoleculeMeasureAngle(mol, &(ap1->r), &(ap2->r), &(ap3->r));
	return (fabs(d - anp->a0 * kRad2Deg) >= 20.0);
}

int
md_check_abnormal_dihedral(MDArena *arena, Molecule *mol, int idx)
{
	TorsionPar *tp;
	Atom *ap1, *ap2, *ap3, *ap4;
	Int idx2;
	Double d;
	if (mol == NULL)
		mol = arena->mol;
	if (arena->par == NULL || idx < 0 || idx >= mol->ndihedrals || (idx2 = arena->dihedral_par_i[idx]) < 0 || idx2 >= arena->par->ndihedralPars)
		return -1;
	ap1 = ATOM_AT_INDEX(mol->atoms, mol->dihedrals[idx * 4]);
	ap2 = ATOM_AT_INDEX(mol->atoms, mol->dihedrals[idx * 4 + 1]);
	ap3 = ATOM_AT_INDEX(mol->atoms, mol->dihedrals[idx * 4 + 2]);
	ap4 = ATOM_AT_INDEX(mol->atoms, mol->dihedrals[idx * 4 + 3]);
	tp = arena->par->dihedralPars + idx2;
	if (tp->k[0] == 0.0)
		return 0;
	d = MoleculeMeasureDihedral(mol, &(ap1->r), &(ap2->r), &(ap3->r), &(ap4->r));
	if (tp->period[0] == 0)
		return (fabs(d - tp->phi0[0] * kRad2Deg) >= 20.0);
	else {
		d = (tp->period[0] * (d - tp->phi0[0] * kRad2Deg)) / 360.0 + 0.5;
		d = (d - floor(d) - 0.5) * 360.0; /* map to [-180, 180] */
		return (fabs(d) >= 20.0);
	}
}

int
md_check_abnormal_improper(MDArena *arena, Molecule *mol, int idx)
{
	TorsionPar *tp;
	Atom *ap1, *ap2, *ap3, *ap4;
	Int idx2;
	Double d;
	if (mol == NULL)
		mol = arena->mol;
	if (arena->par == NULL || idx < 0 || idx >= mol->ndihedrals || (idx2 = arena->improper_par_i[idx]) < 0 || idx2 >= arena->par->nimproperPars)
		return -1;
	ap1 = ATOM_AT_INDEX(mol->atoms, mol->impropers[idx * 4]);
	ap2 = ATOM_AT_INDEX(mol->atoms, mol->impropers[idx * 4 + 1]);
	ap3 = ATOM_AT_INDEX(mol->atoms, mol->impropers[idx * 4 + 2]);
	ap4 = ATOM_AT_INDEX(mol->atoms, mol->impropers[idx * 4 + 3]);
	tp = arena->par->improperPars + idx2;
	if (tp->k[0] == 0.0)
		return 0;
	d = MoleculeMeasureDihedral(mol, &(ap1->r), &(ap2->r), &(ap3->r), &(ap4->r));
	if (tp->period[0] == 0)
		return (fabs(d - tp->phi0[0] * kRad2Deg) >= 20.0);
	else {
		d = (tp->period[0] * (d - tp->phi0[0] * kRad2Deg)) / 360.0 + 0.5;
		d = (d - floor(d) - 0.5) * 360.0; /* map to [-180, 180] */
		return (fabs(d) >= 20.0);
	}
}

static int
s_search_bond(MDArena *arena, int n1, int n2)
{
	int i, *ip;
	if (n1 < 0 || n1 >= arena->mol->natoms || n2 < 0 || n2 >= arena->mol->natoms)
		return -1;
	for (i = 0, ip = arena->mol->bonds; i < arena->mol->nbonds; i++, ip += 2) {
		if ((ip[0] == n1 && ip[1] == n2) || (ip[0] == n2 && ip[1] == n1))
			return i;
	}
	return -1;
}

static int
s_search_angle(MDArena *arena, int n1, int n2, int n3)
{
	int i, *ip;
	if (n1 < 0 || n1 >= arena->mol->natoms || n2 < 0 || n2 >= arena->mol->natoms || n3 < 0 || n3 >= arena->mol->natoms)
		return -1;
	for (i = 0, ip = arena->mol->angles; i < arena->mol->nangles; i++, ip += 3) {
		if (ip[1] == n2 && ((ip[0] == n1 && ip[2] == n3) || (ip[0] == n3 && ip[2] == n1)))
			return i;
	}
	return -1;
}

static int
s_search_dihedral(MDArena *arena, int n1, int n2, int n3, int n4)
{
	int i, *ip;
	if (n1 < 0 || n1 >= arena->mol->natoms || n2 < 0 || n2 >= arena->mol->natoms || n3 < 0 || n3 >= arena->mol->natoms || n4 < 0 || n4 >= arena->mol->natoms)
		return -1;
	for (i = 0, ip = arena->mol->dihedrals; i < arena->mol->ndihedrals; i++, ip += 4) {
		if ((ip[0] == n1 && ip[1] == n2 && ip[2] == n3 && ip[3] == n4) || (ip[0] == n4 && ip[1] == n3 && ip[2] == n2 && ip[3] == n1))
			return i;
	}
	return -1;
}

static int
s_search_improper(MDArena *arena, int n1, int n2, int n3, int n4)
{
	int i, *ip;
	if (n1 < 0 || n1 >= arena->mol->natoms || n2 < 0 || n2 >= arena->mol->natoms || n3 < 0 || n3 >= arena->mol->natoms || n4 < 0 || n4 >= arena->mol->natoms)
		return -1;
	for (i = 0, ip = arena->mol->impropers; i < arena->mol->nimpropers; i++, ip += 4) {
		if ((ip[0] == n1 && ip[1] == n2 && ip[2] == n3 && ip[3] == n4) || (ip[0] == n4 && ip[1] == n3 && ip[2] == n2 && ip[3] == n1))
			return i;
	}
	return -1;
}

/*  Find vdw parameters and build the in-use list */
static int
s_find_vdw_parameters(MDArena *arena)
{
	Int idx, i, j, type1, type2, t1, t2, nmissing;
	Double cutoff6, cutoff12;
	Parameter *par = arena->par;
	Molecule *mol = arena->mol;
	Atom *ap;
	VdwPar *vp;

	if (arena->vdw_par_i != NULL)
		free(arena->vdw_par_i);
	arena->vdw_par_i = (Int *)calloc(sizeof(Int), mol->natoms);
	if (arena->vdw_par_i == NULL)
		md_panic(arena, ERROR_out_of_memory);

	/*  Find the vdw parameter; priority: (1) variant-aware in mol->par, (2) mol->par ignoring variants,
	    (3) variant->aware in global par, (4) global par ignoring variants  */
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (mol->par != NULL) {
			vp = ParameterLookupVdwPar(mol->par, ap->type, 0);
		/*	vp = ParameterLookupVdwPar(mol->par, i, kParameterLookupLocal | kParameterLookupNoWildcard);
			if (vp == NULL)
				vp = ParameterLookupVdwPar(mol->par, ap->type, kParameterLookupLocal | kParameterLookupGlobal); */
			if (vp != NULL) {
				arena->vdw_par_i[i] = (vp - mol->par->vdwPars) + ATOMS_MAX_NUMBER * 2;
				continue;
			}
		}
		vp = ParameterLookupVdwPar(gBuiltinParameters, ap->type, 0);
		if (vp != NULL) {
			arena->vdw_par_i[i] = (vp - gBuiltinParameters->vdwPars) + ATOMS_MAX_NUMBER;
			continue;
		}
		/*  Record as missing  */
		vp = ParameterLookupVdwPar(par, ap->type, kParameterLookupMissing | kParameterLookupNoBaseAtomType);
		if (vp == NULL) {
			vp = AssignArray(&par->vdwPars, &par->nvdwPars, sizeof(VdwPar), par->nvdwPars, NULL);
			vp->src = -1;
			vp->type1 = ap->type;
		}
		arena->vdw_par_i[i] = (vp - par->vdwPars);
	}
	nmissing = par->nvdwPars;

	/*  Copy the vdw parameters  */
	/*  Atoms with the same vdw parameters point to one common entry  */
	/*  Except when the atom appears in one of the pair-specific vdw parameters with
	    atom index specification; in that case, that atom is given a separate entry  */
	for (i = 0, idx = par->nvdwPars; i < mol->natoms; i++) {
		t1 = arena->vdw_par_i[i];
		if (t1 < ATOMS_MAX_NUMBER)
			continue;
		arena->vdw_par_i[i] = idx;
		if (mol->par != NULL) {
			/*  Look up the pair-specific vdw parameters with atom index specification  */
			for (j = mol->par->nvdwpPars - 1; j >= 0; j--) {
				VdwPairPar *vpp = mol->par->vdwpPars + j;
				if (vpp->type1 == i || vpp->type2 == i)
					break;
			}
		} else j = -1;
		if (j < 0) {
			/*  Not found: all other entries with the same vdw parameters share the entry  */
			for (j = i + 1; j < mol->natoms; j++) {
				if (arena->vdw_par_i[j] == t1)
					arena->vdw_par_i[j] = idx;
			}
		}
		if (t1 >= ATOMS_MAX_NUMBER * 2)
			vp = mol->par->vdwPars + (t1 - ATOMS_MAX_NUMBER * 2);
		else vp = gBuiltinParameters->vdwPars + (t1 - ATOMS_MAX_NUMBER);
		AssignArray(&par->vdwPars, &par->nvdwPars, sizeof(VdwPar), idx, vp);
		idx++;
	}
	
	/*  Debug output  */
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		fprintf(arena->debug_result, "\n  Atom van der Waals parameters\n");
		fprintf(arena->debug_result, "  No. type vdw_i   sigma    eps    sigma14  eps14\n");
		for (i = 0; i < mol->natoms; i++) {
			char s[8];
			idx = arena->vdw_par_i[i];
			fprintf(arena->debug_result, "%5d %-4s %5d ", i+1, AtomTypeDecodeToString(mol->atoms[i].type, s), idx);
			if (idx >= 0) {
				Double sigma, eps, sigma2, eps2;
				sigma = par->vdwPars[idx].A;
				eps = par->vdwPars[idx].B;
				if (sigma != 0.0) {
					eps = 0.25 * eps * eps / sigma;
					sigma = pow(sigma / eps * 0.25, 1.0/12.0);
				} else {
					eps = sigma = 0.0;
				}
				sigma2 = par->vdwPars[idx].A14;
				eps2 = par->vdwPars[idx].B14;
				if (sigma2 != 0.0) {
					eps2 = 0.25 * eps2 * eps2 / sigma2;
					sigma2 = pow(sigma2 / eps2 * 0.25, 1.0/12.0);
				} else {
					eps2 = sigma2 = 0.0;
				}
				fprintf(arena->debug_result, " %7.3f %7.3f %7.3f %7.3f\n", sigma, eps*INTERNAL2KCAL, sigma2, eps2*INTERNAL2KCAL);
			} else {
				fprintf(arena->debug_result, " (not available)\n");
			}
		}
	}

	/*  Build cache  */
	if (arena->vdw_cache != NULL)
		free(arena->vdw_cache);
	arena->vdw_cache = (MDVdwCache *)calloc(sizeof(MDVdwCache), par->nvdwPars * par->nvdwPars + 1);
	if (arena->vdw_cache == NULL)
		md_panic(arena, ERROR_out_of_memory);
	cutoff6 = arena->cutoff * arena->cutoff;
	cutoff6 = 1.0 / (cutoff6 * cutoff6 * cutoff6);
	cutoff12 = cutoff6 * cutoff6;
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		fprintf(arena->debug_result, "\n  van der Waals parameter cache table\n");
		fprintf(arena->debug_result, " idx1  idx2  type1 type2  sigma    eps    sigma14  eps14\n");
	}
	for (i = 0; i < par->nvdwPars; i++) {
		t1 = arena->vdw_par_i[i];
		type1 = par->vdwPars[t1].type1;
		for (j = i; j < par->nvdwPars; j++) {
			VdwPairPar newpar, *vpp;
			Double sigma1, sigma2, eps1, eps2;
			vpp = NULL;
			t2 = arena->vdw_par_i[j];
			type2 = par->vdwPars[t2].type1;  /*  Not type2  */
			if (i != j) {
				/*  Look up the pair-specific van der Waals parameters  */
				if (mol->par != NULL) {
					vpp = ParameterLookupVdwPairPar(mol->par, type1, type2, 0);
				}
				if (vpp == NULL)
					vpp = ParameterLookupVdwPairPar(gBuiltinParameters, type1, type2, 0);
			}
			if (vpp != NULL) {
				newpar = *vpp;
			} else {
				/*  Create a pair-specific VdwPar for this pair  */
				VdwPar *p1, *p2;
				p1 = par->vdwPars + i;
				p2 = par->vdwPars + j;
				if (p1->A < 1e-15 && p1->A > -1e-15) {
					sigma1 = 1.0;
					eps1 = 0.0;
				} else {
					eps1 = 0.25 * p1->B * p1->B / p1->A;
					sigma1 = pow(p1->B * 0.25 / eps1, 1.0/6.0);
				}
				if (p2->A < 1e-15 && p2->A > -1e-15) {
					sigma2 = 1.0;
					eps2 = 0.0;
				} else {
					eps2 = 0.25 * p2->B * p2->B / p2->A;
					sigma2 = pow(p2->B * 0.25 / eps2, 1.0/6.0);
				}
				sigma1 = (sigma1 + sigma2) * 0.5;
				eps1 = sqrt(eps1 * eps2);
				sigma1 = sigma1 * sigma1 * sigma1;
				sigma1 = sigma1 * sigma1;
				newpar.B = 4.0 * sigma1 * eps1;
				newpar.A = newpar.B * sigma1;
				if (p1->A14 < 1e-15 && p1->A14 > -1e-15) {
					sigma1 = 1.0;
					eps1 = 0.0;
				} else {
					eps1 = 0.25 * p1->B14 * p1->B14 / p1->A14;
					sigma1 = pow(p1->B14 * 0.25 / eps1, 1.0/6.0);
				}
				if (p2->A14 < 1e-15 && p2->A14 > -1e-15) {
					sigma2 = 1.0;
					eps2 = 0.0;
				} else {
					eps2 = 0.25 * p2->B14 * p2->B14 / p2->A14;
					sigma2 = pow(p2->B14 * 0.25 / eps2, 1.0/6.0);
				}
				sigma1 = (sigma1 + sigma2) * 0.5;
				eps1 = sqrt(eps1 * eps2);
				sigma1 = sigma1 * sigma1 * sigma1;
				sigma1 = sigma1 * sigma1;
				newpar.B14 = 4.0 * sigma1 * eps1;
				newpar.A14 = newpar.B14 * sigma1;
				newpar.type1 = type1;
				newpar.type2 = type2;
			}
			idx = i * par->nvdwPars + j;
			arena->vdw_cache[idx].par = newpar;
			/*  Cache the value at the cutoff distance  */
			arena->vdw_cache[idx].vcut = newpar.A * cutoff12 - newpar.B * cutoff6;
			arena->vdw_cache[idx].vcut14 = newpar.A14 * cutoff12 - newpar.B * cutoff6;
			arena->vdw_cache[j * par->nvdwPars + i] = arena->vdw_cache[idx];
			if (arena->debug_result != NULL && arena->debug_output_level > 0) {
				char s1[8], s2[8];
				eps1 = 0.25 * newpar.B * newpar.B / newpar.A;
				sigma1 = pow(newpar.B * 0.25 / eps1, 1.0/6.0);
				eps2 = 0.25 * newpar.B14 * newpar.B14 / newpar.A14;
				sigma2 = pow(newpar.B14 * 0.25 / eps2, 1.0/6.0);				
				fprintf(arena->debug_result, "%5d %5d %-4s  %-4s  %7.3f %7.3f %7.3f %7.3f\n", idx, i, AtomTypeDecodeToString(type1, s1), AtomTypeDecodeToString(type2, s2), sigma1, eps1*INTERNAL2KCAL, sigma2, eps2*INTERNAL2KCAL);
			}
		} /* end loop j */
	} /* end loop i */
	arena->nmissing += nmissing;
	return nmissing;
}

/*  Find bond parameters  */
static int
s_find_bond_parameters(MDArena *arena)
{
	Int i, j, idx, type1, type2, t1, nmissing = 0;
	Molecule *mol = arena->mol;
	Parameter *par = arena->par;
	BondPar *bp;

	if (mol->nbonds > 0) {
		if (arena->bond_par_i != NULL)
			free(arena->bond_par_i);
		arena->bond_par_i = (Int *)calloc(sizeof(Int), mol->nbonds);
		if (arena->bond_par_i == NULL)
			md_panic(arena, ERROR_out_of_memory);
		if (arena->anbond_r0 != NULL)
			free(arena->anbond_r0);
		arena->anbond_r0 = (Double *)calloc(sizeof(Double), mol->nbonds);
		if (arena->anbond_r0 == NULL)
			md_panic(arena, ERROR_out_of_memory);

		/*  Find the bond parameter; priority: (1) atom index specific in mol->par, (2)
		 atom type specific in mol->par, (3) atom type specific in built-in  */
		for (i = 0; i < mol->nbonds; i++) {
			Int i1, i2;
			Atom *ap1, *ap2;
			i1 = mol->bonds[i * 2];
			i2 = mol->bonds[i * 2 + 1];
			ap1 = ATOM_AT_INDEX(mol->atoms, i1);
			ap2 = ATOM_AT_INDEX(mol->atoms, i2);
			type1 = ap1->type;
			type2 = ap2->type;
			bp = NULL;
			if (mol->par != NULL) {
				bp = ParameterLookupBondPar(mol->par, type1, type2, 0);
				if (bp != NULL) {
					arena->bond_par_i[i] = (bp - mol->par->bondPars) + ATOMS_MAX_NUMBER * 2;
					continue;
				}
			}
			bp = ParameterLookupBondPar(gBuiltinParameters, type1, type2, 0);
			if (bp != NULL) {
				arena->bond_par_i[i] = (bp - gBuiltinParameters->bondPars) + ATOMS_MAX_NUMBER;
				continue;
			}
			bp = ParameterLookupBondPar(par, type1, type2, kParameterLookupMissing | kParameterLookupNoBaseAtomType);
			if (bp == NULL) {
				/*  Record as missing  */
				bp = AssignArray(&par->bondPars, &par->nbondPars, sizeof(BondPar), par->nbondPars, NULL);
				bp->src = -1;
				bp->type1 = type1;
				bp->type2 = type2;
			}
			arena->bond_par_i[i] = (bp - par->bondPars);
		}
		nmissing = par->nbondPars;
	
		/*  Copy the bond parameters  */
		for (i = 0, idx = par->nbondPars; i < mol->nbonds; i++) {
			t1 = arena->bond_par_i[i];
			if (t1 < ATOMS_MAX_NUMBER)
				continue;
			arena->bond_par_i[i] = idx;
			for (j = i + 1; j < mol->nbonds; j++) {
				if (arena->bond_par_i[j] == t1)
					arena->bond_par_i[j] = idx;
			}
			if (t1 >= ATOMS_MAX_NUMBER * 2)
				bp = mol->par->bondPars + (t1 - ATOMS_MAX_NUMBER * 2);
			else bp = gBuiltinParameters->bondPars + (t1 - ATOMS_MAX_NUMBER);
			AssignArray(&par->bondPars, &par->nbondPars, sizeof(BondPar), idx, bp);
			idx++;
		}
	}
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		char s1[8], s2[8];
		fprintf(arena->debug_result, "\n  Bond parameters\n");
		fprintf(arena->debug_result, "  No. atom1 atom2 type1 type2     r0      k\n");
		for (i = 0; i < mol->nbonds; i++) {
			idx = arena->bond_par_i[i];
			fprintf(arena->debug_result, "%5d %5d %5d ", i+1, mol->bonds[i*2]+1, mol->bonds[i*2+1]+1);
			if (idx < 0) {
				fprintf(arena->debug_result, "<< missing >>\n");
				continue;
			}
			fprintf(arena->debug_result, " %-4s  %-4s %7.3f %7.3f\n", AtomTypeDecodeToString(par->bondPars[idx].type1, s1), AtomTypeDecodeToString(par->bondPars[idx].type2, s2), par->bondPars[idx].r0, par->bondPars[idx].k*INTERNAL2KCAL);
		}
	}
	arena->nmissing += nmissing;
	return nmissing;
}

/*  Find angle parameters  */
static int
s_find_angle_parameters(MDArena *arena)
{
	Int i, j, idx, type1, type2, type3, t1, nmissing = 0;
	Molecule *mol = arena->mol;
	Parameter *par = arena->par;
	AnglePar *ap;
	
	if (mol->nangles > 0) {
		if (arena->angle_par_i != NULL)
			free(arena->angle_par_i);
		arena->angle_par_i = (Int *)calloc(sizeof(Int), mol->nangles);
		if (arena->angle_par_i == NULL)
			md_panic(arena, ERROR_out_of_memory);
		
		/*  Find the angle parameter; priority: (1) atom index specific in mol->par, (2)
		 atom type specific in mol->par, (3) atom type specific in built-in  */
		for (i = 0; i < mol->nangles; i++) {
			Int i1, i2, i3;
			i1 = mol->angles[i * 3];
			i2 = mol->angles[i * 3 + 1];
			i3 = mol->angles[i * 3 + 2];
			type1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
			type2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
			type3 = ATOM_AT_INDEX(mol->atoms, i3)->type;
			if (mol->par != NULL) {
			/*	ap = ParameterLookupAnglePar(mol->par, i1, i2, i3, kParameterLookupLocal | kParameterLookupNoWildcard);
				if (ap == NULL)
					ap = ParameterLookupAnglePar(mol->par, type1, type2, type3, kParameterLookupLocal | kParameterLookupGlobal); */
				ap = ParameterLookupAnglePar(mol->par, type1, type2, type3, 0);
				if (ap != NULL) {
					arena->angle_par_i[i] = (ap - mol->par->anglePars) + ATOMS_MAX_NUMBER * 2;
					continue;
				}
			}
			ap = ParameterLookupAnglePar(gBuiltinParameters, type1, type2, type3, 0);
			if (ap != NULL) {
				arena->angle_par_i[i] = (ap - gBuiltinParameters->anglePars) + ATOMS_MAX_NUMBER;
				continue;
			}
			/*  Record as missing  */
			ap = ParameterLookupAnglePar(par, type1, type2, type3, kParameterLookupMissing | kParameterLookupNoBaseAtomType);
			if (ap == NULL) {
				ap = AssignArray(&par->anglePars, &par->nanglePars, sizeof(AnglePar), par->nanglePars, NULL);
				ap->src = -1;
				ap->type1 = type1;
				ap->type2 = type2;
				ap->type3 = type3;
			}
			arena->angle_par_i[i] = (ap - par->anglePars);
		}
		nmissing = par->nanglePars;
		
		/*  Copy the angle parameters  */
		for (i = 0, idx = par->nanglePars; i < mol->nangles; i++) {
			t1 = arena->angle_par_i[i];
			if (t1 < ATOMS_MAX_NUMBER)
				continue;
			arena->angle_par_i[i] = idx;
			for (j = i + 1; j < mol->nangles; j++) {
				if (arena->angle_par_i[j] == t1)
					arena->angle_par_i[j] = idx;
			}
			if (t1 >= ATOMS_MAX_NUMBER * 2)
				ap = mol->par->anglePars + (t1 - ATOMS_MAX_NUMBER * 2);
			else ap = gBuiltinParameters->anglePars + (t1 - ATOMS_MAX_NUMBER);
			AssignArray(&par->anglePars, &par->nanglePars, sizeof(AnglePar), idx, ap);
			idx++;
		}
	}
	
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		char s1[8], s2[8], s3[8];
		fprintf(arena->debug_result, "\n  Angle parameters\n");
		fprintf(arena->debug_result, "  No. atom1 atom2 atom3 type1 type2 type3      a0      k\n");
		for (i = 0; i < mol->nangles; i++) {
			idx = arena->angle_par_i[i];
			fprintf(arena->debug_result, "%5d %5d %5d %5d ", i+1, mol->angles[i*3]+1, mol->angles[i*3+1]+1, mol->angles[i*3+2]+1);
			if (idx < 0) {
				fprintf(arena->debug_result, "<< missing >>\n");
				continue;
			}
			fprintf(arena->debug_result, " %-4s  %-4s  %-4s %7.3f %7.3f\n",  AtomTypeDecodeToString(par->anglePars[idx].type1, s1), AtomTypeDecodeToString(par->anglePars[idx].type2, s2), AtomTypeDecodeToString(par->anglePars[idx].type3, s3), par->anglePars[idx].a0 * 180.0 / 3.1415927, par->anglePars[idx].k*INTERNAL2KCAL);
		}
	}
	arena->nmissing += nmissing;
	return nmissing;
}

/*  Find dihedral parameters  */
static int
s_find_dihedral_parameters(MDArena *arena)
{
	Int i, j, idx, type1, type2, type3, type4, t1, nmissing = 0;
	Molecule *mol = arena->mol;
	Parameter *par = arena->par;
	TorsionPar *tp;
	
	if (mol->ndihedrals > 0) {
		if (arena->dihedral_par_i != NULL)
			free(arena->dihedral_par_i);
		arena->dihedral_par_i = (Int *)calloc(sizeof(Int), mol->ndihedrals);
		if (arena->dihedral_par_i == NULL)
			md_panic(arena, ERROR_out_of_memory);
		
		/*  Find the dihedral parameter; priority: (1) atom index specific in mol->par, (2)
		 atom type specific in mol->par, (3) atom type specific in built-in  */
		for (i = 0; i < mol->ndihedrals; i++) {
			Int i1, i2, i3, i4;
			i1 = mol->dihedrals[i * 4];
			i2 = mol->dihedrals[i * 4 + 1];
			i3 = mol->dihedrals[i * 4 + 2];
			i4 = mol->dihedrals[i * 4 + 3];
			type1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
			type2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
			type3 = ATOM_AT_INDEX(mol->atoms, i3)->type;
			type4 = ATOM_AT_INDEX(mol->atoms, i4)->type;
			if (mol->par != NULL) {
			/*	tp = ParameterLookupDihedralPar(mol->par, i1, i2, i3, i4, kParameterLookupLocal | kParameterLookupNoWildcard);
				if (tp == NULL)
					tp = ParameterLookupDihedralPar(mol->par, type1, type2, type3, type4, kParameterLookupLocal | kParameterLookupGlobal); */
				tp = ParameterLookupDihedralPar(mol->par, type1, type2, type3, type4, 0);
				if (tp != NULL) {
					arena->dihedral_par_i[i] = (tp - mol->par->dihedralPars) + ATOMS_MAX_NUMBER * 2;
					continue;
				}
			}
			tp = ParameterLookupDihedralPar(gBuiltinParameters, type1, type2, type3, type4, 0);
			if (tp != NULL) {
				arena->dihedral_par_i[i] = (tp - gBuiltinParameters->dihedralPars) + ATOMS_MAX_NUMBER;
				continue;
			}
			/*  Record as missing  */
			tp = ParameterLookupDihedralPar(par, type1, type2, type3, type4, kParameterLookupMissing | kParameterLookupNoBaseAtomType);
			if (tp == NULL) {
				tp = AssignArray(&par->dihedralPars, &par->ndihedralPars, sizeof(TorsionPar), par->ndihedralPars, NULL);
				tp->src = -1;
				tp->type1 = type1;
				tp->type2 = type2;
				tp->type3 = type3;
				tp->type4 = type4;
			}
			arena->dihedral_par_i[i] = (tp - par->dihedralPars);
		}
		nmissing = par->ndihedralPars;
	
		/*  Copy the dihedral parameters  */
		for (i = 0, idx = par->ndihedralPars; i < mol->ndihedrals; i++) {
			t1 = arena->dihedral_par_i[i];
			if (t1 < ATOMS_MAX_NUMBER)
				continue;
			arena->dihedral_par_i[i] = idx;
			for (j = i + 1; j < mol->ndihedrals; j++) {
				if (arena->dihedral_par_i[j] == t1)
					arena->dihedral_par_i[j] = idx;
			}
			if (t1 >= ATOMS_MAX_NUMBER * 2)
				tp = mol->par->dihedralPars + (t1 - ATOMS_MAX_NUMBER * 2);
			else tp = gBuiltinParameters->dihedralPars + (t1 - ATOMS_MAX_NUMBER);
			AssignArray(&par->dihedralPars, &par->ndihedralPars, sizeof(TorsionPar), idx, tp);
			idx++;
		}
	}
	
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		Int j;
		char s1[8], s2[8], s3[8], s4[8];
		fprintf(arena->debug_result, "\n  Dihedral parameters\n");
		fprintf(arena->debug_result, "  No. atom1 atom2 atom3 atom4 type1 type2 type3 type4 mult  phi0      k     per\n");
		for (i = 0; i < mol->ndihedrals; i++) {
			idx = arena->dihedral_par_i[i];
			fprintf(arena->debug_result, "%5d %5d %5d %5d %5d ", i+1, mol->dihedrals[i*4]+1, mol->dihedrals[i*4+1]+1, mol->dihedrals[i*4+2]+1, mol->dihedrals[i*4+3]+1);
			if (idx < 0) {
				fprintf(arena->debug_result, "<< missing >>\n");
				continue;
			}
			fprintf(arena->debug_result, "%-4s  %-4s  %-4s  %-4s  %2d ", AtomTypeDecodeToString(par->dihedralPars[idx].type1, s1), AtomTypeDecodeToString(par->dihedralPars[idx].type2, s2), AtomTypeDecodeToString(par->dihedralPars[idx].type3, s3), AtomTypeDecodeToString(par->dihedralPars[idx].type4, s4), par->dihedralPars[idx].mult);
			for (j = 0; j < par->dihedralPars[idx].mult; j++) {
				fprintf(arena->debug_result, "%7.3f %7.3f %1d ", par->dihedralPars[idx].phi0[j]*180/PI, par->dihedralPars[idx].k[j]*INTERNAL2KCAL, par->dihedralPars[idx].period[j]);
			}
			fprintf(arena->debug_result, "\n");
		}
	}
	arena->nmissing += nmissing;
	return nmissing;
}

/*  Find improper parameters  */
static int
s_find_improper_parameters(MDArena *arena)
{
	Int i, j, idx, type1, type2, type3, type4, t1, nmissing = 0;
	Molecule *mol = arena->mol;
	Parameter *par = arena->par;
	TorsionPar *tp;
	
	if (mol->nimpropers > 0) {
		if (arena->improper_par_i != NULL)
			free(arena->improper_par_i);
		arena->improper_par_i = (Int *)calloc(sizeof(Int), mol->nimpropers);
		if (arena->improper_par_i == NULL)
			md_panic(arena, ERROR_out_of_memory);
		
		/*  Find the improper parameter; priority: (1) atom index specific in mol->par, (2)
		 atom type specific in mol->par, (3) atom type specific in built-in  */
		for (i = 0; i < mol->nimpropers; i++) {
			Int i1, i2, i3, i4;
			i1 = mol->impropers[i * 4];
			i2 = mol->impropers[i * 4 + 1];
			i3 = mol->impropers[i * 4 + 2];
			i4 = mol->impropers[i * 4 + 3];
			type1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
			type2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
			type3 = ATOM_AT_INDEX(mol->atoms, i3)->type;
			type4 = ATOM_AT_INDEX(mol->atoms, i4)->type;
			if (mol->par != NULL) {
			/*	tp = ParameterLookupImproperPar(mol->par, i1, i2, i3, i4, kParameterLookupLocal | kParameterLookupNoWildcard);
				if (tp == NULL)
					tp = ParameterLookupImproperPar(mol->par, type1, type2, type3, type4, kParameterLookupLocal | kParameterLookupGlobal); */
				tp = ParameterLookupImproperPar(mol->par, type1, type2, type3, type4, 0);
				if (tp != NULL) {
					arena->improper_par_i[i] = (tp - mol->par->improperPars) + ATOMS_MAX_NUMBER * 2;
					continue;
				}
			}
			tp = ParameterLookupImproperPar(gBuiltinParameters, type1, type2, type3, type4, 0);
			if (tp != NULL) {
				arena->improper_par_i[i] = (tp - gBuiltinParameters->improperPars) + ATOMS_MAX_NUMBER;
				continue;
			}
			/*  Record as missing  */
			tp = ParameterLookupImproperPar(par, type1, type2, type3, type4, kParameterLookupMissing | kParameterLookupNoBaseAtomType);
			if (tp == NULL) {
				tp = AssignArray(&par->improperPars, &par->nimproperPars, sizeof(TorsionPar), par->nimproperPars, NULL);
				tp->src = -1;
				tp->type1 = type1;
				tp->type2 = type2;
				tp->type3 = type3;
				tp->type4 = type4;
			}
			arena->improper_par_i[i] = (tp - par->improperPars);
		}
		nmissing = par->nimproperPars;
	
		/*  Copy the improper parameters  */
		for (i = 0, idx = par->nimproperPars; i < mol->nimpropers; i++) {
			t1 = arena->improper_par_i[i];
			if (t1 < ATOMS_MAX_NUMBER)
				continue;
			arena->improper_par_i[i] = idx;
			for (j = i + 1; j < mol->nimpropers; j++) {
				if (arena->improper_par_i[j] == t1)
					arena->improper_par_i[j] = idx;
			}
			if (t1 >= ATOMS_MAX_NUMBER * 2)
				tp = mol->par->improperPars + (t1 - ATOMS_MAX_NUMBER * 2);
			else tp = gBuiltinParameters->improperPars + (t1 - ATOMS_MAX_NUMBER);
			AssignArray(&par->improperPars, &par->nimproperPars, sizeof(TorsionPar), idx, tp);
			idx++;
		}
	}
	
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		Int j;
		char s1[8], s2[8], s3[8], s4[8];
		fprintf(arena->debug_result, "\n  Improper parameters\n");
		fprintf(arena->debug_result, "  No. atom1 atom2 atom3 atom4 type1 type2 type3 type4 mult  phi0      k   per\n");
		for (i = 0; i < mol->nimpropers; i++) {
			idx = arena->improper_par_i[i];
			fprintf(arena->debug_result, "%5d %5d %5d %5d %5d ", i+1, mol->impropers[i*4]+1, mol->impropers[i*4+1]+1, mol->impropers[i*4+2]+1, mol->impropers[i*4+3]+1);
			if (idx < 0) {
				fprintf(arena->debug_result, "<< missing >>\n");
				continue;
			}
			fprintf(arena->debug_result, "%-4s  %-4s  %-4s  %-4s  %2d ", AtomTypeDecodeToString(par->improperPars[idx].type1, s1), AtomTypeDecodeToString(par->improperPars[idx].type2, s2), AtomTypeDecodeToString(par->improperPars[idx].type3, s3), AtomTypeDecodeToString(par->improperPars[idx].type4, s4), par->improperPars[idx].mult);
			for (j = 0; j < par->improperPars[idx].mult; j++) {
				fprintf(arena->debug_result, "%7.3f %7.3f %1d", par->improperPars[idx].phi0[j]*180/PI, par->improperPars[idx].k[j]*INTERNAL2KCAL, par->improperPars[idx].period[j]);
			}
			fprintf(arena->debug_result, "\n");
		}
	}	
	arena->nmissing += nmissing;
	return nmissing;
}

/*  Find one fragment, starting from start_index  */
static void
s_find_fragment_sub(MDArena *arena, Int start_index, Int fragment_index)
{
	Atom *atoms = arena->mol->atoms;
	int i, j;
	for (i = 0; i < atoms[start_index].nconnects; i++) {
		j = atoms[start_index].connects[i];
		if (j >= 0 && j < arena->natoms_uniq) {
			int n = arena->fragment_indices[j];
			if (n < 0) {
				arena->fragment_indices[j] = fragment_index;
				s_find_fragment_sub(arena, j, fragment_index);
			} else if (n != fragment_index) {
				if (arena->log_result != NULL)
					fprintf(arena->log_result, "Warning: internal inconsistency in finding fragment at atom %d\n", j + 1);
			}
		}
	}
}

/*  Find fragments  */
void
md_find_fragments(MDArena *arena)
{
	int i, idx, nuniq;
	if (arena->fragment_indices != NULL)
		free(arena->fragment_indices);
	arena->fragment_indices = (Int *)malloc(sizeof(Int) * arena->natoms_uniq);
	if (arena->fragment_indices == NULL)
		md_panic(arena, "Low memory in md_find_fragments");
	nuniq = arena->natoms_uniq;
	for (i = 0; i < nuniq; i++)
		arena->fragment_indices[i] = -1;
	idx = 0;
	for (i = 0; i < nuniq; i++) {
		if (arena->fragment_indices[i] >= 0)
			continue;
		s_find_fragment_sub(arena, i, idx);
		idx++;
	}
	arena->nfragments = idx;
	if (arena->fragment_info != NULL)
		free(arena->fragment_info);
	arena->fragment_info = (struct MDFragmentInfo *)calloc(sizeof(struct MDFragmentInfo), idx);
	if (arena->fragment_info == NULL)
		md_panic(arena, "Low memory in md_find_fragments");
}

/*  Calculate center-of-mass  */
void
md_center_of_mass(MDArena *arena, Vector *cp)
{
	int i, n;
	Atom *ap;
	Double sumw;
	VecZero(*cp);
	sumw = 0.0;
	n = arena->mol->natoms;
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		VecScaleInc(*cp, ap->r, ap->weight);
		sumw += ap->weight;
	}
	VecScaleSelf(*cp, 1.0 / sumw);
}

/*  Relocate center-of-mass  */
int
md_relocate_center(MDArena *arena)
{
	int i, n;
	Atom *ap;
	Vector cm;
	n = arena->mol->natoms;
	md_center_of_mass(arena, &cm);
	VecDec(cm, arena->initial_center);
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		VecDec(ap->r, cm);
	}
	return 0;
}

/*  Quench translational and rotational momenta  */
int
md_quench_momenta(MDArena *arena)
{
	int i, n;
	Atom *ap;
	Vector am, cm, m, v, v1, v2;
	Double w, sumw, kinetic;
	Mat33 mi;
	
	if (!arena->quench_angular_momentum && !arena->quench_translational_momentum)
		return 0;

	n = arena->mol->natoms;

	/*  Center of mass  */
	md_center_of_mass(arena, &cm);
	
	/*  Initial kinetic energy (times 2)  */
	kinetic = 0.0;
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		kinetic += ap->weight * VecLength2(ap->v);
	}

	/*  Quench the angular momentum  */
	/*  Calculate the total angular momentum  */
	
	if (arena->quench_angular_momentum) {
		VecZero(am);
		for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
			VecSub(v, ap->r, cm);
			VecCross(v1, v, ap->v);
			VecScaleInc(am, v1, ap->weight);
		}
		/*  Moment of inertia  */
		memset(mi, 0, sizeof(mi));
		for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
			Double x = ap->r.x - cm.x;
			Double y = ap->r.y - cm.y;
			Double z = ap->r.z - cm.z;
			w = ap->weight;
			mi[0] += w * (y * y + z * z);
			mi[4] += w * (z * z + x * x);
			mi[8] += w * (x * x + y * y);
			mi[1] -= w * x * y;
			mi[5] -= w * y * z;
			mi[2] -= w * z * x;
		}
		mi[3] = mi[1];
		mi[6] = mi[2];
		mi[7] = mi[5];
		/*  Calculate the angular velocity as a rigid body  */
		/*  I * w = L; I, moment of inertia; w, angular velocity; L, angular momentum  */
		MatrixInvert(mi, mi);
		MatrixVec(&v, mi, &am);
		/*  Subtract the velocity at the atom position derived from the angular velocity  */
		for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
			VecSub(v1, ap->r, cm);
			VecCross(v2, v, v1);
			VecDec(ap->v, v2);
		}
	}
	
	/*  Quench the translational momentum  */
	if (arena->quench_translational_momentum) {
		VecZero(m);
		sumw = 0.0;
		for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
			VecScaleInc(m, ap->v, ap->weight);
			sumw += ap->weight;
		}
		VecScaleSelf(m, 1.0 / sumw);
		for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++)
			VecDec(ap->v, m);
	}

	/*  Current kinetic energy (times 2) */
	w = 0.0;
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		w += ap->weight * VecLength2(ap->v);
	}
	w = sqrt(kinetic / w);
	
	/*  Scale the velocities to keep the kinetic energy  */
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		VecScaleSelf(ap->v, w);
	}
	return 0;
}

/*  Initilize the runtime fields that are affected by the atom positions  */
/*  (Should be called after coordinates are updated without changing structure info)  */
void
md_init_for_positions(MDArena *arena)
{
	Molecule *mol;
	Parameter *par;
	int i, j, idx;

	if (arena == NULL || (mol = arena->mol) == NULL || (par = arena->par) == NULL)
		return;
	
	/*  Initialize fix positions  */
	for (i = j = 0; i < mol->natoms; i++) {
		Atom *ap = &(mol->atoms[i]);
		if (ap->fix_force != 0.0) {
			ap->fix_pos = ap->r;
			j++;
		}
	}
	if (j > 0)
		md_log(arena, "Number of fixed atoms = %d\n", j);
/*	if (arena->fix_atoms != NULL) {
		for (i = 0; i < arena->nfix_atoms; i++) {
			j = arena->fix_atoms[i].index;
			if (j >= 0 && j < mol->natoms)
				arena->fix_atoms[i].pos = mol->atoms[j].r;
		}
		printf("Number of fixed atoms = %d\n", arena->nfix_atoms);
	} */

	/*  Abnormal bonds  */
	if (arena->anbond_thres > 0.0) {
		Vector r12;
		Double r;
		for (i = 0; i < mol->nbonds; i++) {
			idx = arena->bond_par_i[i];
			if (idx < 0)
				continue;
			VecSub(r12, mol->atoms[mol->bonds[i*2]].r, mol->atoms[mol->bonds[i*2+1]].r);
			r = VecLength(r12);
			if (r >= (1 + arena->anbond_thres) * par->bondPars[idx].r0)
				arena->anbond_r0[i] = r - par->bondPars[idx].r0;
			else
				arena->anbond_r0[i] = 0.0;
		}
	}
	
	/*  Center of mass  */
	md_center_of_mass(arena, &(arena->initial_center));
}

const char *
md_prepare(MDArena *arena, int check_only)
{
	Int i, idx;
	Int t1, t2;
	Molecule *mol;

	if (arena->xmol->needsMDRebuild) {
		/*  Set molecule again  */
		md_arena_set_molecule(arena, arena->xmol);
		arena->xmol->needsMDRebuild = 0;
	} else {
		md_copy_coordinates_to_internal(arena);
	}
	mol = arena->mol;

	arena->errmsg[0] = 0;
	arena->setjmp_buf = NULL;

	/*  Table of missing parameters  */
/*	Int *missing = NULL; */
/*	Int nmissing = 0;  */
	
	if (mol == NULL)
		return "molecule is not defined";
	if (mol->natoms == 0 || mol->atoms == NULL)
		return "molecule has no atoms";

	/*  Open files  */
	if (!check_only) {
		char *cwd = NULL;
		const char *err = NULL;
		if (mol->path != NULL) {
			/*  Temporarily change to the document directory  */
			char *p;
			char *fname = strdup(mol->path);
			if ((p = strrchr(fname, '/')) != NULL
				|| (p = strrchr(fname, '\\')) != NULL
				) {
				*p = 0;
				cwd = getcwd(NULL, 0);
				chdir(fname);
			}
			free(fname);
		}
		if (arena->log_result_name != NULL && arena->log_result == NULL) {
			arena->log_result = fopen(arena->log_result_name, "wb");
			if (arena->log_result == NULL)
				err = "cannot create log file";
		}
		if (err == NULL && arena->coord_result_name != NULL && arena->coord_result == NULL) {
			arena->coord_result = fopen(arena->coord_result_name, "wb");
			if (arena->coord_result == NULL)
				err = "cannot create coord file";
		}
		if (err == NULL && arena->vel_result_name != NULL && arena->vel_result == NULL) {
			arena->vel_result = fopen(arena->vel_result_name, "wb");
			if (arena->vel_result == NULL)
				err = "cannot create vel file";
		}
		if (err == NULL && arena->force_result_name != NULL && arena->force_result == NULL) {
			arena->force_result = fopen(arena->force_result_name, "wb");
			if (arena->force_result == NULL)
				err = "cannot create force file";
		}
		if (err == NULL && arena->extend_result_name != NULL && arena->extend_result == NULL) {
			arena->extend_result = fopen(arena->extend_result_name, "wb");
			if (arena->extend_result == NULL)
				err = "cannot create extend file";
		}
		if (err == NULL && arena->debug_result_name != NULL && arena->debug_result == NULL) {
			arena->debug_result = fopen(arena->debug_result_name, "wb");
			if (arena->debug_result == NULL)
				err = "cannot create debug file";
		}
		if (cwd != NULL) {
			chdir(cwd);
			free(cwd);
		}
		if (err != NULL)
			return err;
	}
	
	/*  Count symmetry unique atoms  */
	/*  Atoms should be ordered so that all symmetry related atoms come after symmetry unique atoms  */
	for (i = t1 = 0, t2 = -1; i < mol->natoms; i++) {
		Atom *ap = mol->atoms + i;
		Symop symop = ap->symop;
		if (SYMOP_ALIVE(symop)) {
			if (t2 == -1)
				t2 = i;  /*  The index of the first non-unique atom  */
		} else {
			t1++;        /*  The number of unique atoms  */
		}
	}
	if (t2 >= 0 && t1 != t2)
		return "all symmetry related atoms should be after symmetry unique atoms";
/*	if (t1 > mol->natoms && mol->box == NULL)
		return "symmetry operation is used but no periodic box is defined"; */

	if (mol->cell != NULL) {
		arena->periodic_a = (mol->cell->flags[0] != 0);
		arena->periodic_b = (mol->cell->flags[1] != 0);
		arena->periodic_c = (mol->cell->flags[2] != 0);
	}

	arena->natoms_uniq = t1;
	
	if (arena->debug_result != NULL) {
		time_t loc_time;
		time(&loc_time);
		fprintf(arena->debug_result, "---- LWMD Started at %s", ctime(&loc_time));
	}

	/*  Recalc angle/dihedral/improper tables from the bond table and parameters  */
/*	s_rebuild_angles(arena);
	s_rebuild_dihedrals(arena);
	s_rebuild_impropers(arena); */

	/*  Statistics of the molecule  */
	md_log(arena,
		   "Number of bonds = %d\n"
		   "Number of angles = %d\n"
		   "Number of dihedrals = %d\n"
		   "Number of impropers = %d\n", mol->nbonds, mol->nangles, mol->ndihedrals, mol->nimpropers);
	
	t1 = t2 = 0;
	for (i = 0; i < mol->natoms; i++) {
		if (mol->atoms[i].fix_force > 0)
			t1++;
		if (mol->atoms[i].fix_force < 0)
			t2++;
	}
	if (t1 > 0)
		md_log(arena, "Number of constrained atoms = %d\n", t1);
	if (t2 > 0)
		md_log(arena, "Number of fixed atoms = %d\n", t2);


	if (arena->natoms_uniq < mol->natoms) {
		md_log(arena, "Number of symmetry-unique atoms = %d\n", arena->natoms_uniq);
		t2 = 0;
		for (i = 0; i < arena->natoms_uniq; i++) {
			if (mol->atoms[i].fix_force < 0)
				t2++;
		}
	}
		
	arena->degree_of_freedom = 3 * (arena->natoms_uniq - t2);
	md_log(arena, "Degree of freedom = %d\n", arena->degree_of_freedom);

	/*  Build local cache of the used parameters  */
	if (arena->par != NULL)
		ParameterRelease(arena->par);
	arena->par = ParameterNew();
	arena->nmissing = arena->nsuspicious = 0;
	s_find_vdw_parameters(arena);
	s_find_bond_parameters(arena);
	s_find_angle_parameters(arena);
	s_find_dihedral_parameters(arena);
	s_find_improper_parameters(arena);
	
	if (arena->nmissing > 0) {
	/*	for (i = 0; i < nmissing; i++) {
			char s1[6], s2[6], s3[6], s4[6];
			type1 = missing[i * 5];
			AtomTypeDecodeToString(missing[i * 5 + 1], s1);
			AtomTypeDecodeToString(missing[i * 5 + 2], s2);
			AtomTypeDecodeToString(missing[i * 5 + 3], s3);
			AtomTypeDecodeToString(missing[i * 5 + 4], s4);
			if (type1 == 0)
				md_warning(arena, "Missing vdw parameter for %s\n", s1, s2);
			else if (type1 == 1)
				md_warning(arena, "Missing bond parameter for %s-%s\n", s1, s2);
			else if (type1 == 2)
				md_warning(arena, "Missing angle parameter for %s-%s-%s\n", s1, s2, s3);
			else if (type1 == 3)
				md_warning(arena, "Missing dihedral parameter for %s-%s-%s-%s\n", s1, s2, s3, s4);
			else if (type1 == 4)
				md_warning(arena, "Missing improper parameter for %s-%s-%s-%s\n", s1, s2, s3, s4);
		}
		free(missing); */
		return "some parameters are missing";
	}
	
	/*  Build the exclusion table  */
	s_make_exclusion_list(arena);
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		fprintf(arena->debug_result, "\n  MDExclusion table for each atom\n");
		fprintf(arena->debug_result, "  No.  {self and special} {1-2} {1-3} {1-4}\n");
		for (i = 0; i < mol->natoms; i++) {
			fprintf(arena->debug_result, "%5d ", i+1);
			for (idx = arena->exinfo[i].index0; idx <= arena->exinfo[i+1].index0; idx++) {
				int n = 0;
				if (idx == arena->exinfo[i].index0)
					n += fprintf(arena->debug_result, "{");
				if (idx == arena->exinfo[i].index1)
					n += fprintf(arena->debug_result, "} {");
				if (idx == arena->exinfo[i].index2)
					n += fprintf(arena->debug_result, "} {");
				if (idx == arena->exinfo[i].index3)
					n += fprintf(arena->debug_result, "} {");
				if (idx < arena->exinfo[i+1].index0) {
					if (n == 0)
						fprintf(arena->debug_result, " ");
					fprintf(arena->debug_result, "%d", arena->exlist[idx]+1);
				}
			}
			fprintf(arena->debug_result, "}\n");
		}
	}

	/*  Parameter checking only  */
	if (check_only) {
		arena->is_initialized = 1;  /*  Only static fields are ready  */
		arena->mol->needsMDRebuild = 0;
		return NULL;
	}
	
	/*  Allocate storage for Verlet list  */
	arena->max_nverlets = mol->natoms;
	if (arena->verlets != NULL)
		free(arena->verlets);
	arena->verlets = (MDVerlet *)calloc(sizeof(MDVerlet), arena->max_nverlets);
	arena->nverlets = 0;
	if (arena->verlets_dr != NULL)
		free(arena->verlets_dr);
	arena->verlets_dr = (Vector *)calloc(sizeof(Vector), mol->natoms);
	if (arena->verlets == NULL || arena->verlets_dr == NULL)
		md_panic(arena, ERROR_out_of_memory);
	arena->last_verlet_step = -1;
	if (arena->verlet_i != NULL)
		free(arena->verlet_i);
	arena->verlet_i = (Int *)calloc(sizeof(Int), mol->natoms + 1);
	
	/*  Allocate storage for partial energy/force  */
	if (arena->energies != NULL)
		free(arena->energies);
	arena->energies = (Double *)calloc(sizeof(Double), kEndIndex);
	if (arena->forces != NULL)
		free(arena->forces);
	arena->forces = (Vector *)calloc(sizeof(Vector), kKineticIndex * mol->natoms);
	if (arena->energies == NULL || arena->forces == NULL)
		md_panic(arena, ERROR_out_of_memory);

	/*  Initialize storage for pair interaction calculations */
	/*  (The storage will be allocated when necessary)  */
/*	if (arena->group_flags != NULL) {
		free(arena->group_flags);
		arena->group_flags = NULL;
	} */
	if (arena->pair_forces != NULL) {
		free(arena->pair_forces);
		arena->pair_forces = NULL;
	}

	/*  Allocate ring buffer   */
	if (arena->ring != NULL) {
		free(arena->ring->buf);
		free(arena->ring);
	}
	arena->ring = (MDRing *)calloc(sizeof(MDRing), 1);
	if (arena->ring == NULL)
		md_panic(arena, ERROR_out_of_memory);
	arena->ring->size = mol->natoms;
	if (arena->pressure != NULL && arena->pressure->disabled == 0)
		arena->ring->size += 4;
	arena->ring->nframes = 2000 / arena->ring->size;
	if (arena->ring->nframes < 2)
		arena->ring->nframes = 2;
	arena->ring->buf = (Vector *)calloc(sizeof(Vector), arena->ring->size * arena->ring->nframes);
	if (arena->ring->buf == NULL)
		md_panic(arena, ERROR_out_of_memory);
	arena->ring->next = 0;
	arena->ring->count = 0;

	/*  Initialize temperature statistics  */
	arena->sum_temperature = 0.0;
	arena->nsum_temperature = 0;
	
	/*  Initialize unique bond/angle/dihedral/improper table  */
/*	s_find_symmetry_unique(arena); */

	/*  Initialize the position-dependent fields  */
	md_init_for_positions(arena);

	/*  Clear the snapshot storage  */
/*	if (arena->snapshots != NULL) {
		for (i = 0; i < arena->nsnapshots; i++) {
			if (arena->snapshots[i] != NULL)
				free(arena->snapshots[i]);
		}
		free(arena->snapshots);
		arena->snapshots = NULL;
	}
	arena->nsnapshots = 0;
*/
	
	/*  Random number seed  */
	{
		unsigned int seed = arena->random_seed;
		if (seed == 0)
			seed = (unsigned int)time(NULL);
		md_srand(seed);
		md_log(arena, "Random number seed = %u\n", seed);
	}
	
	/*  Clear the surface area record (will be reallocated within calc_surface_force() if necessary) */
	if (arena->sp_arena != NULL) {
		clear_sp_arena(arena);
	}
/*	if (arena->sp2_arena != NULL) {
		clear_sp2_arena(arena);
	}
 */

	/*  Graphite potential  */
	if (arena->use_graphite) {
		if (arena->graphite != NULL)
			graphite_release(arena->graphite);
		arena->graphite = graphite_new();
	}

	if (arena->velocities_read == 0)
		md_init_velocities(arena);

	/*  Fragment analysis  */
	md_find_fragments(arena);
	md_log(arena, "%d fragments found\n", arena->nfragments);
	
	/*  Pressure control statistics  */
	if (arena->pressure != NULL) {
		if (mol->cell == NULL)
			return "pressure control is requested but no unit cell is defined";
		pressure_prepare(arena);
	}

	arena->is_initialized = 2;   /*  Runtime fields are ready  */
	arena->mol->needsMDRebuild = 0;
	arena->request_abort = 0;

	return NULL;
}

#if __WXMSW__
#define random rand
#define srandom srand
#define kRandMax ((double)RAND_MAX)
#else
#define kRandMax 2147483648.0
#endif

/*  A uniform random number in [0,1]  */
Double
md_rand(void)
{
	return (double)random() / kRandMax;
}

/*  A random number with gaussian distribution  */
Double
md_gaussian_rand(void)
{
	static int waiting = 0;
	static Double next_value = 0;
	Double f, r, v1, v2;
	if (waiting) {
		waiting = 0;
		return next_value;
	}
	r = 2.0;
	while (r >= 1.0 || r < 1.523e-8) {
		v1 = 2.0 * md_rand() - 1.0;
		v2 = 2.0 * md_rand() - 1.0;
		r = v1 * v1 + v2 * v2;
	}
	f = sqrt(-2.0 * log(r) / r);
	waiting = 1;
	next_value = v1 * f;
	return v2 * f;
}

/*  Seed the random number  */
void
md_srand(unsigned int seed)
{
	srandom(seed);
}

/*  Scale velocities to match the given temperature  */
void
md_scale_velocities(MDArena *arena)
{
	int i;
	Atom *ap;
	Double ttemp, scale;
	if (arena->nsum_temperature == 0) {
		Double kinetic = 0.0;
		for (i = 0, ap = arena->mol->atoms; i < arena->natoms_uniq; i++, ap++) {
			if (ap->fix_force < 0)
				continue;
			kinetic += ap->weight * VecLength2(ap->v);
		}
		kinetic *= 0.5;
		ttemp = 2.0 * kinetic / (arena->degree_of_freedom * BOLTZMANN);
	} else {
		ttemp = arena->sum_temperature / arena->nsum_temperature;
	}
	scale = sqrt(arena->temperature / ttemp);
	for (i = 0, ap = arena->mol->atoms; i < arena->mol->natoms; i++, ap++) {
		if (ap->fix_force < 0)
			continue;
		VecScaleSelf(ap->v, scale);
	}
	arena->sum_temperature = 0;
	arena->nsum_temperature = 0;
}

/*  Give random velocities that matches the given temperature  */
void
md_init_velocities(MDArena *arena)
{
	int i, n;
	Double w;
/*	Double temp = arena->temperature; */
	Atom *ap = arena->mol->atoms;
	n = arena->mol->natoms;
	for (i = 0; i < n; i++, ap++) {
		if (ap->fix_force < 0 || fabs(ap->weight) < 1e-6) {
			ap->v.x = ap->v.y = ap->v.z = 0;
		} else {
			w = sqrt(arena->temperature * BOLTZMANN / ap->weight);
			ap->v.x = w * md_gaussian_rand();
			ap->v.y = w * md_gaussian_rand();
			ap->v.z = w * md_gaussian_rand();
		/*	ap->v.x = md_rand() - 0.5;
			ap->v.y = md_rand() - 0.5;
			ap->v.z = md_rand() - 0.5; */
		}
	}
	arena->sum_temperature = 0;
	arena->nsum_temperature = 0;

	/*  Quench the total momentum  */
	md_quench_momenta(arena);

	/*  Adjust the temperature  */
	md_scale_velocities(arena);
	
	arena->velocities_read = 0;
}

void
md_bootstrap(MDArena *arena)
{
	/*  Set the initial velocities and calculate the current force  */
	md_init_velocities(arena);
	calc_force(arena);
}

static int
s_md_output_extend_info(MDArena *arena, FILE *fp)
{
	Vector *ap, *bp, *cp, *op;
	static Vector zero = {0.0, 0.0, 0.0};
	if (arena->coord_result_frame % 10 == 0 || ftello(fp) == 0) {
		fprintf(fp, "EXLABEL: %7s %7s %11s ", "FRAME", "STEP", "POT_ENERGY");
		fprintf(fp, "%7s %7s %7s ", "CELL_AX", "CELL_AY", "CELL_AZ");
		fprintf(fp, "%7s %7s %7s ", "CELL_BX", "CELL_BY", "CELL_BZ");
		fprintf(fp, "%7s %7s %7s ", "CELL_CX", "CELL_CY", "CELL_CZ");
		fprintf(fp, "%7s %7s %7s\n", "CELL_OX", "CELL_OY", "CELL_OZ");
	}
		
	fprintf(fp, "EXINFO:  %7d %7d %11.5f ", arena->coord_result_frame, arena->step, arena->total_energy * INTERNAL2KCAL);
	if (arena->mol->cell != NULL) {
		ap = &(arena->mol->cell->axes[0]);
		bp = &(arena->mol->cell->axes[1]);
		cp = &(arena->mol->cell->axes[2]);
		op = &(arena->mol->cell->origin);
	} else {
		ap = bp = cp = op = &zero;
	}
	fprintf(fp, "%7.3f %7.3f %7.3f ", ap->x, ap->y, ap->z);
	fprintf(fp, "%7.3f %7.3f %7.3f ", bp->x, bp->y, bp->z);
	fprintf(fp, "%7.3f %7.3f %7.3f ", cp->x, cp->y, cp->z);
	fprintf(fp, "%7.3f %7.3f %7.3f\n", op->x, op->y, op->z);

	fflush(fp);
	
	return 0;
}

int
md_output_results(MDArena *arena)
{
	int i, j, natoms;
	Atom *ap;
	Vector *av, *bv, *cv;

	natoms = (arena->output_expanded_atoms ? arena->mol->natoms : arena->natoms_uniq);

#if 0
	if (arena->coord_result == NULL && arena->coord_result_name != NULL) {
		arena->coord_result = fopen(arena->coord_result_name, "wb");
		if (arena->coord_result == NULL) {
			snprintf(sErrBuf, sizeof sErrBuf, "Cannot open coordinate result file %s", arena->coord_result_name);
			arena->errmsg = sErrBuf;
			return 1;
		}
		fprintf(arena->coord_result, "TITLE: coordinates for %d atoms\n", natoms);
		arena->coord_result_frame = 0;
	}
	if (arena->vel_result == NULL && arena->vel_result_name != NULL) {
		arena->vel_result = fopen(arena->vel_result_name, "wb");
		if (arena->vel_result == NULL) {
			snprintf(sErrBuf, sizeof sErrBuf, "Cannot open velocity result file %s", arena->vel_result_name);
			arena->errmsg = sErrBuf;
			return 1;
		}
		fprintf(arena->vel_result, "TITLE: velocities for %d atoms\n", natoms);
	}
	if (arena->force_result == NULL && arena->force_result_name != NULL) {
		arena->force_result = fopen(arena->force_result_name, "wb");
		if (arena->force_result == NULL) {
			snprintf(sErrBuf, sizeof sErrBuf, "Cannot open force result file %s", arena->force_result_name);
			arena->errmsg = sErrBuf;
			return 1;
		}
		fprintf(arena->force_result, "TITLE: force for %d atoms\n", natoms);
	}
	if (arena->extend_result == NULL && arena->extend_result_name != NULL) {
		arena->extend_result = fopen(arena->extend_result_name, "wb");
		if (arena->extend_result == NULL) {
			snprintf(sErrBuf, sizeof sErrBuf, "Cannot open extend info file %s", arena->extend_result_name);
			arena->errmsg = sErrBuf;
			return 1;
		}
		fprintf(arena->extend_result, "# Frame step energy ax ay az bx by bz cx cy cz ox oy oz\n");
	}
#endif

	ap = arena->mol->atoms;
	if (arena->wrap_coordinates && arena->mol->cell != NULL) {
		av = &(arena->mol->cell->axes[0]);
		bv = &(arena->mol->cell->axes[1]);
		cv = &(arena->mol->cell->axes[2]);
	} else {
		av = bv = cv = NULL;
	}
	if (arena->coord_result != NULL) {
		if (arena->wrap_coordinates)
			md_wrap_coordinates(arena);
		for (i = j = 0; i < natoms; i++, j = (j + 3) % 10) {
			Vector r = ap[i].r;
			if (av != NULL) {
				VecScaleInc(r, *av, ap[i].wrap_dx);
				VecScaleInc(r, *bv, ap[i].wrap_dy);
				VecScaleInc(r, *cv, ap[i].wrap_dz);
			}
			fprintf(arena->coord_result, " %7.3f%s %7.3f%s %7.3f%s",
				r.x, (j == 9 ? "\n" : ""),
				r.y, (j == 8 ? "\n" : ""),
				r.z, (j == 7 ? "\n" : ""));
		}
		if (j != 0)
			fprintf(arena->coord_result, "\n");
		fflush(arena->coord_result);
	}
	if (arena->vel_result != NULL) {
		for (i = j = 0; i < natoms; i++, j = (j + 3) % 10) {
			fprintf(arena->vel_result, " %7.3f%s %7.3f%s %7.3f%s",
				ap[i].v.x*1000, (j == 9 ? "\n" : ""),
				ap[i].v.y*1000, (j == 8 ? "\n" : ""),
				ap[i].v.z*1000, (j == 7 ? "\n" : ""));
		}
		if (j != 0)
			fprintf(arena->vel_result, "\n");
		fflush(arena->vel_result);
	}
	if (arena->force_result != NULL) {
		for (i = j = 0; i < natoms; i++, j = (j + 3) % 10) {
			fprintf(arena->force_result, " %7.3f%s %7.3f%s %7.3f%s",
				ap[i].f.x, (j == 9 ? "\n" : ""),
				ap[i].f.y, (j == 8 ? "\n" : ""),
				ap[i].f.z, (j == 7 ? "\n" : ""));
		}
		if (j != 0)
			fprintf(arena->force_result, "\n");
		fflush(arena->force_result);
	}
	if (arena->extend_result != NULL) {
		s_md_output_extend_info(arena, arena->extend_result);
	} else if (arena->log_result != NULL) {
		s_md_output_extend_info(arena, arena->log_result);
	}
	arena->coord_result_frame++;
	return 0;
}

void
md_output_energies(MDArena *arena)
{
	int i;
	int periodic = (arena->mol->cell != NULL) && (arena->periodic_a || arena->periodic_b || arena->periodic_c);
/*	if (arena->log_result == NULL)
		return; */
	if ((arena->step / arena->energy_output_freq) % 10 == 0) {
		md_log(arena, "ELABEL:  %11s %11s %11s %11s %11s %11s", "STEP", "TOTAL_POT", "BOND", "ANGLE", "DIHEDRAL", "IMPROPER");
		md_log(arena, " %11s %11s %11s %11s %11s %11s %11s %11s", "VDW", "ELECT", "AUX", "SURFACE", "KINETIC", "NET", "TEMP", "TEMP_AVG");
		if (periodic)
			md_log(arena, " %11s", "VOLUME");
		md_log(arena, "\n");
	}
	md_log(arena, "ENERGY:  %11d %11.5f", arena->step, arena->total_energy * INTERNAL2KCAL);
	for (i = 0; i < kEndIndex; i++)
		md_log(arena, " %11.5f", arena->energies[i] * INTERNAL2KCAL);
	md_log(arena, " %11.5f", (arena->energies[kKineticIndex] + arena->total_energy) * INTERNAL2KCAL);
	md_log(arena, " %11.5f", arena->transient_temperature);
	md_log(arena, " %11.5f", (arena->nsum_temperature > 0 ? arena->sum_temperature / arena->nsum_temperature : 0.0));
	if (periodic) {
		Vector v, *av, *bv, *cv;
		av = &(arena->mol->cell->axes[0]);
		bv = &(arena->mol->cell->axes[1]);
		cv = &(arena->mol->cell->axes[2]);
		VecCross(v, *av, *bv);
		md_log(arena, " %11.5f", VecDot(v, *cv));
	}
	md_log(arena, "\n");
	md_log(arena, NULL);
}

void
md_update_velocities(MDArena *arena)
{
	int i, natoms;
	Byte use_sym;
	Atom *ap;
	Double w, wt;
/*	Double wsum; */
/*	Vector m1, m2; */
	Double halftimestep = arena->timestep * 0.5;
	ap = arena->mol->atoms;
	natoms = arena->mol->natoms;
	use_sym = (arena->mol->nsyms > 0);
	Double limit = arena->velocity_limit;
	limit *= limit;
/*	VecZero(m1);
	VecZero(m2);
	wsum = 0.0; */
	for (i = 0; i < natoms; i++, ap++) {
		if (use_sym && ap->symop.alive)
			continue;
		if (ap->fix_force < 0)
			continue;
		wt = ap->weight;
		if (fabs(wt) < 1e-6)
			continue;
		w = halftimestep / wt;
	/*	VecScaleInc(m1, ap->v, wt); */
		VecScaleInc(ap->v, ap->f, w);
		w = VecLength2(ap->v);
		if (w > limit || !isfinite(w))  /*  Is isfinite() available in other platforms?? */
			md_panic(arena, "the velocity for atom %d exceeded limit at step %d", i+1, arena->step);
	/*	VecScaleInc(m2, ap->v, wt); */
	/*	wsum += wt; */
	}
/*
	VecDec(m2, m1);
	w = 1.0 / wsum;
	VecScaleSelf(m2, w);
	for (i = 0, ap = arena->mol->atoms; i < natoms; i++, ap++) {
		if (use_sym && ap->symop.alive)
			continue;
		if (ap->fix_force < 0)
			continue;
		VecDec(ap->v, m2);
	}
*/
	arena->velocities_read = 0;
	/*  For debug  */
/*	VecZero(m2);
	for (i = 0, ap = arena->mol->atoms; i < natoms; i++, ap++) {
		if (use_sym && ap->sym_op > 0)
			continue;
		wt = ap->weight;
		VecScaleInc(m2, ap->v, wt);
	} */
}

void
md_update_positions(MDArena *arena)
{
	int i, natoms;
	Atom *ap;
	Double timestep = arena->timestep;
	Vector *vdr = arena->verlets_dr;
	Vector dr;
/*	Double w, limit;
	Double kinetic, kinetic_uniq; */
	Byte use_sym;
	ap = arena->mol->atoms;
	natoms = arena->mol->natoms;
/*	limit = arena->velocity_limit;
	limit *= limit; */
	use_sym = (arena->mol->nsyms > 0);
/*	kinetic = kinetic_uniq = 0.0; */
	for (i = 0; i < natoms; i++, ap++, vdr++) {
		if (use_sym && ap->symop.alive)
			continue;
		if (ap->fix_force < 0)
			continue;
		VecScale(dr, ap->v, timestep);
		VecInc(ap->r, dr);
		VecInc(*vdr, dr);
	}

	/*  Update the abnormal bond parameters  */
	if (arena->anbond_thres > 0.0) {
		Double *fp = arena->anbond_r0;
		for (i = 0; i < arena->mol->nbonds; i++) {
			if (fp[i] <= 0.0)
				continue;
			fp[i] -= arena->anbond_anneal_rate;
			if (fp[i] < 0.0)
				fp[i] = 0.0;
		}
	}
	
	/*  Relocate the center of mass (if necessary)  */
	if (arena->relocate_center) {
		md_relocate_center(arena);
	}
}

void
md_calc_kinetic_energy(MDArena *arena)
{
	int i, n, nuniq;
	Double w, kinetic, kinetic_uniq;
	Atom *ap;
	n = arena->mol->natoms;
	nuniq = arena->natoms_uniq;
	kinetic = kinetic_uniq = 0;
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		if (ap->fix_force < 0)
			continue;
		w = VecLength2(ap->v) * ap->weight;
		kinetic += w;
		if (i < nuniq)
			kinetic_uniq += w; 
	}
	arena->energies[kKineticIndex] = kinetic * 0.5;
	arena->transient_temperature = kinetic_uniq / (arena->degree_of_freedom * BOLTZMANN);
	arena->sum_temperature += arena->transient_temperature;
	arena->nsum_temperature++;
	arena->average_temperature = arena->sum_temperature / arena->nsum_temperature;
}

/*  Andersen thermostat  */
void
md_andersen_thermostat(MDArena *arena)
{
	int n = arena->natoms_uniq;
	int i;
	Double w, de;
	Atom *ap;
	Double q = arena->andersen_thermo_coupling;
	de = 0.0;
	for (i = 0, ap = arena->mol->atoms; i < n; i++, ap++) {
		if (ap->fix_force < 0 || fabs(ap->weight) < 1e-6) {
			ap->v.x = ap->v.y = ap->v.z = 0;
		} else if (md_rand() < q) {
			de -= VecLength2(ap->v) * ap->weight;
			w = sqrt(arena->temperature * BOLTZMANN / ap->weight);
			ap->v.x = w * md_gaussian_rand();
			ap->v.y = w * md_gaussian_rand();
			ap->v.z = w * md_gaussian_rand();
			de += VecLength2(ap->v) * ap->weight;
		}
	}
	de *= 0.5;
	arena->energies[kKineticIndex] += de;
}

void
md_transform_vec_by_symmetry(MDArena *arena, Vector *dst, const Vector *src, Symop rec, int no_translation)
{
	Transform temp;
	if (!rec.alive)
		return;
	if (rec.sym > 0 && rec.sym <= arena->mol->nsyms) {
		memmove(temp, arena->cellsyms[rec.sym - 1], sizeof(temp));
		if (no_translation)
			temp[9] = temp[10] = temp[11] = 0.0;
		TransformVec(dst, temp, src);
	} else *dst = *src;
	if (!no_translation) {
		Vector *vp = arena->mol->cell->axes;
		VecScaleInc(*dst, vp[0], rec.dx);
		VecScaleInc(*dst, vp[1], rec.dy);
		VecScaleInc(*dst, vp[2], rec.dz);
	}
}

void
md_wrap_coordinates(MDArena *arena)
{
	/*  Calculate the offset for each fragment to wrap into the unit cell (0,0,0)-(1,1,1) */
	int i, n;
	int last_n;
	Symop last_symop;
	Molecule *mol = arena->mol;
	Atom *ap;

	if (mol == NULL || mol->natoms == 0)
		return;
		
	if (arena->nfragments == 0)
		md_find_fragments(arena);

	/*  Calculate the center of mass for each fragment  */
	for (n = 0; n < arena->nfragments; n++) {
		VecZero(arena->fragment_info[n].pos);
		arena->fragment_info[n].mass = 0.0;
	}
	for (i = 0, ap = mol->atoms; i < arena->natoms_uniq; i++, ap++) {
		n = arena->fragment_indices[i];
		VecScaleInc(arena->fragment_info[n].pos, ap->r, ap->weight);
		arena->fragment_info[n].mass += ap->weight;
	}
	for (n = 0; n < arena->nfragments; n++) {
		VecScaleSelf(arena->fragment_info[n].pos, 1.0/arena->fragment_info[n].mass);
	}
	
	/*  Calculate the offset  */
	/*  A trick: atoms belonging to one fragment are almost always located in a consequent
		positions. So we cache the 'last' information to avoid duplicate calculation
		efficiently.  */
	last_n = -1;
	last_symop.alive = 0;
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap++) {
		if (ap->symop.alive)
			n = arena->fragment_indices[ap->symbase];
		else n = arena->fragment_indices[i];
		if (n == last_n && ap->symop.alive == last_symop.alive
		&& (!last_symop.alive || (ap->symop.dx == last_symop.dx && ap->symop.dy == last_symop.dy && ap->symop.dz == last_symop.dz && ap->symop.sym == last_symop.sym))) {
			/*  Belongs to the 'last' fragment  */
			ap->wrap_dx = ap[-1].wrap_dx;
			ap->wrap_dy = ap[-1].wrap_dy;
			ap->wrap_dz = ap[-1].wrap_dz;
		} else {
			/*  Calculate the offset from the position of the center of mass */
			Vector r = arena->fragment_info[n].pos;
			TransformVec(&r, arena->mol->cell->rtr, &r);
			if (ap->symop.alive && ap->symop.sym > 0 && ap->symop.sym <= arena->mol->nsyms) {
				/*  The translational components of symop are not included  */
				TransformVec(&r, arena->mol->syms[ap->symop.sym - 1], &r);
			}
			ap->wrap_dx = (arena->periodic_a ? -floor(r.x) : 0);
			ap->wrap_dy = (arena->periodic_b ? -floor(r.y) : 0);
			ap->wrap_dz = (arena->periodic_c ? -floor(r.z) : 0);
		}
		last_n = n;
		last_symop = ap->symop;
	}
}

void
md_amend_by_symmetry(MDArena *arena)
{
	int i, natoms;
	Atom *ap, *ap0;
	if (arena->mol->nsyms == 0)
		return;
	ap = ap0 = arena->mol->atoms;
	natoms = arena->mol->natoms;
	for (i = 0; i < natoms; i++, ap++) {
		Symop symop;
		if (!ap->symop.alive)
			continue;
		symop = ap->symop;
		ap0 = arena->mol->atoms + ap->symbase;
		md_transform_vec_by_symmetry(arena, &(ap->r), &(ap0->r), symop, 0);
		md_transform_vec_by_symmetry(arena, &(ap->f), &(ap0->f), symop, 1);
		md_transform_vec_by_symmetry(arena, &(ap->v), &(ap0->v), symop, 1);
	}
	
}

void
md_snapshot(MDArena *arena, int idx)
{
	Atom *ap;
	Vector *vp;
	MDSnapshot **spp;
	size_t size;
	int i, natoms;
	if (idx < 0)
		return;
	spp = (MDSnapshot **)AssignArray(&arena->snapshots, &arena->nsnapshots, sizeof(*spp), idx, NULL);
	if (spp == NULL)
		goto low_memory;
	natoms = arena->mol->natoms;
	size = sizeof(**spp) + sizeof(Vector) * (natoms - 1) * 3;
	if (*spp == NULL)
		*spp = (MDSnapshot *)malloc(size);
	else
		*spp = (MDSnapshot *)realloc(*spp, size);
	if (*spp == NULL)
		goto low_memory;
	memset(*spp, 0, size);
	(*spp)->step = arena->step;
	(*spp)->natoms = natoms;
	vp = (*spp)->rvf;
	for (i = 0, ap = arena->mol->atoms; i < natoms; i++, ap++) {
		*vp++ = ap->r;
		*vp++ = ap->v;
		*vp++ = ap->f;
	}
	return;
  low_memory:
	md_panic(arena, "Low memory while saving snapshot");
}

void
md_restore(MDArena *arena, int idx)
{
	int i, natoms, natoms1;
	Atom *ap;
	Vector *vp;
/*	MDSnapshot **spp; */

	if (idx < 0 || idx >= arena->nsnapshots)
		return;
	vp = arena->snapshots[idx]->rvf;
	natoms = arena->mol->natoms;
	natoms1 = arena->snapshots[idx]->natoms;
	if (natoms != natoms1) {
		md_warning(arena, "restore: the number of atoms in snapshot (%d) is different from the current number of atoms (%d)\n", natoms1, natoms);
	}
	for (i = 0, ap = arena->mol->atoms; i < natoms && i < natoms1; i++, ap++) {
		ap->r = *vp++;
		ap->v = *vp++;
		ap->f = *vp++;
	}
	arena->last_verlet_step = -1;  /*  The Verlet list needs update  */
}

int
md_step(MDArena *arena)
{
	md_update_velocities(arena);
	md_update_positions(arena);
	/* md_rattle_coordinate(arena); */
	md_amend_by_symmetry(arena);
	calc_force(arena);
/*	md_calc_kinetic_energy(arena); */
	md_update_velocities(arena);
	/* md_rattle_velocity(arena); */
	return 0;
}

void
md_minimize_init(MDArena *arena)
{
	int i;
	static const Vector zerov = {0, 0, 0};
	if (arena->old_forces != NULL)
		free(arena->old_forces);
	arena->old_forces = (Vector *)calloc(sizeof(Vector), arena->mol->natoms * 2);
	if (arena->old_forces == NULL)
		md_panic(arena, ERROR_out_of_memory);
	arena->old_pos = arena->old_forces + arena->mol->natoms;
	arena->f_len2 = arena->old_f_len2 = arena->max_gradient = 0.0;
	for (i = 0; i < arena->mol->natoms; i++) {
		arena->mol->atoms[i].v = arena->mol->atoms[i].f = zerov;
	}
	arena->conv_flag = 0;
}

int
md_minimize_step(MDArena *arena)
{
	Double bk, w1, w2, w3, dump;
	Double low, mid, high, low_energy, mid_energy, high_energy, lambda;
	Double low_limit, high_limit;
	Int i, j, retval, natoms_movable;
	Atom *atoms = arena->mol->atoms;
	Atom *ap;
	Int natoms = arena->mol->natoms;
	Vector r, *vp, *vdr;
	const Double phi = 0.618033988749895;  /*  The golden ratio  */

	md_amend_by_symmetry(arena);

	w1 = w2 = 0.0;
	retval = 0;
	natoms_movable = 0;
	for (i = 0, ap = atoms, vp = arena->old_forces; i < natoms; i++, ap++, vp++) {
		if (ap->fix_force < 0)
			continue;
		w1 += VecLength2(ap->f);
		w2 += VecDot(ap->f, *vp);
		natoms_movable++;
	}

	arena->f_len2 = w1;
	if (arena->old_f_len2 == 0.0) {
		/*  New direction  */
		bk = 0.0;
	} else {
		bk = (w1 - w2) / arena->old_f_len2;
		if (bk < 0.0)
			bk = 0.0;  /*  New direction  */
	}
	/*  Update the search direction  */
	arena->old_f_len2 = arena->f_len2;
	w2 = w3 = 0.0;
	dump = 1.0;
	for (i = 0, ap = atoms; i < natoms; i++, ap++) {
		if (ap->fix_force < 0)
			continue;
		w1 = VecLength2(ap->f) * dump;
		if (!isfinite(w1))  /*  Is isfinite() available in other platforms?? */
			md_panic(arena, "the gradient at atom %d exceeded limit at step %d", i+1, arena->step);
		if (w1 > 1e4)
			dump *= 1e4 / w1;
	}
	dump = sqrt(dump);
	for (i = 0, ap = atoms, vp = arena->old_forces; i < natoms; i++, ap++, vp++) {
		if (ap->fix_force < 0)
			continue;
		*vp = ap->f;
		*(vp + natoms) = ap->r;
		ap->v.x = ap->v.x * bk + ap->f.x * dump;
		ap->v.y = ap->v.y * bk + ap->f.y * dump;
		ap->v.z = ap->v.z * bk + ap->f.z * dump;
		w1 = VecLength2(ap->v);
		w2 += w1;
	/*	if (w1 > 1e4 || !isfinite(w1))
			md_panic(arena, "the gradient at atom %d exceeded limit at step %d", i+1, arena->step); */
		if (w1 > w3)
			w3 = w1;
	}
	w3 = sqrt(w3);
	arena->max_gradient = w3;
	arena->v_len2 = w2;
/*	printf("f_len2 = %g, v_len2 = %g, bk = %g, max_grad = %g\n", arena->f_len2, arena->v_len2, bk, w3); */
	if (bk == 0.0 && w3 < arena->gradient_convergence)
		return 1;  /*  Gradient is sufficiently small  */

	/*  Proceed along ap->v until the energy increases  */
	low_limit = arena->coordinate_convergence / arena->max_gradient;
	high_limit = 0.1 / arena->max_gradient;
	low = 0.0;
	low_energy = arena->total_energy;
/*	lambda = 1e-3 * arena->f_len2 / w2; */
	lambda = high_limit;
	high = lambda;
	while (1) {
		for (j = 0, ap = atoms, vp = arena->old_pos, vdr = arena->verlets_dr; j < natoms; j++, ap++, vp++, vdr++) {
			if (ap->fix_force < 0)
				continue;
			r = ap->r;
			ap->r.x = vp->x + ap->v.x * lambda;
			ap->r.y = vp->y + ap->v.y * lambda;
			ap->r.z = vp->z + ap->v.z * lambda;
			VecDec(r, ap->r);
			VecInc(*vdr, r);
		}
		calc_force(arena);
		mid = lambda;
		mid_energy = arena->total_energy;
		if (mid_energy < low_energy) {
			/*  mid is the 'sufficiently large' step to give lower total energy  */
			if (mid == high) {
				/*  Higher limit: move by this amount  */
				retval = 0;
				goto cleanup;
			}
			break;
		}
		high = mid;
		high_energy = mid_energy;
		lambda *= 0.25;
		if (lambda < low_limit) {
			/*  Cannot find point with lower energy than the starting point  */
			/*  Restore the original position  */
			for (j = 0, ap = atoms, vp = arena->old_pos; j < natoms; j++, ap++, vp++) {
				if (ap->fix_force < 0)
					continue;
				r = ap->r;
				ap->r = *vp;
				VecDec(r, ap->r);
				VecInc(*vdr, r);
			}
			calc_force(arena);
			lambda = 0.0;
			if (bk == 0.0)
				retval = 2;  /*  Atom movement is sufficiently small  */
			goto cleanup;
		}
	}
/*	printf("Line minimization [%g, %g, %g] (energy [%g, %g, %g]) ", low, mid, high, low_energy, mid_energy, high_energy); */
	/*  low_energy >= mid_energy < high_energy  */
	/*  Binary search for minimum  */
	for (i = 0; i < 10; i++) {
		if (high - mid > mid - low) {
			lambda = high - (high - mid) * phi;
			for (j = 0, ap = atoms, vp = arena->old_pos, vdr = arena->verlets_dr; j < natoms; j++, ap++, vp++, vdr++) {
				if (ap->fix_force < 0)
					continue;
				r = ap->r;
				ap->r.x = vp->x + ap->v.x * lambda;
				ap->r.y = vp->y + ap->v.y * lambda;
				ap->r.z = vp->z + ap->v.z * lambda;
				VecDec(r, ap->r);
				VecInc(*vdr, r);
			}	
			calc_force(arena);
			if (arena->total_energy < mid_energy) {
				low = mid;
				low_energy = mid_energy;
				mid = lambda;
				mid_energy = arena->total_energy;
			} else {
				high = lambda;
				high_energy = arena->total_energy;
			}
		} else {
			lambda = mid - (mid - low) * phi;
			for (j = 0, ap = atoms, vp = arena->old_pos, vdr = arena->verlets_dr; j < natoms; j++, ap++, vp++, vdr++) {
				if (ap->fix_force < 0)
					continue;
				r = ap->r;
				ap->r.x = vp->x + ap->v.x * lambda;
				ap->r.y = vp->y + ap->v.y * lambda;
				ap->r.z = vp->z + ap->v.z * lambda;
				VecDec(r, ap->r);
				VecInc(*vdr, r);
			}	
			calc_force(arena);
			if (arena->total_energy < mid_energy) {
				high = mid;
				high_energy = mid_energy;
				mid = lambda;
				mid_energy = arena->total_energy;
			} else {
				low = lambda;
				low_energy = arena->total_energy;
			}
		}
		if (bk == 0.0 && (high - low) * arena->max_gradient < arena->coordinate_convergence) {
			retval = 2;  /*  Atom movement is sufficiently small  */
			break;
		}
	}
  cleanup:
/*	printf("Final lambda = %g (%g)\n", lambda, arena->total_energy); */
	return retval;
}

void
md_cell_recalculate(MDArena *arena)
{
	Molecule *mol = arena->mol;
	XtalCell *cell = mol->cell;
	int i;
	MoleculeCalculateCellFromAxes(mol->cell, 1);
	for (i = 0; i < mol->nsyms; i++) {
		Transform temp;
		TransformMul(temp, mol->syms[i], cell->rtr);
		TransformMul(arena->cellsyms[i], cell->tr, temp);
	}
}

/*  scale_atoms = 1: symmetry-unique atoms are transformed to keep the fractional coordinate constant */
/*  scale_atoms = 2: same as scale_atoms = 1, except that the center of mass of each fragment of connected
	atoms are scaled and the atoms within one fragment are moved by the same amount */
/*  NOTE: only symmetry-unique atoms are moved, so don't forget to do md_symmetry_amend()!  */
void
md_scale_cell(MDArena *arena, const Transform tf, int scale_atoms)
{
	Transform grad, grad_inv, grad_tr, grad_inv_tr;
	int i;
	Double volume;
	Vector v;
	Atom *ap;
	XtalCell *cell = arena->mol->cell;
	memmove(grad, cell->rtr, sizeof(grad));
	VecCross(v, cell->axes[0], cell->axes[1]);
	volume = VecDot(v, cell->axes[2]);
	TransformVec(&cell->axes[0], tf, &cell->axes[0]);
	TransformVec(&cell->axes[1], tf, &cell->axes[1]);
	TransformVec(&cell->axes[2], tf, &cell->axes[2]);
	md_cell_recalculate(arena);

	/*  Deformation gradient (temp) = celltr * old_rcelltr  */
	TransformMul(grad, cell->tr, grad);
/*	grad[9] = grad[10] = grad[11] = 0.0;  */
	TransformInvert(grad_inv, grad);
	memmove(grad_tr, grad, sizeof(grad));
	MatrixTranspose(grad_tr, grad_tr);
	memmove(grad_inv_tr, grad_inv, sizeof(grad_inv));
	MatrixTranspose(grad_inv_tr, grad_inv_tr);

	if (scale_atoms == 1) {
		/*  Scale atom positions and velocities  */
		for (i = 0, ap = arena->mol->atoms; i < arena->natoms_uniq; i++, ap++) {
			if (ap->periodic_exclude)
				continue;
			TransformVec(&ap->r, grad, &ap->r);
			MatrixVec(&ap->f, grad, &ap->f);         /* No translational component */
			MatrixVec(&ap->v, grad_inv_tr, &ap->v);  /* No translational component */
		}
	} else if (scale_atoms == 2) {
		int j;
		/*  Scale atom positions and velocities by fragments  */
		for (i = 0; i < arena->nfragments; i++) {
			VecZero(arena->fragment_info[i].pos);
			arena->fragment_info[i].mass = 0.0;
		}
		/*  Get the center of mass for each fragment  */
		for (i = 0, ap = arena->mol->atoms; i < arena->natoms_uniq; i++, ap++) {
			j = arena->fragment_indices[i];
			if (j < 0 || j >= arena->nfragments)
				continue;
			VecInc(arena->fragment_info[j].pos, ap->r);
			arena->fragment_info[j].mass += ap->weight;
		}
		/*  Calculate the offset for each fragment  */
		for (j = 0; j < arena->nfragments; j++) {
			Vector *vp = &(arena->fragment_info[j].pos);
			if (arena->fragment_info[j].mass > 0.0) {
				VecScaleSelf(*vp, 1.0 / arena->fragment_info[j].mass);
				TransformVec(&v, grad, vp);
				VecSub(*vp, v, *vp);
			} else VecZero(*vp);
		}
		/*  Transform  */
		for (i = 0, ap = arena->mol->atoms; i < arena->natoms_uniq; i++, ap++) {
			if (ap->periodic_exclude)
				continue;
			MatrixVec(&ap->f, grad, &ap->f);  /*  This is not exact  */
			MatrixVec(&ap->v, grad_inv_tr, &ap->v);
			j = arena->fragment_indices[i];
			if (j < 0 || j >= arena->nfragments)
				continue;
			VecInc(ap->r, arena->fragment_info[j].pos);
		}
	}
	
	arena->last_verlet_step = -1;
}

void
md_finish(MDArena *arena)
{
	if (arena->coord_result != NULL) {
		fflush(arena->coord_result);
	}
	if (arena->vel_result != NULL) {
		fflush(arena->vel_result);
	}
	if (arena->force_result != NULL) {
		fflush(arena->force_result);
	}
	if (arena->extend_result != NULL) {
		fflush(arena->extend_result);
	}
	
	if (arena->debug_result != NULL) {
		fflush(arena->debug_result);
	}
}

int
md_copy_coordinates_from_internal(MDArena *arena)
{
	/*  Copy the internal r/v/f to xmol  */
	int i;
	Atom *ap1, *ap2;
	if (arena->mol == NULL || arena->xmol == NULL)
		return -1;  /*  Not initialized  */
	if (arena->mol->natoms != arena->xmol->natoms)
		return -2;  /*  Number of atoms does not match  */
	for (i = 0, ap1 = arena->mol->atoms, ap2 = arena->xmol->atoms; i < arena->mol->natoms; i++, ap1 = ATOM_NEXT(ap1), ap2 = ATOM_NEXT(ap2)) {
		ap2->r = ap1->r;
		ap2->v = ap1->v;
		ap2->f = ap1->f;
	}
	if (arena->mol->cell != NULL && arena->xmol->cell != NULL)
		memmove(arena->xmol->cell, arena->mol->cell, sizeof(XtalCell));
	return 0;
}

int
md_copy_coordinates_to_internal(MDArena *arena)
{
	/*  Copy the xmol coordinates to internal  */
	int i;
	Atom *ap1, *ap2;
	if (arena->mol == NULL || arena->xmol == NULL)
		return -1;  /*  Not initialized  */
	if (arena->mol->natoms != arena->xmol->natoms)
		return -2;  /*  Number of atoms does not match  */
	for (i = 0, ap1 = arena->mol->atoms, ap2 = arena->xmol->atoms; i < arena->mol->natoms; i++, ap1 = ATOM_NEXT(ap1), ap2 = ATOM_NEXT(ap2)) {
		ap1->r = ap2->r;
	}
	if (arena->mol->cell != NULL && arena->xmol->cell != NULL)
		memmove(arena->mol->cell, arena->xmol->cell, sizeof(XtalCell));
	arena->xmol->needsMDCopyCoordinates = 0;
	return 0;
}

int
md_is_running(MDArena *arena)
{
	if (arena != NULL && arena->is_running)
		return 1;
	else return 0;
}

int
md_main(MDArena *arena, int minimize)
{
//	extern int do_callback(MDArena *);
	jmp_buf env;
	const char *msg;
	int retval = 0;
	int (*md_step_func)(MDArena *);

	if (arena->is_initialized < 2 || arena->xmol->needsMDRebuild) {
		/*  Prepare MD parameters and runtime fields  */
		msg = md_prepare(arena, 0);
		if (msg != NULL) {
			snprintf(arena->errmsg, sizeof(arena->errmsg), "%s", msg);
			return 1;
		}
		arena->xmol->needsMDCopyCoordinates = 1;  /*  Coordinates will be copied below  */
	}
	
	if (arena->xmol->needsMDCopyCoordinates) {
		MoleculeLock(arena->xmol);
		retval = md_copy_coordinates_to_internal(arena);
		MoleculeUnlock(arena->xmol);
		if (retval != 0)
			return retval;
		arena->last_verlet_step = -1;  /*  The Verlet list needs update  */
	}
	
	arena->is_running = 1;
	
	if (setjmp(env) == 0) {
		arena->setjmp_buf = &env;
		if (minimize) {
			md_minimize_init(arena);
			md_step_func = md_minimize_step;
			arena->minimize_complete = 0;
		} else {
			md_step_func = md_step;
		}

		/*  Calculate initial energies and forces  */
		arena->step = arena->start_step;
		md_amend_by_symmetry(arena);
		calc_force(arena);
		md_calc_kinetic_energy(arena);
		if (arena->step == 0 || arena->start_step == arena->end_step) {
			md_output_results(arena);
			md_output_energies(arena);
		}
		if (arena->step == 0 && arena->md_callback_func != NULL) {
			retval = (*(arena->md_callback_func))(arena);
			if (retval != 0) {
				snprintf(arena->errmsg, sizeof(arena->errmsg), "MD Interrupt");
				goto cleanup;
			}
		}

		/*  Run simulation  */
		for (arena->step = arena->start_step + 1; arena->step <= arena->end_step; arena->step++) {

			/*  Molecules may be modified from the callback procedure  */
			if (arena->is_initialized < 2 || arena->xmol->needsMDRebuild) {
				msg = md_prepare(arena, 0);
				if (msg != NULL) {
					snprintf(arena->errmsg, sizeof(arena->errmsg), "%s", msg);
					retval = 1;
					goto cleanup;
				}
			}

			retval = (*md_step_func)(arena);
			md_calc_kinetic_energy(arena);
			if (arena->rescale_temp_freq > 0 && arena->step % arena->rescale_temp_freq == 0)
				md_scale_velocities(arena);
			if (arena->reinit_temp_freq > 0 && arena->step % arena->reinit_temp_freq == 0)
				md_init_velocities(arena);
			if (arena->andersen_thermo_freq > 0 && arena->step % arena->andersen_thermo_freq == 0)
				md_andersen_thermostat(arena);
			if (arena->pressure != NULL)
				pressure_control(arena);

			if (arena->coord_output_freq > 0 && arena->step % arena->coord_output_freq == 0)
				md_output_results(arena);
			if (arena->energy_output_freq > 0 && arena->step % arena->energy_output_freq == 0)
				md_output_energies(arena);

			if (retval != 0) {
				if (minimize) {
					switch (retval) {
						case 1:
							md_log(arena, "Minimize: minimization converged because the maximum gradient becomes less than the threshold.\n");
							retval = 0;
							arena->minimize_complete = 1;
							break;
						case 2:
							md_log(arena, "Minimize: minimization converged because the maximum movement of atoms becomes less than the threshold.\n");
							retval = 0;
							arena->minimize_complete = 1;
							break;
					}
				}
				goto cleanup;
			}
			if (arena->request_abort != 0) {
				snprintf(arena->errmsg, sizeof(arena->errmsg), "MD Abort by request");
				retval = 1;
				goto cleanup;
			}

			if (arena->md_callback_func != NULL && (arena->callback_freq == 0 || arena->step % arena->callback_freq == 0)) {
				retval = (*(arena->md_callback_func))(arena);
				if (retval != 0) {
					snprintf(arena->errmsg, sizeof(arena->errmsg), "MD Interrupt");
					goto cleanup;
				}
			}
			
		}
		arena->step--;

	} else {
		/*  Return from md_panic()  */
		retval = -1;  /*  Some fatal error  */
	}

cleanup:
	if (retval == 0) {
		arena->start_step = arena->step;  /*  Prepare for next run  */
	/*	MoleculeLock(arena->xmol);
		retval = md_copy_coordinates_from_internal(arena);
		MoleculeUnlock(arena->xmol); */
	}
	
	arena->setjmp_buf = NULL;

	if (arena->is_running) {
		arena->is_running = 0;
		md_finish(arena);
	}

	return retval;
}

void
md_set_default(MDArena *arena)
{
	arena->start_step = 0;
	arena->end_step = 1000;
	arena->timestep = 1;

	arena->coord_output_freq = 10;
	arena->energy_output_freq = 10;
	arena->cutoff = 9.0;
	arena->electro_cutoff = 9.0;
	arena->pairlist_distance = 10.0;
	arena->use_xplor_shift = 1;
	arena->scale14_vdw = 0.5;
	arena->scale14_elect = 0.83;
	arena->temperature = 300.0;
	arena->velocities_read = 0;
	arena->dielectric = 4.8;
/*	arena->probe_radius = 3.2; */
	arena->probe_radius = 0.0;
	arena->surface_tension = -0.005;
	arena->surface_potential_freq = 5;

	arena->velocity_limit = 100.0;
	arena->sym_tolerance = 5e-3;
	arena->gradient_convergence = 1e-6;
	arena->coordinate_convergence = 1e-8;
	arena->relocate_center = 1;
	arena->andersen_thermo_freq = 50;
	arena->andersen_thermo_coupling = 0.1;
	arena->pressure_freq = 0;
/*	arena->pressure_coupling = 0.4;
	arena->pressure_trial_width = 0.01;
	arena->pressure_control_algorithm = 0;
	arena->pressure_fluctuate_cell_origin = 0.01;
	arena->pressure_fluctuate_cell_orientation = 0.01;
	arena->pressure[0] = 1;
	arena->pressure[1] = 0;
	arena->pressure[2] = 0;
	arena->pressure[3] = 0;
	arena->pressure[4] = 1;
	arena->pressure[5] = 0;
	arena->pressure[6] = 0;
	arena->pressure[7] = 0;
	arena->pressure[8] = 1;
	arena->cell_flexibility[0] = -1;
	arena->cell_flexibility[1] = -1;
	arena->cell_flexibility[2] = -1;
	arena->cell_flexibility[3] = 0;
	arena->cell_flexibility[4] = 0;
	arena->cell_flexibility[5] = 0;
	arena->cell_flexibility[6] = 0;
	arena->cell_flexibility[7] = 0; */
/*	arena->cella.x = 1;
	arena->cella.y = 0;
	arena->cella.z = 0;
	arena->cellb.x = 0;
	arena->cellb.y = 1;
	arena->cellb.z = 0;
	arena->cellc.x = 0;
	arena->cellc.y = 0;
	arena->cellc.z = 1; */
/*	if (arena->mol != NULL && arena->mol->box != NULL)
		md_cell_recalculate(arena); */
}

MDArena *
md_arena_new(Molecule *xmol)
{
	MDArena *arena;
	arena = (MDArena *)calloc(sizeof(MDArena), 1);
	if (arena == NULL)
		return NULL;
	arena->refCount = 1;
	md_set_default(arena);
	if (xmol != NULL) {
		md_arena_set_molecule(arena, xmol);
	}
	return arena;
}

MDArena *
md_arena_set_molecule(MDArena *arena, Molecule *xmol)
{
	/*  xmol is set to mol  */
	if (arena == NULL)
		return NULL;
	if (arena->xmol != xmol) {
		if (arena->xmol != NULL) {
			if (arena->xmol->arena == arena)
				arena->xmol->arena = NULL;
			MoleculeRelease(arena->xmol);
		}
		arena->xmol = xmol;
		if (xmol != NULL) {
			MoleculeRetain(arena->xmol);
			arena->xmol->arena = arena;
		}
	}
	
	/*  Dispose the internal cache  */
	if (arena->mol != NULL) {
		MoleculeRelease(arena->mol);
		arena->mol = NULL;
	}
	
	if (xmol != NULL) {
		/*  Create an internal copy  */
		Molecule *mol = MoleculeNew();
		Atom *ap;
		int i;
		memset(mol, 0, sizeof(Molecule));
		NewArray(&mol->atoms, &mol->natoms, gSizeOfAtomRecord, xmol->natoms);
		memmove(mol->atoms, xmol->atoms, gSizeOfAtomRecord * xmol->natoms);
		/*  Note: aniso and frames are unnecessary  */
		for (i = 0, ap = mol->atoms; i < xmol->natoms; i++, ap = ATOM_NEXT(ap)) {
			ap->aniso = NULL;
			ap->frames= NULL;
			ap->nframes = 0;
		}
		NewArray(&mol->bonds, &mol->nbonds, sizeof(Int) * 2, xmol->nbonds);
		memmove(mol->bonds, xmol->bonds, sizeof(Int) * 2 * xmol->nbonds);
		NewArray(&mol->angles, &mol->nangles, sizeof(Int) * 3, xmol->nangles);
		memmove(mol->angles, xmol->angles, sizeof(Int) * 3 * xmol->nangles);
		NewArray(&mol->dihedrals, &mol->ndihedrals, sizeof(Int) * 4, xmol->ndihedrals);
		memmove(mol->dihedrals, xmol->dihedrals, sizeof(Int) * 4 * xmol->ndihedrals);
		NewArray(&mol->impropers, &mol->nimpropers, sizeof(Int) * 4, xmol->nimpropers);
		memmove(mol->impropers, xmol->impropers, sizeof(Int) * 4 * xmol->nimpropers);
		NewArray(&mol->syms, &mol->nsyms, sizeof(Transform), xmol->nsyms);
		memmove(mol->syms, xmol->syms, sizeof(Transform) * xmol->nsyms);
		if (xmol->cell != NULL) {
			mol->cell = (XtalCell *)malloc(sizeof(XtalCell));
			memmove(mol->cell, xmol->cell, sizeof(XtalCell));
		}
		if (xmol->path != NULL)
			mol->path = strdup(xmol->path);
/*		if (xmol->box != NULL) {
			mol->box = (PeriodicBox *)malloc(sizeof(PeriodicBox));
			memmove(mol->box, xmol->box, sizeof(PeriodicBox));
		} */
		mol->arena = arena;
		mol->par = xmol->par;
		if (mol->par != NULL)
			ParameterRetain(mol->par);
		arena->mol = mol;
	}
	return arena;
}

MDArena *
md_arena_retain(MDArena *arena)
{
	if (arena != NULL)
		arena->refCount++;
	return arena;
}

void
md_arena_release(MDArena *arena)
{
	int i;
	if (arena == NULL)
		return;
	if (--arena->refCount != 0)
		return;
	if (arena->mol != NULL) {
		if (arena->mol->arena == arena)
			arena->mol->arena = NULL;
		MoleculeRelease(arena->mol);
	}
	if (arena->xmol != NULL) {
		if (arena->xmol->arena == arena)
			arena->xmol->arena = NULL;
		MoleculeRelease(arena->xmol);
	}
	if (arena->par != NULL)
		ParameterRelease(arena->par);
	if (arena->log_result_name != NULL)
		free((void *)arena->log_result_name);
	if (arena->coord_result_name != NULL)
		free((void *)arena->coord_result_name);
	if (arena->vel_result_name != NULL)
		free((void *)arena->vel_result_name);
	if (arena->force_result_name != NULL)
		free((void *)arena->force_result_name);
	if (arena->extend_result_name != NULL)
		free((void *)arena->extend_result_name);
	if (arena->debug_result_name != NULL)
		free((void *)arena->debug_result_name);
	if (arena->log_result != NULL)
		fclose(arena->log_result);
	if (arena->coord_result != NULL)
		fclose(arena->coord_result);
	if (arena->vel_result != NULL)
		fclose(arena->vel_result);
	if (arena->force_result != NULL)
		fclose(arena->force_result);
	if (arena->extend_result != NULL)
		fclose(arena->extend_result);
	if (arena->debug_result != NULL)
		fclose(arena->debug_result);
	if (arena->custom_bond_pars != NULL)
		free(arena->custom_bond_pars);
	if (arena->custom_pars != NULL)
		free(arena->custom_pars);
	if (arena->bond_par_i != NULL)
		free(arena->bond_par_i);
	if (arena->angle_par_i != NULL)
		free(arena->angle_par_i);
	if (arena->dihedral_par_i != NULL)
		free(arena->dihedral_par_i);
	if (arena->improper_par_i != NULL)
		free(arena->improper_par_i);
	if (arena->vdw_par_i != NULL)
		free(arena->vdw_par_i);
	if (arena->vdw_cache != NULL)
		free(arena->vdw_cache);
	if (arena->energies != NULL)
		free(arena->energies);
	if (arena->forces != NULL)
		free(arena->forces);
	if (arena->pair_forces != NULL)
		free(arena->pair_forces);
	if (arena->cellsyms != NULL)
		free(arena->cellsyms);
	if (arena->pressure != NULL)
		pressure_release(arena->pressure);
	if (arena->fragment_indices != NULL)
		free(arena->fragment_indices);
	if (arena->fragment_info != NULL)
		free(arena->fragment_info);
#warning "TODO: Is arena->exatoms really necessary? "
	if (arena->exatoms != NULL)
		free(arena->exatoms);
	if (arena->special_positions != NULL)
		free(arena->special_positions);
	if (arena->exlist != NULL)
		free(arena->exlist);
	if (arena->exinfo != NULL)
		free(arena->exinfo);
	if (arena->verlets != NULL)
		free(arena->verlets);
	if (arena->verlets_dr != NULL)
		free(arena->verlets_dr);
	if (arena->verlet_i != NULL)
		free(arena->verlet_i);
	if (arena->snapshots != NULL) {
		for (i = 0; i < arena->nsnapshots; i++) {
			if (arena->snapshots[i] != NULL)
				free(arena->snapshots[i]);
		}
		free(arena->snapshots);
	}
	if (arena->old_forces != NULL)
		free(arena->old_forces);
	if (arena->old_pos != NULL)
		free(arena->old_pos);
	if (arena->graphite != NULL)
		graphite_release(arena->graphite);
	if (arena->ring != NULL)
		free(arena->ring);
	free(arena);
}

/*
PROGRAM
{
	call force(f)
	do loop=1,nstep
	  v = v + (dt/2) * (f/m)
	  r = r + dt*v;
	  call rattle_coodinate
	  call force(f)
	  v = v + (dt/2) * (f/m)
	  call rattle_velocity
    enddo
}
*/
