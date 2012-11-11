/*
 *  MDEwald.c
 *
 *  Created by Toshi Nagata on 2012/11/03.
 *  Copyright 2012 Toshi Nagata. All rights reserved.
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
#include "MDEwald.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>

/*  Temporary structure for Particle Mesh Ewald (Some fields are also used for Direct Ewald)  */
struct MDPME {
	Int grid_x, grid_y, grid_z;  /*  Grid size for each cell direction  */
	fftw_complex *qarray;        /*  Q array; size grid_x * grid_y * grid_z  */
	fftw_complex *qfarray;       /*  Finv(Q) (Finv is inverse FFT); size grid_x * grid_y * grid_z  */
	fftw_complex *qqarray;       /*  F(B*C*Finv(Q)) (F is FFT); size grid_x * grid_y * grid_z  */
	double       *marray;        /*  M array; size (grid_x + grid_y + grid_z) * 2 * natoms  */
	double       *barray;        /*  B array; size grid_x + grid_y + grid_z  */
	double       *bcarray;       /*  B*C array; size grid_x * grid_y * grid_z  */
	fftw_plan plan1;             /*  Inverse FFT of Q */
	fftw_plan plan2;             /*  FFT of (B*C*Finv(Q))  */
	Int       order;             /*  Order of cardinal B spline (must be even and >=4 and <= MAX_DIM_SPLINE )  */
	Double direct_limit;         /*  Maximum (k/beta) for direct Ewald summation (default = 1.29; this is where erfc(lim*PI) = 1e-8  */
	Int *kxyz;                   /*  Reciprocal lattice indices for direct Ewald summation  */
	Int nkxyz;
	Double self_correction;      /*  Self correction term, -beta/sqrt(PI)*sigma(charge**2)  */
};

static Double *s_spline_coeffs;

#define SPLINE_COEFF(_n, _m, _k) s_spline_coeffs[(_k) + MAX_DIM_SPLINE * ((_m) + MAX_DIM_SPLINE * ((_n) - 2))]

/*  Calculate the coefficients of the cardinal B-spline functions  */
static void
s_initialize_spline_coeffs(void)
{
	Int n, m, k;
	if (s_spline_coeffs != NULL && SPLINE_COEFF(2, 0, 1) == 1.0 && SPLINE_COEFF(2, 1, 1) == -1.0)
		return;  /*  Already initialized  */
	s_spline_coeffs = (Double *)calloc(sizeof(Double), (MAX_DIM_SPLINE - 1) * MAX_DIM_SPLINE * MAX_DIM_SPLINE);
	SPLINE_COEFF(2, 0, 0) = 0.0;
	SPLINE_COEFF(2, 0, 1) = 1.0;
	SPLINE_COEFF(2, 1, 0) = 2.0;
	SPLINE_COEFF(2, 1, 1) = -1.0;
	for (n = 3; n <= MAX_DIM_SPLINE; n++) {
		for (m = 0; m < n; m++) {
			Double a, b, b_prev;
			b_prev = 0.0;
			for (k = 0; k < n; k++) {
				Int sgn, j, bin;
				b = 0.0;
				bin = 1;
				sgn = 1;
				if (m > 0) {
					for (j = k; j < n; j++) {
						b += SPLINE_COEFF(n - 1, m - 1, j) * bin * sgn;
						sgn *= -1;
						bin = bin * (j + 1) / (j + 1 - k); /* Always dividable; this is a binomial coefficient (j+1 k)  */
					}
				}
				if (k > 0)
					a = SPLINE_COEFF(n - 1, m, k - 1);
				else a = 0.0;
				a = (a + n * b - b_prev) / (n - 1);
				b_prev = b;
				SPLINE_COEFF(n, m, k) = a;
			}
		}
	}
}

static Double
s_calc_spline(Int order, Double x)
{
	Int m, i;
	Double y;
	if (x <= 0 || x >= order)
		return 0.0;
	m = floor(x);
	y = SPLINE_COEFF(order, m, order - 1);
	for (i = order - 2; i >= 0; i--)
		y = y * x + SPLINE_COEFF(order, m, i);
	return y;
}

inline static void
s_complex_add(fftw_complex dst, fftw_complex src1, fftw_complex src2)
{
	dst[0] = src1[0] + src2[0];
	dst[1] = src1[1] + src2[1];
}

inline static void
s_complex_scale_add(fftw_complex dst, fftw_complex src1, fftw_complex src2, double scale)
{
	dst[0] = src1[0] + src2[0] * scale;
	dst[1] = src1[1] + src2[1] * scale;
}

inline static void
s_complex_mul(fftw_complex dst, fftw_complex src1, fftw_complex src2)
{
	double real;
	real = src1[0] * src2[0] - src1[1] * src2[1];
	dst[1] = src1[0] * src2[1] + src1[1] * src2[0];
	dst[0] = real;
}

/*  For debug  */
static void
s_output_array_as_cube(MDArena *arena, fftw_complex *ary, const char *fname)
{
	/*  For debug  */
	/*  Output array for visualize  */
	int i, j, k, n;
	double max, min, d;
	FILE *fp = fopen(fname, "w");
	Molecule *mol = arena->xmol;
	XtalCell *cp = mol->cell;
	MDPME *pme = arena->pme;

	fprintf(fp, "PME test\n");
	fprintf(fp, "Output of the array\n");
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", -(mol->natoms), cp->origin.x, cp->origin.y, cp->origin.z);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", pme->grid_x, cp->axes[0].x / pme->grid_x, cp->axes[0].y / pme->grid_x, cp->axes[0].z / pme->grid_x);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", pme->grid_y, cp->axes[1].x / pme->grid_y, cp->axes[1].y / pme->grid_y, cp->axes[1].z / pme->grid_y);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", pme->grid_z, cp->axes[2].x / pme->grid_z, cp->axes[2].y / pme->grid_z, cp->axes[2].z / pme->grid_z);
	
	/*  Atomic information  */
	for (i = 0; i < mol->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(mol->atoms, i);
		fprintf(fp, "%5d %11.6f %11.6f %11.6f %11.6f\n", ap->atomicNumber, (double)ap->charge, ap->r.x, ap->r.y, ap->r.z);
	}
	fprintf(fp, "%5d%5d\n", 1, 1);
	max = -1e10;
	min = 1e10;
	for (i = 0; i < pme->grid_x * pme->grid_y * pme->grid_z; i++) {
		d = sqrt(ary[i][0] * ary[i][0] + ary[i][1] * ary[i][1]);
		if (d > max)
			max = d;
		if (d < min)
			min = d;
	}
	max = (max > -min ? max : -min);
	for (i = n = 0; i < pme->grid_x; i++) {
		for (j = 0; j < pme->grid_y; j++) {
			for (k = 0; k < pme->grid_z; k++) {
				fprintf(fp, " %12.5e", sqrt(ary[n][0] * ary[n][0] + ary[n][1] * ary[n][1]) / max);
				n++;
				if (k == pme->grid_z - 1 || k % 6 == 5)
					fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
}

/*  Prepare kxyz for direct Ewald sum  */
void
s_prepare_direct_ewald(MDArena *arena)
{
	Vector raxes[3];
	XtalCell *cell = arena->mol->cell;
	Int pass, count, dc, n, ninit, nstep;
	raxes[0].x = cell->rtr[0];
	raxes[0].y = cell->rtr[3];
	raxes[0].z = cell->rtr[6];
	raxes[1].x = cell->rtr[1];
	raxes[1].y = cell->rtr[4];
	raxes[1].z = cell->rtr[7];
	raxes[2].x = cell->rtr[2];
	raxes[2].y = cell->rtr[5];
	raxes[2].z = cell->rtr[8];
	ninit = 3;
	nstep = 1;
	count = 0;
	for (pass = 0; pass < 2; pass++) {
		Int *ip;
		if (pass == 1)
			ip = arena->pme->kxyz;
		for (n = ninit; n <= 50; n += nstep) {
			Int nx, ny, nz;
			Vector vx, vy, vm;
			Double len;
			dc = 0;
			for (nx = -n + 1; nx < n; nx++) {
				VecScale(vx, raxes[0], nx);
				for (ny = -n + 1; ny < n; ny++) {
					VecScale(vy, raxes[1], ny);
					for (nz = -n + 1; nz < n; nz++) {
						if (n > ninit && nx > -n + nstep && nx < n - nstep && ny > -n + nstep && ny < n - nstep && nz > -n + nstep && nz < n - nstep)
							continue;
						if (nx == 0 && ny == 0 && nz == 0)
							continue;
						VecScale(vm, raxes[2], nz);
						VecInc(vm, vx);
						VecInc(vm, vy);
						len = VecLength(vm);
						if (len / arena->ewald_beta <= arena->pme->direct_limit) {
							if (pass == 0) {
								dc++;
								count++;
							} else {
								*ip++ = nx;
								*ip++ = ny;
								*ip++ = nz;
								count--;
							}
						}
					}
				}
			}
			if ((pass == 0 && dc == 0) || (pass == 1 && count <= 0))
				break;
		}
		if (pass == 0) {
			arena->pme->kxyz = (Int *)realloc(arena->pme->kxyz, sizeof(Int) * 3 * count);
			arena->pme->nkxyz = count;
		}
	}
}		

/*  For debug  */
void
calc_direct_electrostatic(MDArena *arena)
{
	Molecule *mol = arena->mol;
	XtalCell *cell = arena->mol->cell;
	Vector *forces, *f;
	Double energy;
	Double coul;
	Int i, j, n, nx, ny, nz, zero;
	Atom *ap, *apj;
	
	forces = (Vector *)calloc(sizeof(Vector), mol->natoms);
	energy = 0.0;
	coul = COULOMBIC / arena->dielectric;
	for (n = 5; n < 50; n += 5) {
		Vector vx, vy, vz;
		Double de;
		de = 0.0;
		for (nx = -n + 1; nx < n; nx++) {
			if (nx > -n + 5 && nx < n - 5)
				continue;
			VecScale(vx, cell->axes[0], nx);
			for (ny = -n + 1; ny < n; ny++) {
				if (ny > -n + 5 && ny < n - 5)
					continue;
				VecScale(vy, cell->axes[1], ny);
				for (nz = -n + 1; nz < n; nz++) {
					if (nz > -n + 5 && nz < n - 5)
						continue;
					VecScale(vz, cell->axes[2], nz);
					zero = (nx == 0 && ny == 0 && nz == 0);
					for (i = 0, ap = mol->atoms, f = forces; i < mol->natoms; i++, ap = ATOM_NEXT(ap), f++) {
						Vector ri, rj, rij;
						ri = ap->r;
						for (j = 0, apj = mol->atoms; j < mol->natoms; j++, apj = ATOM_NEXT(apj)) {
							Double w1, w2, w3, qiqj;
							if (i == j && zero)
								continue;
							qiqj = ap->charge * apj->charge;
							rj = apj->r;
							VecSub(rij, ri, rj);
							VecInc(rij, vx);
							VecInc(rij, vy);
							VecInc(rij, vz);
							w1 = 1.0 / VecLength2(rij);
							w2 = sqrt(w1);
							de += w2 * qiqj;
							w3 = w1 * w2 * qiqj;
							f->x += rij.x * w3;
							f->y += rij.y * w3;
							f->z += rij.z * w3;
						}
					}
				}
			}
		}
		de *= 0.5 * coul;
		energy += de;
		if (arena->debug_result != NULL) {
			fprintf(arena->debug_result, "Direct electrostatic energy\n");
			fprintf(arena->debug_result, "N = %d energy = %.9g %.9g\n", n, energy * INTERNAL2KCAL, de * INTERNAL2KCAL);
			for (i = 0; i < mol->natoms; i++) {
				fprintf(arena->debug_result, "  %d: %.9g %.9g %.9g\n", i, forces[i].x * coul * INTERNAL2KCAL, forces[i].y * coul * INTERNAL2KCAL, forces[i].z * coul * INTERNAL2KCAL);
			}
		}
		if (fabs(de) < 1e-7 || fabs(de / energy) < 1e-6)
			break;
	}
}

/*  Direct Ewald Sum  */
void
calc_direct_ewald_force(MDArena *arena)
{
	Molecule *mol = arena->mol;
	XtalCell *cell = arena->mol->cell;
	Vector raxes[3], *forces;
	Double energy_rec, coul_v;
	Int i, n, nx, ny, nz, *ip;
	Atom *ap;

	/*  Reciprocal cell axes  */
	raxes[0].x = cell->rtr[0];
	raxes[0].y = cell->rtr[3];
	raxes[0].z = cell->rtr[6];
	raxes[1].x = cell->rtr[1];
	raxes[1].y = cell->rtr[4];
	raxes[1].z = cell->rtr[7];
	raxes[2].x = cell->rtr[2];
	raxes[2].y = cell->rtr[5];
	raxes[2].z = cell->rtr[8];
	energy_rec = 0.0;
	coul_v = COULOMBIC / arena->dielectric / (2 * PI * TransformDeterminant(cell->tr));
	if (arena->debug_result != NULL)
		fprintf(arena->debug_result, "Direct Ewald sum\n");
	forces = &arena->forces[kPMEIndex * arena->mol->natoms];
	for (n = 0, ip = arena->pme->kxyz; n < arena->pme->nkxyz; n++, ip += 3) {
		Vector vm, f;
		Double de, w, ww, m2;
		Double s_real, s_imag, *mp;
		nx = ip[0];
		ny = ip[1];
		nz = ip[2];
		de = 0.0;
		VecScale(vm, raxes[0], nx);
		VecScaleInc(vm, raxes[1], ny);
		VecScaleInc(vm, raxes[2], nz);
		m2 = VecLength2(vm);
		s_real = s_imag = 0;
		mp = arena->pme->marray;
		for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
			Double mr = VecDot(vm, ap->r) * 2 * PI;
			mp[0] = ap->charge * cos(mr);
			mp[1] = ap->charge * sin(mr);
			s_real += mp[0];
			s_imag += mp[1];
			mp += 2;
		}
		w = PI * PI * m2 / (arena->ewald_beta * arena->ewald_beta);
		ww = exp(-w) / m2;
		de += ww * (s_real * s_real + s_imag * s_imag);
		mp = arena->pme->marray;
		for (i = 0; i < mol->natoms; i++) {
			Double w2 = 4 * PI * ww * (s_imag * mp[0] - s_real * mp[1]);
			f.x = w2 * (nx * raxes[0].x + ny * raxes[1].x + nz * raxes[2].x);
			f.y = w2 * (nx * raxes[0].y + ny * raxes[1].y + nz * raxes[2].y);
			f.z = w2 * (nx * raxes[0].z + ny * raxes[1].z + nz * raxes[2].z);
			VecInc(forces[i], f);
			mp += 2;
		}
		if (arena->debug_result != NULL && arena->debug_output_level > 0) {
			fprintf(arena->debug_result, "(%d %d %d) %.9g %.9g %.6g\n", nx, ny, nz, (ww * (s_real * s_real + s_imag * s_imag)) *coul_v * INTERNAL2KCAL, de * coul_v * INTERNAL2KCAL, sqrt(w) / PI);
		}
		de *= coul_v;
		energy_rec += de;
	}
	for (i = 0; i < mol->natoms; i++) {
		VecScaleSelf(forces[i], coul_v);
	}
	arena->energies[kPMEIndex] += energy_rec;
	if (arena->debug_result != NULL && arena->debug_output_level > 0) {
		for (i = 0; i < mol->natoms; i++) {
			fprintf(arena->debug_result, "  reciprocal force %d = {%g,%g,%g}\n", i, forces[i].x * INTERNAL2KCAL, forces[i].y * INTERNAL2KCAL, forces[i].z * INTERNAL2KCAL);
		}
	}
}

/*  Particle Mesh Ewald Sum  */
void
calc_pme_force(MDArena *arena)
{
	Molecule *mol;
	MDPME *pme;
	Atom *ap;
	Vector r;
	Int grid_x, grid_y, grid_z, order;
	Int i, j, k, kx, ky, kz, n, n0, n1, n2;
	double m0, m1, t, *mp, *bp, *bcp;
	fftw_complex *cp;
	double volume;
	double beta2;
	double dielec_r;
	Double *energy = &arena->energies[kPMEIndex];
	Vector *forces = &arena->forces[kPMEIndex * arena->mol->natoms];

	if (arena == NULL || (mol = arena->mol) == NULL || (pme = arena->pme) == NULL)
		return;
	
	forces = &arena->forces[kPMEIndex * arena->mol->natoms];
	dielec_r = COULOMBIC / arena->dielectric;
	
	/*  Fill M and Q  */
	grid_x = pme->grid_x;
	grid_y = pme->grid_y;
	grid_z = pme->grid_z;
	order = pme->order;
	cp = pme->qarray;
	beta2 = arena->ewald_beta * arena->ewald_beta;
	memset(cp, 0, sizeof(fftw_complex) * grid_x * grid_y * grid_z);
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		/*  Calculate fractional coordinate  */
		TransformVec(&r, mol->cell->rtr, &ap->r);
		r.x *= grid_x;
		r.y *= grid_y;
		r.z *= grid_z;
		mp = pme->marray + (grid_x + grid_y + grid_z) * 2 * i;
		for (k = 0; k < grid_x; k++) {
			n0 = floor((r.x - order - k) / grid_x) + 1;
			n1 = floor((r.x - k) / grid_x);
			m0 = m1 = 0.0;
			for (n = n0; n <= n1; n++) {
				t = r.x - k - n * grid_x;
				m0 += s_calc_spline(order, t);
				m1 += (s_calc_spline(order - 1, t) - s_calc_spline(order - 1, t - 1)) * grid_x;
			}
			mp[k * 2] = m0;
			mp[k * 2 + 1] = m1;
		}
		for (k = 0; k < grid_y; k++) {
			n0 = floor((r.y - order - k) / grid_y) + 1;
			n1 = floor((r.y - k) / grid_y);
			m0 = m1 = 0.0;
			for (n = n0; n <= n1; n++) {
				t = r.y - k - n * grid_y;
				m0 += s_calc_spline(order, t);
				m1 += (s_calc_spline(order - 1, t) - s_calc_spline(order - 1, t - 1)) * grid_y;
			}
			mp[(k + grid_x) * 2] = m0;
			mp[(k + grid_x) * 2 + 1] = m1;
		}
		for (k = 0; k < grid_z; k++) {
			n0 = floor((r.z - order - k) / grid_z) + 1;
			n1 = floor((r.z - k) / grid_z);
			m0 = m1 = 0.0;
			for (n = n0; n <= n1; n++) {
				t = r.z - k - n * grid_z;
				m0 += s_calc_spline(order, t);
				m1 += (s_calc_spline(order - 1, t) - s_calc_spline(order - 1, t - 1)) * grid_z;
			}
			mp[(k + grid_x + grid_y) * 2] = m0;
			mp[(k + grid_x + grid_y) * 2 + 1] = m1;
		}
		for (kx = 0; kx < grid_x; kx++) {
			for (ky = 0; ky < grid_y; ky++) {
				for (kz = 0; kz < grid_z; kz++) {
					cp[kz + grid_z * (ky + grid_y * kx)][0] += ap->charge * mp[kx * 2] * mp[(ky + grid_x) * 2] * mp[(kz + grid_x + grid_y) * 2];
				}
			}
		}
	}
	
	/*  Q is transformed by inverse FFT into qfarray */
	fftw_execute(pme->plan1);
	
	/*  B and C are built and multiply into qfarray  */
	volume = TransformDeterminant(mol->cell->tr);
	bp = pme->barray;
	bcp = pme->bcarray;
	for (i = 0; i < 3; i++) {
		int idx_ofs, grid_i;
		switch (i) {
			case 0: grid_i = grid_x; idx_ofs = 0; break;
			case 1: grid_i = grid_y; idx_ofs = grid_x; break;
			case 2: grid_i = grid_z; idx_ofs = grid_x + grid_y; break;
		}
		for (j = 0; j < grid_i; j++) {
			fftw_complex b1, b2;
			double t;
			if (j == grid_i / 2) {
				bp[idx_ofs + j] = 0.0;
				continue;
			}
			t = 2 * PI * (order - 1) * j / grid_i;
			b1[0] = cos(t);
			b1[1] = sin(t);
			b2[0] = b2[1] = 0.0;
			for (k = 0; k < order - 1; k++) {
				double s;
				t = 2 * PI * j * k / grid_i;
				s = s_calc_spline(order, k + 1);
				b2[0] += cos(t) * s;
				b2[1] -= sin(t) * s;
			}
			t = 1.0 / (b2[0] * b2[0] + b2[1] * b2[1]);
			s_complex_mul(b1, b1, b2);
			bp[idx_ofs + j] = (b1[0] * b1[0] + b1[1] * b1[1]) * t;
		}
	}
	for (kx = 0; kx < grid_x; kx++) {
		n0 = (kx <= grid_x / 2 ? kx : kx - grid_x);
		for (ky = 0; ky < grid_y; ky++) {
			n1 = (ky <= grid_y / 2 ? ky : ky - grid_y);
			for (kz = 0; kz < grid_z; kz++) {
				Transform *rtr = &(mol->cell->rtr);
				Vector mv;
				double d;
				n2 = (kz <= grid_z / 2 ? kz : kz - grid_z);
				/*  rtr[0,3,6], rtr[1,4,7], rtr[2,5,8] are the reciprocal unit cell vectors  */
				if (n0 == 0 && n1 == 0 && n2 == 0)
					d = 0.0;
				else {
					mv.x = (*rtr)[0] * n0 + (*rtr)[1] * n1 + (*rtr)[2] * n2;
					mv.y = (*rtr)[3] * n0 + (*rtr)[4] * n1 + (*rtr)[5] * n2;
					mv.z = (*rtr)[6] * n0 + (*rtr)[7] * n1 + (*rtr)[8] * n2;
					d = VecLength2(mv);
					d = exp(-PI * PI * d / beta2) / (PI * volume * d);
				}
				k = kz + grid_z * (ky + grid_y * kx);
				bcp[k] = bp[kx] * bp[grid_x + ky] * bp[grid_x + grid_y + kz] * d;
				pme->qfarray[k][0] *= bcp[k];
				pme->qfarray[k][1] *= bcp[k];
			}
		}
	}

	/*  qfarray is transformed by FFT into qqarray  */
	fftw_execute(pme->plan2);
	
	/*  Calculate energy and force  */
	*energy = 0.0;
	memset(forces, 0, sizeof(Vector) * mol->natoms);
	for (kx = 0; kx < grid_x; kx++) {
		for (ky = 0; ky < grid_y; ky++) {
			for (kz = 0; kz < grid_z; kz++) {
				k = kz + grid_z * (ky + grid_y * kx);
				*energy += pme->qarray[k][0] * pme->qqarray[k][0];
				mp = pme->marray;
				for (i = 0; i < mol->natoms; i++) {
					t = ATOM_AT_INDEX(mol->atoms, i)->charge * pme->qqarray[k][0];
					forces[i].x += t * mp[kx * 2 + 1] * mp[(ky + grid_x) * 2] * mp[(kz + grid_x + grid_y) * 2];
					forces[i].y += t * mp[kx * 2] * mp[(ky + grid_x) * 2 + 1] * mp[(kz + grid_x + grid_y) * 2];
					forces[i].z += t * mp[kx * 2] * mp[(ky + grid_x) * 2] * mp[(kz + grid_x + grid_y) * 2 + 1];
					mp += (grid_x + grid_y + grid_z) * 2;
				}
			}
		}
	}
	{
		/*  Scale the force vectors and transform into cartesian  */
		Transform tf;
		memmove(&tf, &mol->cell->rtr, sizeof(Transform));
		tf[9] = tf[10] = tf[11] = 0.0;
		for (i = 0; i < mol->natoms; i++) {
			Vector v;
			VecScale(v, forces[i], dielec_r);
			TransformVec(&forces[i], tf, &v);
		}
	}
	*energy *= dielec_r * 0.5;
	
	if (arena->debug_result) {
		fprintf(arena->debug_result, "PME energy = %.9g\n", *energy * INTERNAL2KCAL);
		if (arena->debug_output_level > 0) {
			fprintf(arena->debug_result, "PME forces:\n");
			for (i = 0; i < mol->natoms; i++)
				fprintf(arena->debug_result, "%d %.9g %.9g %.9g\n", i, forces[i].x * INTERNAL2KCAL, forces[i].y * INTERNAL2KCAL, forces[i].z * INTERNAL2KCAL);
			if (arena->debug_output_level > 1) {
				fprintf(arena->debug_result, "Grid values:\n");		
				for (i = 0; i < (grid_x + grid_y + grid_z); i++) {
					fprintf(arena->debug_result, "%d %g %g\n", i, pme->marray[i * 2], pme->marray[i * 2 + 1]);
				}
			}
		}
	}
}

void
calc_ewald_force(MDArena *arena)
{
	if (arena->use_ewald != 0) {
		if (arena->use_ewald == 1)
			calc_pme_force(arena);
		else if (arena->use_ewald == 2)
			calc_direct_ewald_force(arena);
		else return;
		arena->energies[kPMEIndex] += arena->pme->self_correction;
		if (arena->debug_result != NULL && arena->debug_output_level > 0) {
			fprintf(arena->debug_result, "%s Ewald reciprocal energy = %.9g\n", (arena->use_ewald == 1 ? "Particle Mesh" : "Direct"), (arena->energies[kPMEIndex] - arena->pme->self_correction) * INTERNAL2KCAL);
			fprintf(arena->debug_result, "Self correction energy = %.9g\n", arena->pme->self_correction * INTERNAL2KCAL);
			fprintf(arena->debug_result, "Net Ewald indirect energy = %.9g\n", arena->energies[kPMEIndex] * INTERNAL2KCAL);
		}
	}
}

/*  For debug  */
int
pme_test(MDArena *arena)
{
	s_initialize_spline_coeffs();
	calc_pme_force(arena);
	return 0;
}

void
pme_release(MDArena *arena)
{
	if (arena == NULL || arena->pme == NULL)
		return;
	if (arena->pme->plan1 != 0)
		fftw_destroy_plan(arena->pme->plan1);
	if (arena->pme->plan2 != 0)
		fftw_destroy_plan(arena->pme->plan2);
	if (arena->pme->qarray != NULL)
		fftw_free(arena->pme->qarray);
	if (arena->pme->qqarray != NULL)
		fftw_free(arena->pme->qqarray);
	if (arena->pme->qfarray != NULL)
		fftw_free(arena->pme->qfarray);
	if (arena->pme->marray != NULL)
		free(arena->pme->marray);
	if (arena->pme->barray != NULL)
		free(arena->pme->barray);
	if (arena->pme->bcarray != NULL)
		free(arena->pme->bcarray);
	if (arena->pme->kxyz)
		free(arena->pme->kxyz);
	free(arena->pme);
	arena->pme = NULL;
}

void
pme_init(MDArena *arena)
{
	MDPME *pme;
	Int xyz, x_y_z;
	if (arena == NULL || arena->mol == NULL || arena->mol->cell == NULL) {
		if (arena->pme != NULL)
			pme_release(arena);
		return;
	}
	if (arena->use_ewald == 2) {
		/*  Direct Ewald  */
		if (arena->pme != NULL)
			pme_release(arena);
		arena->pme = (MDPME *)calloc(sizeof(MDPME), 1);
		if (arena->pme == NULL)
			return;
		/*  Only marray and kxyz is allocated (for force calculation)  */
		arena->pme->marray = (double *)malloc(sizeof(double) * 2 * arena->mol->natoms);
		arena->pme->direct_limit = 1.29;
		/*  kxyz is initialized  */
		s_prepare_direct_ewald(arena);
	} else {	
		s_initialize_spline_coeffs();
		xyz = arena->ewald_grid_x * arena->ewald_grid_y * arena->ewald_grid_z;
		x_y_z = arena->ewald_grid_x + arena->ewald_grid_y + arena->ewald_grid_z;
		if (xyz == 0) {
			if (arena->pme != NULL)
				pme_release(arena);
			return;
		}
		if (arena->pme == NULL) {
			arena->pme = (MDPME *)calloc(sizeof(MDPME), 1);
			if (arena->pme == NULL)
				return;
		}
		pme = arena->pme;
		if (pme->grid_x != arena->ewald_grid_x || pme->grid_y != arena->ewald_grid_y || pme->grid_z != arena->ewald_grid_z) {
			/*  Need to rebuild internal table  */
			if (pme->grid_x * pme->grid_y * pme->grid_z < xyz) {
				if (pme->qarray != NULL)
					fftw_free(pme->qarray);
				pme->qarray = fftw_alloc_complex(xyz);
				if (pme->qfarray != NULL)
					fftw_free(pme->qfarray);
				pme->qfarray = fftw_alloc_complex(xyz);
				if (pme->qqarray != NULL)
					fftw_free(pme->qqarray);
				pme->qqarray = fftw_alloc_complex(xyz);
				pme->bcarray = (double *)realloc(pme->bcarray, sizeof(double) * xyz);
			}
			pme->marray = (double *)realloc(pme->marray, sizeof(double) * x_y_z * 2 * arena->mol->natoms);
			pme->barray = (double *)realloc(pme->barray, sizeof(double) * x_y_z);
			pme->grid_x = arena->ewald_grid_x;
			pme->grid_y = arena->ewald_grid_y;
			pme->grid_z = arena->ewald_grid_z;
			if (xyz != 0 && x_y_z != 0) {
				if (pme->plan1 != NULL)
					fftw_destroy_plan(pme->plan1);
				if (pme->plan2 != NULL)
					fftw_destroy_plan(pme->plan2);
				pme->plan1 = fftw_plan_dft_3d(pme->grid_x, pme->grid_y, pme->grid_z, pme->qarray, pme->qfarray, FFTW_BACKWARD, FFTW_ESTIMATE);
				pme->plan2 = fftw_plan_dft_3d(pme->grid_x, pme->grid_y, pme->grid_z, pme->qfarray, pme->qqarray, FFTW_FORWARD, FFTW_ESTIMATE);
			}
			pme->order = 8;
		}
	}
	/*  Calculate self-correction term  */
	{
		Double en;
		Int i;
		Atom *ap;
		en = 0.0;
		for (i = 0, ap = arena->mol->atoms; i < arena->mol->natoms; i++, ap = ATOM_NEXT(ap)) {
			en += ap->charge * ap->charge;
		}
		arena->pme->self_correction = -(COULOMBIC/arena->dielectric) * arena->ewald_beta * PI2R * en;
	}
}
