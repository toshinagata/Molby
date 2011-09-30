/*
 *  MDGraphite.c
 *  Molby
 *
 *  Created by Toshi Nagata on 07/10/07.
 *  Copyright 2007 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "MDGraphite.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*  Table entry  */
/*  The table should be created for each type of atom of different vdW parameters.
	So we need to keep the vdW parameters and the table address.  */
typedef struct MDGraphiteTable {
	Double  A, B;    /*  A = 4*pow(r0, 12)*eps, B = 4*pow(r0, 6)*eps  */
	Double  offset;  /*  Offset for cutoff truncation */
	Double  cutoff2; /*  Square of the cutoff distance  */
	Double *cache;	/*  Double[Sx*Sy*Sz*4]  */
} MDGraphiteTable;

struct MDGraphiteArena {
	Int natoms_uniq;   /*  Must be the same as arena->natoms_uniq; otherwise, all internal information is rebuilt. */

	Int needs_update;   /* If non-zero, then the internal information is rebuilt at the next force calculation  */
	Int *table_indices; /* Index to tables[] for each symmetry-unique atom  */
	Double *energies;    /* Cache energies of symmetry-unique atoms */
	Double last_energy;
	Vector *last_forces;
	Int ntables;
	MDGraphiteTable *tables;
	Vector origin;               /*  Origin of the graphite plane in cartesian coordinates */
	Vector xaxis, yaxis, zaxis;  /*  The axes of the graphite plane in cartesian coordinates */
	Mat33 rcell;                 /*  Convert cartesian to graphite internal coordinates  */

	Double R;  /* C-C bond length. Default: 1.23  */
	Double C_eps;  /* LJ epsilon for graphite C atoms. Default: 0.0860 * KCAL2INTERNAL  */
	Double C_sig;  /* LJ sigma for graphite C atoms. Default: 3.3997  */

/*  The energies and forces at selected points in the asymmetric unit is
	cached in graphite->table[]. Sx, Sy, Sz and Sz0 specifies how these points
	are selected.
	The asymmetric unit is a rectangle (0, 0)-(0.75R, sqrt(3)*0.5R). This rectangle
	area (including borders) is split into an (Sx * Sy) lattice, and the lattice 
	points (i, j) are selected as sample points. For each point, the z-coordinate
	in the range [Z_min, Z_lim] is sampled, so that there are Sz0 points span 
	[Z_min, Z_switch) and (Sz - Sz0) points span [Z_switch, Z_lim]. Thus, the
	cartesian coordinates of a sample point (i, j, k) are:
		x = 0.75*R*i/(Sx-1),
		y = sqrt(3)*0.5R*j/(Sy-1),
		z = (k >= Sz0 ? Z_switch+(Z_lim-Z_switch)*(k-Sz0)/(Sz-Sz0-1)
					  : Z_min+(Z_switch-Z_min)*k/Sz0)
	The cached value can be accessed by:
		table[((i*Sy + j)*Sz + k)*4 + n]
	where n is 0, 1, 2, 3 for the energy, force.x, force.y, force.z, respectively.
*/
	Int use_table;
	Int Sx;  /*  Default: 16  */
	Int Sy;  /*  Default: 16  */
	Int Sz;  /*  Default: 64  */
	Int Sz0; /*  Default: 52  */
	Double Z_min;     /*  Default: 1.0  */
	Double Z_switch;  /*  Default: 5.0  */
	Double Z_lim;     /*  Default: 10.0  */
	
	/*  If periodic_average > 0 and periodic boundary conditions are enabled, then
		the graphite potential will be averaged for periodic_average number of image
		atoms for each periodic cell axes. */
	Int periodic_average;

#if DEBUG
	FILE *debug_fp;
#endif

};

/*  A = 4*pow(r0, 12)*eps, B = 4*pow(r0, 6)*eps  */
static Double
s_graphite_lj(Vector r, const MDGraphiteTable *tr, Vector *force)
{
	Double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
	Double w2, w6, w12, k0, k1;
	if (r2 >= tr->cutoff2) {
		force->x = force->y = force->z = 0.0;
		return 0.0;
	}
	w2 = 1.0 / r2;
	w6 = w2 * w2 * w2;
	w12 = w6 * w6;
	k0 = tr->A * w12 - tr->B * w6 + tr->offset;
	k1 = 6 * w2 * (2.0 * tr->A * w12 - tr->B * w6);
	if (force != NULL) {
		force->x = r.x * k1;
		force->y = r.y * k1;
		force->z = r.z * k1;
	}
	return k0;
}

/*  Calculate the total potential and force between an atom and the graphite sheet.
	The graphite sheet is assumed as a perfect hexagonal array of carbon atoms, and
	only van der Waals force is considered.
	The origin is at the center of a hexagon, and all carbon atoms lie on the xy plane.
	The x-axis points from the origin to one carbon atom. The carbon atoms can be
	indexed as:
		(0.5R*(3*i+1), sqrt(3)*0.5R*j), where i+j is odd
	This function calculates the vdW potential and force from the graphite 
	atoms in the range (i,j)=(-10,-15) to (+10,+15).
	It is expected that the target atom is close to the z-axis. In the present code,
	this function is only invoked for the position in the asymmetric unit, (0,0)..
	(1.5R, sqrt(3)*0.5R).
*/
static Double
s_graphite_asymunit(MDGraphiteArena *graphite, Vector r, const MDGraphiteTable *tr, Vector *force)
{
	Int i, j;
	Double e0, en;
	Vector r0, f0;
	force->x = force->y = force->z = 0.0;
	en = 0.0;
	for (i = -10; i <= 10; i++) {
		for (j = -15; j <= 15; j++) {
			if ((i + j) % 2 == 0)
				continue;
			r0.x = r.x - 0.5 * graphite->R * (3 * i + 1);
			r0.y = r.y - graphite->R * 0.866025403784439 * j;
			r0.z = r.z;
			e0 = s_graphite_lj(r0, tr, &f0);
			VecInc(*force, f0);
			en += e0;
			r0.x = r.x - 0.5 * graphite->R * (3 * i - 1);
			e0 = s_graphite_lj(r0, tr, &f0);
			VecInc(*force, f0);
			en += e0;
		}
	}
	return en;
}

/*  Nearest graphite atoms; x, y coordinates should be multiplied by
R/2 and (sqrt(3)/2)R, respectively  */
/*
static Int sNearestTable[] = {
	-1, -3,  1, -3,
	-4, -2, -2, -2,  2, -2,  4, -2,
	-5, -1, -1, -1,  1, -1,  5, -1,
	-4,  0, -2,  0,  2,  0,  4,  0,
	-5,  1, -1,  1,  1,  1,  5,  1,
	-4,  2, -2,  2,  2,  2,  4,  2,
	-1,  3,  1,  3
}; */
static Int sNearestTable[] = {
	-2, -4,  2, -4,
	-5, -3, -1, -3,  1, -3,  5, -3,
	-4, -2, -2, -2,  2, -2,  4, -2,
	-7, -1, -5, -1, -1, -1,  1, -1,  5, -1,  7, -1,
	-8,  0, -4,  0, -2,  0,  2,  0,  4,  0,  8,  0,
	-7,  1, -5,  1, -1,  1,  1,  1,  5,  1,  7,  1,
	-4,  2, -2,  2,  2,  2,  4,  2,
	-5,  3, -1,  3,  1,  3,  5,  3,
	-2,  4,  2, -4
};
#define kSizeOfNearestTable (sizeof(sNearestTable)/sizeof(sNearestTable[0])/2)

/*  The nearest 24 graphite atoms are always explicitly calculated  */
Double
s_graphite_nearest_asymunit(MDGraphiteArena *graphite, Vector r, const MDGraphiteTable *tr, Vector *force)
{
	Int i;
	Double en;
	Vector q, dr, f;
	en = force->x = force->y = force->z = 0.0;
	for (i = 0; i < kSizeOfNearestTable; i++) {
		q.x = sNearestTable[i*2] * graphite->R * 0.5;
		q.y = sNearestTable[i*2+1] * 0.866025403784439 * graphite->R;
		dr.x = r.x - q.x;
		dr.y = r.y - q.y;
		dr.z = r.z;
		en += s_graphite_lj(dr, tr, &f);
		force->x += f.x;
		force->y += f.y;
		force->z += f.z;
	}
	return en;
}

static void
s_graphite_make_table(MDGraphiteArena *graphite, MDGraphiteTable *tr)
{
	Int i, j, k;
	Vector v, f, f0;
	Double en, en0, *fp;
	k = graphite->Sx * graphite->Sy * graphite->Sz * sizeof(Double) * 4;
	if (tr->cache == NULL)
		tr->cache = (Double *)malloc(k);
	else
		tr->cache = (Double *)realloc(tr->cache, k);
	memset(tr->cache, 0, k);
	for (k = 0; k < graphite->Sz; k++) {
		if (k >= graphite->Sz0)
			v.z = graphite->Z_switch + (graphite->Z_lim - graphite->Z_switch) * (k - graphite->Sz0) / (graphite->Sz - graphite->Sz0 - 1);
		else
			v.z = graphite->Z_min + (graphite->Z_switch - graphite->Z_min) * k / graphite->Sz0;  /*  Not (Sz0 - 1)  */
		for (i = 0; i < graphite->Sx; i++) {
			for (j = 0; j < graphite->Sy; j++) {
				v.x = 0.75 * graphite->R * i / (graphite->Sx - 1.0);
				v.y = 0.866025403784439 * graphite->R * j / (graphite->Sy - 1.0);
				en = s_graphite_asymunit(graphite, v, tr, &f);
				en0 = s_graphite_nearest_asymunit(graphite, v, tr, &f0);
				fp = tr->cache + ((i * graphite->Sy + j) * graphite->Sz + k) * 4;
				*fp++ = en - en0;
				*fp++ = f.x - f0.x;
				*fp++ = f.y - f0.y;
				*fp++ = f.z - f0.z;
			}
		}
	}
}

static inline Double
s_linear_interpolate(MDGraphiteArena *graphite, Double *base, Double dx, Double dy, Double dz)
{
	Int Sy = graphite->Sy;
	Int Sz = graphite->Sz;
	Double a0 = base[0];
	Double a1 = base[Sy * Sz * 4] - a0;
	Double a2 = base[Sz * 4] - a0;
	Double a3 = base[4] - a0;
	Double a4 = base[(Sy + 1) * Sz * 4] - a0 - a1 - a2;
	Double a5 = base[(Sz + 1) * 4] - a0 - a2 - a3;
	Double a6 = base[(Sy * Sz + 1) * 4] - a0 - a1 - a3;
	Double a7 = base[((Sy + 1) * Sz + 1) * 4] - a0 - a1 - a2 - a3 - a4 - a5 - a6;
	return a0 + a1 * dx + a2 * dy + a3 * dz + a4 * dx * dy + a5 * dy * dz + a6 * dz * dx + a7 * dx * dy * dz;
}

static Double
s_graphite(MDGraphiteArena *graphite, Vector v, MDGraphiteTable *tr, Vector *f)
{
	Int i, j, k;
	Double cella, cellb;
	Double en, dx, dy, dz, *fp;
	Int ix, iy, negx, negy, negz, symop;
	Vector v0, f0;
	
	if (v.z < 0) {
		negz = 1;
		v.z = -v.z;
	} else negz = 0;
	if (v.z > graphite->Z_lim) {
		f->x = f->y = f->z = 0.0;
		return 0.0;
	}
	
	cella = 0.75 * graphite->R;
	cellb = 0.866025403784439 * graphite->R;
	dx = v.x / (cella * 4);
	ix = floor(dx + 0.5);
	dx -= ix;
	dy = v.y / (cellb * 4);
	iy = floor(dy + 0.5);
	dy -= iy;
	
	/*  Move vector to the asymmetric unit (0,0)-(0.25,0.25) */
	if (dx < 0) {
		dx = -dx;
		negx = 1;
	} else negx = 0;
	if (dy < 0) {
		dy = -dy;
		negy = 1;
	} else negy = 0;
	if (dy > 0.25) {
		if (dx > 0.25) {
			dx = 0.5 - dx;
			dy -= 0.25;
			symop = 3;
		} else {
			dy = 0.5 - dy;
			symop = 2;
		}
	} else if (dx > 0.25) {
		dx = 0.5 - dx;
		dy = 0.25 - dy;
		symop = 1;
	} else symop = 0;
	v0.x = dx * cella * 4;
	v0.y = dy * cellb * 4;
	v0.z = v.z;
	
	if (tr->cache != NULL) {
		/* 0 <= dx <= 0.25, 0 <= dy <= 0.25, 1.0 <= v.z <= 10.0 */
		/* 0 <= i <= Sx-1, 0 <= j <= Sy-1, 0 <= k <= Sz-1 */
		dx = dx * (graphite->Sx - 1.0) / 0.25;
		dy = dy * (graphite->Sy - 1.0) / 0.25;
		if (v.z >= graphite->Z_switch)
			dz = (v.z - graphite->Z_switch) * (graphite->Sz - graphite->Sz0 - 1) / (graphite->Z_lim - graphite->Z_switch) + graphite->Sz0;
		else
			dz = (v.z - graphite->Z_min) * graphite->Sz0 / (graphite->Z_switch - graphite->Z_min);
		i = floor(dx);
		dx -= i;
		j = floor(dy);
		dy -= j;
		k = floor(dz);
		dz -= k;
		if (i >= graphite->Sx - 1) {
			i = graphite->Sx - 2;
			dx = 1.0;
		}
		if (j >= graphite->Sy - 1) {
			j = graphite->Sy - 2;
			dy = 1.0;
		}
		if (k < 0) {
			k = 0;
			dz = 0.0;
		} else if (k >= graphite->Sz - 1) {
			k = graphite->Sz - 2;
			dz = 1.0;
		}
		fp = tr->cache + ((i*graphite->Sy + j)*graphite->Sz + k)*4;
		en = s_linear_interpolate(graphite, fp, dx, dy, dz);
		f0.x = s_linear_interpolate(graphite, fp + 1, dx, dy, dz);
		f0.y = s_linear_interpolate(graphite, fp + 2, dx, dy, dz);
		f0.z = s_linear_interpolate(graphite, fp + 3, dx, dy, dz);
		en += s_graphite_nearest_asymunit(graphite, v0, tr, f);
		f->x += f0.x;
		f->y += f0.y;
		f->z += f0.z;
	} else {
		en = s_graphite_asymunit(graphite, v0, tr, f);
	}
	
	/*  Back transform to original axis system  */
	switch (symop) {
		case 3: negx = !negx; break;
		case 2: negy = !negy; break;
		case 1: negx = !negx; negy = !negy; break;
		default: break;
	}
	if (negx)
		f->x = -f->x;
	if (negy)
		f->y = -f->y;
	if (negz)
		f->z = -f->z;
	return en;
}

MDGraphiteArena *
graphite_new(void)
{
	static Vector x = {1, 0, 0}, y = {0, 1, 0}, z = {0, 0, 1};
	static Mat33 u = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	MDGraphiteArena *graphite = (MDGraphiteArena *)calloc(sizeof(MDGraphiteArena), 1);
	graphite->use_table = 1;
	graphite->Sx = 16;
	graphite->Sy = 16;
	graphite->Sz = 64;
	graphite->Sz0 = 52;
	graphite->Z_min = 1.0;
	graphite->Z_switch = 5.0;
	graphite->Z_lim = 10.0;
	graphite->R = 1.23;
	graphite->C_eps = 0.0860 * KCAL2INTERNAL;
	graphite->C_sig = 3.3997;
	graphite->xaxis = x;
	graphite->yaxis = y;
	graphite->zaxis = z;
	memmove(graphite->rcell, u, sizeof(Mat33));
	return graphite;
}

static void
s_graphite_release_tables(MDGraphiteArena *graphite)
{
	int i;
	if (graphite == NULL || graphite->ntables == 0 || graphite->tables == NULL)
		return;
	for (i = 0; i < graphite->ntables; i++) {
		if (graphite->tables[i].cache != NULL)
			free(graphite->tables[i].cache);
	}
	free(graphite->tables);
	graphite->tables = NULL;
	graphite->ntables = 0;
}

void
graphite_release(MDGraphiteArena *graphite)
{
	if (graphite != NULL) {
		if (graphite->table_indices != NULL)
			free(graphite->table_indices);
		if (graphite->tables != NULL)
			s_graphite_release_tables(graphite);
	}
	free(graphite);
}

void
graphite_set_origin(MDGraphiteArena *graphite, const Vector *vp)
{
	if (graphite == NULL)
		return;
	graphite->origin = *vp;
}

/*  (Re)initialize the MDGraphiteArena  */
static void
s_graphite_init(MDGraphiteArena *graphite, MDArena *arena)
{
	MDGraphiteTable *tp;
	Atom *ap;
	int nuniq = arena->natoms_uniq;
	int i, j;
	
	/*  Reallocate the internal storage  */
	if (graphite->table_indices == NULL)
		graphite->table_indices = (Int *)malloc(sizeof(Int) * nuniq);
	else
		graphite->table_indices = (Int *)realloc(graphite->table_indices, sizeof(Int) * nuniq);
	if (graphite->energies == NULL)
		graphite->energies = (Double *)malloc(sizeof(Double) * nuniq);
	else
		graphite->energies = (Double *)realloc(graphite->energies, sizeof(Double) * nuniq);
	if (graphite->last_forces == NULL)
		graphite->last_forces = (Vector *)malloc(sizeof(Vector) * nuniq);
	else
		graphite->last_forces = (Vector *)realloc(graphite->last_forces, sizeof(Vector) * nuniq);
	if (graphite->tables != NULL)
		s_graphite_release_tables(graphite);
	
	/*  Calculate the tables  */
	md_log(arena, "Building the graphite table...\n");
	
	for (i = 0, ap = arena->mol->atoms; i < nuniq; i++, ap++) {
		Int vdw_idx;
		VdwPar *vp;
		Double sig, eps, sigc, epsc, A, B;
		vdw_idx = arena->vdw_par_i[i];
		if (vdw_idx < 0 || vdw_idx >= arena->par->nvdwPars) {
			/*  Out of range  */
			graphite->table_indices[i] = -1;
			continue;
		}
		vp = &(arena->par->vdwPars[vdw_idx]);
		if (fabs(vp->A) < 1e-15) {
			sig = 1.0;
			eps = 0.0;
		} else {
			eps = 0.25 * vp->B * vp->B / vp->A;
			sig = pow(vp->B * 0.25 / eps, 1.0/6.0);
		}
		sigc = (graphite->C_sig + sig) * 0.5;
		epsc = sqrt(graphite->C_eps * eps);
		A =	4.0 * pow(sigc, 12.0) * epsc;
		B = 4.0 * pow(sigc, 6.0) * epsc;
		/*  Are there entries with similar A/B values?  */
		for (j = 0; j < graphite->ntables; j++) {
			tp = graphite->tables + j;
			if (fabs(tp->A - A) < 1e-10 && fabs(tp->B - B) < 1e-10)
				break;
		}
		if (j >= graphite->ntables) {
			/*  Create a new entry with a table  */
			Vector v;
			j = graphite->ntables;
			tp = (MDGraphiteTable *)AssignArray(&graphite->tables, &graphite->ntables, sizeof(MDGraphiteTable), j, NULL);
			tp->A = A;
			tp->B = B;
			/*  Calculate the offset: set 0 to offset temporarily, and call s_graphite_lj */
			tp->cutoff2 = arena->cutoff * arena->cutoff + 1.0;  /* +1.0 to avoid cutoff */
			tp->offset = 0.0;
			v.x = arena->cutoff;
			v.y = v.z = 0.0;
			tp->offset = -s_graphite_lj(v, tp, &v);
			tp->cutoff2 = arena->cutoff * arena->cutoff; /* This is the real value */
			if (graphite->use_table)
				s_graphite_make_table(graphite, tp);
			else tp->cache = NULL;
		}
		graphite->table_indices[i] = j;
	}
	graphite->natoms_uniq = arena->natoms_uniq;
	md_log(arena, "%d table(s) were created %s cached data.\n", graphite->ntables, (graphite->use_table ? "with" : "without"));	
	graphite->needs_update = 0;
	
#if DEBUG
	if (arena->debug_result && arena->debug_output_level > 1) {
		graphite->debug_fp = arena->debug_result;
	} else graphite->debug_fp = NULL;
#endif	
	
}

void
graphite_force(MDGraphiteArena *graphite, MDArena *arena, Double *energy, Vector *forces)
{
	Molecule *mol = arena->mol;
	int i;
	int ix, iy, iz, px, py, pz;
	Double pxpypz_1;
	Vector *avp, *bvp, *cvp;
	if (graphite == NULL || mol == NULL || mol->natoms == 0 || arena->natoms_uniq == 0)
		return;
	if (graphite->natoms_uniq != arena->natoms_uniq || graphite->needs_update)
		s_graphite_init(graphite, arena);
	graphite->last_energy = 0;
	memset(graphite->last_forces, 0, sizeof(Vector) * arena->natoms_uniq);
	if (graphite->periodic_average && arena->mol->cell != NULL) {
		px = (arena->periodic_a ? graphite->periodic_average : 1);
		py = (arena->periodic_b ? graphite->periodic_average : 1);
		pz = (arena->periodic_c ? graphite->periodic_average : 1);
	} else {
		px = py = pz = 1;
	}
	pxpypz_1 = 1.0 / (px * py * pz);
	if (mol->cell != NULL) {
		avp = &(mol->cell->axes[0]);
		bvp = &(mol->cell->axes[1]);
		cvp = &(mol->cell->axes[2]);
	}
	for (i = 0; i < arena->natoms_uniq; i++) {
		Vector f;
		Double en;
		Double lambda, dlambda;
		Vector v = mol->atoms[i].r;
		if (mol->atoms[i].fix_force < 0)
			continue;
		if (mol->atoms[i].occupancy == 0.0)
			continue;
		if (arena->nalchem_flags > 0) {
			char c1 = arena->alchem_flags[i];
			if (c1 == 1) {
				lambda = (1.0 - arena->alchem_lambda);
				dlambda = -arena->alchem_dlambda;
			} else if (c1 == 2) {
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
		en = f.x = f.y = f.z = 0.0;
		for (ix = 0; ix < px; ix++) {
			for (iy = 0; iy < py; iy++) {
				for (iz = 0; iz < pz; iz++) {
					Vector f0, v0;
					Double en0;
					if (mol->cell != NULL) {
						v0.x = v.x + avp->x * ix + bvp->x * iy + cvp->x * iz;
						v0.y = v.y + avp->y * ix + bvp->y * iy + cvp->y * iz;
						v0.z = v.z + avp->z * ix + bvp->z * iy + cvp->z * iz;
					} else v0 = v;
					VecDec(v0, graphite->origin);
					MatrixVec(&v0, graphite->rcell, &v0);
					en0 = s_graphite(graphite, v0, graphite->tables + graphite->table_indices[i], &f0);
					en += en0;
					VecInc(f, f0);
				}
			}
		}
		en *= pxpypz_1;
		VecScaleSelf(f, pxpypz_1 * lambda);
		graphite->energies[i] = en * lambda;
		graphite->last_energy += en * lambda;
		graphite->last_forces[i] = f;
		if (dlambda != 0.0)
			arena->alchem_energy += en * dlambda;
	}
	for (i = arena->natoms_uniq; i < mol->natoms; i++) {
		Symop symop = mol->atoms[i].symop;
		if (mol->atoms[i].fix_force < 0)
			continue;
		if (symop.alive) {
			int n = mol->atoms[i].symbase;
			if (n >= 0 && n < arena->natoms_uniq)
				graphite->last_energy += graphite->energies[n];
		}
	}
	if (energy != NULL)
		*energy += graphite->last_energy;
	if (forces != NULL) {
		for (i = 0; i < arena->natoms_uniq; i++)
			VecInc(forces[i], graphite->last_forces[i]);
	}
}

void
graphite_get_axes(MDGraphiteArena *graphite, Vector *op, Vector *xp, Vector *yp, Vector *zp, Double *rp)
{
	if (graphite != NULL) {
		if (op != NULL)
			*op = graphite->origin;
		if (xp != NULL)
			*xp = graphite->xaxis;
		if (yp != NULL)
			*yp = graphite->yaxis;
		if (zp != NULL)
			*zp = graphite->zaxis;
		if (rp != NULL)
			*rp = graphite->R;
	}
}

