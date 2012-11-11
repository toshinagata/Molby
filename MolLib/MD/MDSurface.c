/*
 *  MDSurface.c
 *
 *  Created by Toshi Nagata on 2005/06/22.
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

#include "MDSurface.h"

#include <stdlib.h>
#include <string.h>

#if DEBUG
#define SP_DEBUG 0   /* Change to non-zero for very verbose debug output  */
#else
#undef SP_DEBUG
#endif

#if SP_DEBUG
#include <stdarg.h>
#endif

#if SP_DEBUG
static void
s_sp_debug(MDArena *arena, const char *fmt,...)
{
	va_list ap;
	if (arena == NULL || arena->debug_result == NULL)
		return;
	fprintf(arena->debug_result, "SP_DEBUG: ");
	va_start(ap, fmt);
	vfprintf(arena->debug_result, fmt, ap);
	va_end(ap);
	fprintf(arena->debug_result, "\n");
	fflush(arena->debug_result);
}
#endif

static void
s_allocate_sp_arena(MDArena *arena)
{
	SPArena *sarena;
	Int natoms, i;
	sarena = (SPArena *)calloc(sizeof(SPArena), 1);
	if (sarena == NULL)
		md_panic(arena, ERROR_out_of_memory);
	arena->sp_arena = sarena;
	natoms = arena->mol->natoms;
	sarena->atom_rad = (Double *)calloc(sizeof(Double), natoms);
	if (sarena->atom_rad == NULL)
		md_panic(arena, ERROR_out_of_memory);
	sarena->atom_pot = (Double *)calloc(sizeof(Double), natoms);
	if (sarena->atom_pot == NULL)
		md_panic(arena, ERROR_out_of_memory);
	sarena->atom_area = (Double *)calloc(sizeof(Double), natoms);
	if (sarena->atom_area == NULL)
		md_panic(arena, ERROR_out_of_memory);
	sarena->atom_buried = (Byte *)calloc(sizeof(Byte), natoms);
	if (sarena->atom_buried == NULL)
		md_panic(arena, ERROR_out_of_memory);
	sarena->sphere_idx = (Int *)calloc(sizeof(Int), natoms*2);
	if (sarena->sphere_idx == NULL)
		md_panic(arena, ERROR_out_of_memory);
	for (i = 0; i < natoms; i++) {
		Int idx;
		Double r;
		sarena->atom_pot[i] = arena->surface_tension * KCAL2INTERNAL;
		idx = arena->vdw_par_i[i];
		if (idx < 0) {
			r = 1.0;  /*  Assume 1.0 angstrom if the vdW parameter is missing  */
		} else {
			VdwPar *vp = &(arena->par->vdwPars[idx]);
			r = pow(vp->A * 2.0 / vp->B, 1.0/6.0) * 0.5;
		}
		sarena->atom_rad[i] = arena->probe_radius + r;
	#if SP_DEBUG
		s_sp_debug(arena, "atom_rad[%d] = %f", i+1, sarena->atom_rad[i]);
	#endif
	}
}

/*  Register a pair of crossing spheres  */
static void
s_register_sphere_record(SPArena *sarena, Int i, Int j, const Vector *vijp, Double rhi, Double rhj, Double dij)
{
	crossing_sphere_record csr;
	Vector vw1, vw2, qij;
	Double w;

	csr.i = i;
	csr.j = j;
	csr.dij = dij;
	csr.vij = *vijp;

	/*  gj is the (signed) distance from the center of sphere i to the crossing plane */
	csr.gj = 0.5 * (rhi * rhi - rhj * rhj + csr.dij * csr.dij) / csr.dij;
	/*  The radius of the intersection circle  */
	csr.rj = sqrt(rhi * rhi - csr.gj * csr.gj);

	/*  Calculate the rotation matrix  */
	/*  (cos(t)*cos(p), cos(t)*sin(p), -sin(t)) -> (1,0,0)
	 *  (-sin(p), cos(p), 0) -> (0,1,0)
	 *  (sin(t)*cos(p), sin(t)*sin(p), cos(t)) -> (0,0,1)
	 *  where vij=dij*(sin(t)*cos(p), sin(t)*sin(p), cos(t) */
	w = 1.0 / csr.dij;
	VecScale(vw1, csr.vij, w);
	w = vw1.x * vw1.x + vw1.y * vw1.y;
	if (w > 1e-10) {
		w = sqrt(w);
		vw2.x = -vw1.y / w;
		vw2.y = vw1.x / w;
	} else {
		vw2.x = 0.0;
		vw2.y = 1.0;
	}
	vw2.z = 0.0;
	VecCross(qij, vw2, vw1);
	MatrixGeneralRotation(csr.rot_c, &qij, &vw2, &vw1);
	MatrixTranspose(csr.rot, csr.rot_c);			
	csr.int_start = csr.int_end = 0;
	csr.skip = 0;
	AssignArray(&sarena->spheres, &sarena->maxspheres, sizeof(crossing_sphere_record), sarena->nspheres++, &csr);
}

/*  Comparator for sorting crossing_sphere_record  */
static int
s_compare_sphere_records(const void *a, const void *b)
{
	const crossing_sphere_record *csa = (const crossing_sphere_record *)a;
	const crossing_sphere_record *csb = (const crossing_sphere_record *)b;
	if (csa->i > csb->i)
		return 1;
	else if (csa->i < csb->i)
		return -1;
	else if (csa->gj > csb->gj)
		return 1;
	else if (csa->gj < csb->gj)
		return -1;
	else return 0;
}

#if SP_DEBUG
/*  Check the sanity of the crossing sphere list  */
static int
s_sanity_check_sphere_records(MDArena *arena)
{
	Int natoms = arena->mol->natoms;
	SPArena *sarena = arena->sp_arena;
	crossing_sphere_record *spheres = sarena->spheres;
	Int i, j, k, start, end;
	end = 0;
	for (i = 0; i < natoms; i++) {
		start = sarena->sphere_idx[i*2];
		if (end != start)
			md_warning(arena, "Sanity check: sphere_idx[%d*2] != sphere_idx[%d*2-1]", i+1, i+1);
		end = sarena->sphere_idx[i*2+1];
		for (j = start; j < end; j++) {
			/*  Check the partner record  */
			Int jstart, jend;
			jstart = sarena->sphere_idx[spheres[j].j*2];
			jend = sarena->sphere_idx[spheres[j].j*2+1];
			for (k = jstart; k < jend; k++) {
				if (spheres[k].j == i)
					break;
			}
			if (k == jend)
				md_warning(arena, "Sanity check: (%d,%d) cross but (%d,%d) do not", i+1, j+1, j+1, i+1);
		}
	}
}
#endif

/*  List all pairs of crossing spheres  */
static void
s_list_crossing_spheres(MDArena *arena)
{
	Int i, j;
	Double rhi, rhj, w;
	Int natoms = arena->mol->natoms;
	Atom *atoms = arena->mol->atoms;
	SPArena *sarena = arena->sp_arena;
	/*  Check the buried atoms  */
	for (i = 0; i < natoms; i++) {
		rhi = sarena->atom_rad[i];
		for (j = i + 1; j < natoms; j++) {
			Vector vij;
			VecSub(vij, atoms[j].r, atoms[i].r);
			w = VecLength(vij);
			rhj = sarena->atom_rad[j];
			if (w + rhi <= rhj) {
				sarena->atom_buried[i] = 1;
				break;
			}
			if (w + rhj <= rhi)
				sarena->atom_buried[j] = 1;
		}
	}
	/*  Look for the crossing pairs of spheres  */
	for (i = 0; i < natoms; i++) {
		if (sarena->atom_buried[i]) {
			/*  Do not count the buried atoms  */
		/*	sarena->sphere_idx[i * 2] = sarena->sphere_idx[i * 2 + 1] = start; */
			continue;
		}
		rhi = sarena->atom_rad[i];
		for (j = i + 1; j < natoms; j++) {
			Vector vij;
			Double dij;
			if (j == i || sarena->atom_buried[j])
				continue;
			VecSub(vij, atoms[j].r, atoms[i].r);
			dij = VecLength(vij);
			rhj = sarena->atom_rad[j];
			if (dij >= rhi + rhj || dij <= fabs(rhi - rhj))
				continue;  /*  These two spheres do not intersect  */
			s_register_sphere_record(sarena, i, j, &vij, rhi, rhj, dij);
			VecScaleSelf(vij, -1);
			s_register_sphere_record(sarena, j, i, &vij, rhj, rhi, dij);
		}
	/*	if (end - start >= 2)
			qsort(&sarena->spheres[start], end - start, sizeof(crossing_sphere_record), s_compare_sphere_records);
		sarena->sphere_idx[i * 2] = start;
		sarena->sphere_idx[i * 2 + 1] = end; */
	}
	/*  Sort by (1) i, and (2) gj (i.e. the sphere with larger enclosing volume comes earlier) */
	qsort(sarena->spheres, sarena->nspheres, sizeof(crossing_sphere_record), s_compare_sphere_records);
	/*  Set the sphere_idx[]  */
	j = 0;
	for (i = 0; i < natoms; i++) {
		sarena->sphere_idx[i*2] = j;
		while (j < sarena->nspheres && sarena->spheres[j].i <= i)
			j++;
		sarena->sphere_idx[i*2+1] = j;
	}
#if SP_DEBUG
	s_sp_debug(arena, "Number of crossing spheres: %d", sarena->nspheres);
	for (i = 0; i < natoms; i++) {
		for (j = sarena->sphere_idx[i*2]; j < sarena->sphere_idx[i*2+1]; j++) {
			crossing_sphere_record *csp = &sarena->spheres[j];
			s_sp_debug(arena, "Atom %d: entry %d, {i=%d,j=%d,vij={%f,%f,%f},dij=%f,rj=%f,gj=%f,rot={%f,%f,%f,%f,%f,%f,%f,%f,%f},int_start=%d,int_end=%d",
				i+1, j, csp->i+1, csp->j+1, csp->vij.x, csp->vij.y, csp->vij.z, csp->dij,
				csp->rj, csp->gj,
				csp->rot[0], csp->rot[1], csp->rot[2], csp->rot[3], csp->rot[4],
				csp->rot[5], csp->rot[6], csp->rot[7], csp->rot[8],
				csp->int_start, csp->int_end);
		}
	}
	s_sanity_check_sphere_records(arena);
#endif
}

/*  Examine whether a point on the sphere i is within the crossing circle j  */
static int
s_is_inside(const crossing_sphere_record *csp, const Vector *p)
{
	Double w;
	/*  w = rad * cos(t), rad is the radius of sphere i and t is the angle between csp->vij and p */
	w = VecDot(csp->vij, *p) / csp->dij;
/*	if ((w >= 0 && csp->gj > 0 && w > csp->gj) || (w <= 0 && csp->gj < 0 && w < csp->gj)) */
	if (w >= csp->gj)
		return 1;
	else return 0;
}

/*  Register an intersection point  */
static void
s_register_intersect(MDArena *arena, Int i, const crossing_sphere_record *csj, const crossing_sphere_record *csk, Int sgn, const Vector *qp, const Vector *pp, Int point_id)
{
	intersect_record intersect;
	SPArena *sarena = arena->sp_arena;
	intersect.i = i;
	intersect.j = csj->j;
	intersect.k = csk->j;
	intersect.idx = csj - sarena->spheres;
	if (pp != NULL)
		intersect.s = atan2(pp->y, pp->x);
	else {
		Vector p = arena->mol->atoms[i].r;
		VecSub(p, *qp, p);
		MatrixVec(&p, csj->rot, &p);
		intersect.s = atan2(p.y, p.x);
	}
	intersect.flag = 0;
	intersect.point_id = point_id;
	AssignArray(&sarena->intersects, &sarena->maxintersects, sizeof(intersect_record), sarena->nintersects++, &intersect);
#if DEBUG || 1
	if (csj->i != i || csk->i != i) {
		md_warning(arena, "Internal inconsistency: point (%d,%d,%d;%d) is registered but circle indices are (%d,%d) and (%d,%d)", i+1, csj->j+1, csk->j+1, point_id, csj->i+1, csj->j+1, csk->i+1, csk->j+1);
	}
#endif
	if (point_id >= sarena->npoints) {
		intersect_point_record ipr;
		ipr.sign = sgn;
		ipr.i = i;
		ipr.idx_j = csj - sarena->spheres;
		ipr.idx_k = csk - sarena->spheres;
		ipr.v = *qp;
		AssignArray(&sarena->points, &sarena->maxpoints, sizeof(intersect_point_record), sarena->npoints++, &ipr);
	}
}

/*  Comparator for sorting intersect_records  */
static int
s_compare_intersect_records(const void *a, const void *b)
{
	const intersect_record *ia = (const intersect_record *)a;
	const intersect_record *ib = (const intersect_record *)b;
	if (ia->i < ib->i)
		return -1;
	else if (ia->i > ib->i)
		return 1;
	else if (ia->idx < ib->idx)
		return -1;
	else if (ia->idx > ib->idx)
		return 1;
	else if (ia->s < ib->s)
		return -1;
	else if (ia->s > ib->s)
		return 1;
	else return 0;
}

/*  Simpler comparator for sorting intersect_records  */
static int
s_compare_intersect_records_s_only(const void *a, const void *b)
{
	const intersect_record *ia = (const intersect_record *)a;
	const intersect_record *ib = (const intersect_record *)b;
	if (ia->s < ib->s)
		return -1;
	else if (ia->s > ib->s)
		return 1;
	else return 0;
}

/*  Calculate the intersection point of three spheres  */
static int
s_calc_intersection_points(const crossing_sphere_record *csj, const crossing_sphere_record *csk, Vector *p1, Vector *p2)
{
	Vector vk;
	Double w1, w2, w3, w4, w5;
	/*  vk is the center of circle k in (i,j)-local coordinates */
	w1 = csk->gj / csk->dij;
	MatrixVec(&vk, csj->rot, &csk->vij); /*quat_rotate(&vk, &csk->vij, &csj->quat);*/
	VecScaleSelf(vk, w1);
	w1 = vk.x * vk.x + vk.y * vk.y;
	if (w1 < 1e-10) {
		/*  Special case: vk is parallel to (0,0,1)  */
		/*  They never intersect  */
		return 0;
	}
	w2 = csk->gj * csk->gj - vk.z * csj->gj;
	w3 = csj->rj * csj->rj * w1 - w2 * w2;
	/*  The intersection point is given in (i,j)-coordinates
	 *    x = (vk.x*w2 +- vk.y*sqrt(w3))/w1
	 *    y = (vk.y*w2 -+ vk.x*sqrt(w3))/w1
	 *    z = csj->gj  */
	if (w3 <= 0.0) {
		/*  The discriminant is non-positive; no intersection points */
		return 0;
	}
	w2 /= w1;
	w3 = sqrt(w3) / w1;
	w4 = vk.x * w2;
	w5 = vk.y * w3;
	w2 *= vk.y;
	w3 *= vk.x;
	p1->x = w4 + w5;
	p1->y = w2 - w3;
	p1->z = csj->gj;
	p2->x = w4 - w5;
	p2->y = w2 + w3;
	p2->z = csj->gj;
	return 1;
}

/*  List the relevant intersection points of crossing circles  */
/*  Algorithm 2  */
static void
s_list_intersection_points(MDArena *arena)
{
	Int i, j, k;
	Int natoms = arena->mol->natoms;
	Atom *atoms = arena->mol->atoms;
	SPArena *sarena = arena->sp_arena;

	for (i = 0; i < natoms; i++) {
		Int start, end;
		Vector oij = atoms[i].r;  /*  origin of the sphere i  */
		start = sarena->sphere_idx[i * 2];
		end = sarena->sphere_idx[i * 2 + 1];
		/*  Calculate the intersection points of two circles j and k  */
		/*  Only those points with sphere (atom) number i < j < k are 
		 *  explicitly calculated, and other five points are generated by
		 *  transforming coordinates  */
		for (j = start; j < end; j++) {
			crossing_sphere_record *csj = &sarena->spheres[j];
			if (csj->skip || i >= csj->j)
				continue;
			for (k = j + 1; k < end; k++) {
				crossing_sphere_record *csk = &sarena->spheres[k];
				Vector p[2];
				Int ii;
				if (csk->skip || i >= csk->j)
					continue;
				if (s_calc_intersection_points(csj, csk, &p[0], &p[1]) == 0) {
					/*  No intersection points  */
					/*  One circle may include the other, in which case
					 *  the circle k does not need further processing  */
					p[0].x = csj->rj;
					p[0].y = 0.0;
					p[0].z = csj->gj;
					MatrixVec(&p[0], csj->rot_c, &p[0]); /*quat_rotate(&p[0], &p[0], &csj->quat_c);*/
					if (s_is_inside(csk, &p[0])) {
						csj->skip = 1;
					#if SP_DEBUG
						s_sp_debug(arena, "Circle (%d,%d) is buried within circle (%d,%d)", i+1, csj->j+1, i+1, csk->j+1);
					#endif
						break;
					} else {
						p[1].x = csk->rj;
						p[1].y = 0.0;
						p[1].z = csk->gj;
						MatrixVec(&p[1], csk->rot_c, &p[1]); /*quat_rotate(&p[1], &p[1], &csk->quat_c);*/
						if (s_is_inside(csj, &p[1])) {
							csk->skip = 1;
						#if SP_DEBUG
							s_sp_debug(arena, "Circle (%d,%d) is buried within circle (%d,%d)", i+1, csk->j+1, i+1, csj->j+1);
						#endif
						}
					}
					continue;
				}
				for (ii = 0; ii < 2; ii++) {
					Vector q;
					Int jj;
					/*  q is the global coordinates of the crossing points  */
					MatrixVec(&q, csj->rot_c, &p[ii]); /*quat_rotate(&q, &p[ii], &csj->quat_c);*/
					VecInc(q, oij);
					/*  Check whether this point is inside any crossing sphere  */
					for (jj = start; jj < end; jj++) {
						crossing_sphere_record *csjj = &sarena->spheres[jj];
						Vector qk;
						Double rhoj;
						if (j == jj || k == jj || csjj->skip)
							continue;
						VecSub(qk, q, atoms[csjj->j].r);
						rhoj = sarena->atom_rad[csjj->j];
						if (VecLength2(qk) < rhoj * rhoj) {
						#if SP_DEBUG
							VecDec(q, oij);
							s_sp_debug(arena, "Point {%f,%f,%f} is rejected because it is within the circle (%d,%d)", q.x, q.y, q.z, csjj->i+1, csjj->j+1);
							VecInc(q, oij);
						#endif
							break;
						}
					}
					if (jj >= end) {
						/*  Register this point  */
						/*  The same point is registered as 6 different indices  */
						Int point_id = sarena->npoints;
						Int i1;
						crossing_sphere_record *csj1, *csk1, *csstart, *csend;
						/*  Examine whether circles j and k cross
						 * (Mathematically they should, but the (j,k) pair is independently
						 *  calculated from the (i,j) and (i,k) pairs, so the numerical
						 *  roundups may result in inconsistency)  */
						i1 = csj->j;
						csstart = &sarena->spheres[sarena->sphere_idx[i1 * 2]];
						csend = &sarena->spheres[sarena->sphere_idx[i1 * 2 + 1]];
						for (csk1 = csstart; csk1 < csend; csk1++) {
							if (csk1->j == csk->j)
								break;
						}
						if (csk1 == csend)
							continue;  /*  Ignore this point  */
						for (csj1 = csstart; csj1 < csend; csj1++) {
							if (csj1->j == i)
								break;
						}
						if (csj1 == csend)
							md_panic(arena, "BUG ALERT! Internal inconsistency: spheres (%d,%d)=(i,j) cross but (%d,%d) do not; file %s, line %d", i+1, csj->j+1, csj->j+1, i+1, __FILE__, __LINE__);						
						/*  (i,j,k)  */
						s_register_intersect(arena, i, csj, csk, ii, &q, &p[ii], point_id);
						/*  (i,k,j)  */
						s_register_intersect(arena, i, csk, csj, ii, &q, NULL, point_id);
						/*  (j,i,k)  */
						s_register_intersect(arena, i1, csj1, csk1, ii, &q, NULL, point_id);
						/*  (j,k,i)  */
						s_register_intersect(arena, i1, csk1, csj1, ii, &q, NULL, point_id);
						/*  (k,i,j)  */
						i1 = csk->j;
						csstart = &sarena->spheres[sarena->sphere_idx[i1 * 2]];
						csend = &sarena->spheres[sarena->sphere_idx[i1 * 2 + 1]];
						for (csj1 = csstart; csj1 < csend; csj1++) {
							if (csj1->j == i)
								break;
						}
						for (csk1 = csstart; csk1 < csend; csk1++) {
							if (csk1->j == csj->j)
								break;
						}
						if (csj1 == csend)
							md_panic(arena, "BUG ALERT! Internal inconsistency: spheres (%d,%d)=(i,k) cross but (%d,%d) do not; file %s, line %d", i+1, csk->j+1, csk->j+1, i+1, __FILE__, __LINE__);
						if (csk1 == csend)
							md_panic(arena, "BUG ALERT! Internal inconsistency: spheres (%d,%d)=(j,k) cross but (%d,%d) do not; file %s, line %d", csj->j+1, csk->j+1, csk->j+1, csj->j+1, __FILE__, __LINE__);
						s_register_intersect(arena, i1, csj1, csk1, ii, &q, NULL, point_id);
						/*  (k,j,i)  */
						s_register_intersect(arena, i1, csk1, csj1, ii, &q, NULL, point_id);
					}
				} /*  end loop ii  */
			} /*  end loop k  */
		} /*  end loop j  */
	} /*  end loop i  */
	/*  Sort by (1) index i, (2) sphere index idx, (3) arc parameter s  */
	qsort(sarena->intersects, sarena->nintersects, sizeof(intersect_record), s_compare_intersect_records);
	/*  Record the index range for each crossing_sphere_record  */
	j = 0;
	for (i = 0; i < sarena->nspheres; i++) {
		crossing_sphere_record *csi = &sarena->spheres[i];
		csi->int_start = j;
		while (j < sarena->nintersects && sarena->intersects[j].idx <= i)
			j++;
		csi->int_end = j;
		if (csi->int_start == csi->int_end) {
			/*  A circle with no intersection points  */
			/*  Check if this circle is buried  */
			Vector p1;
			Int idx, idx_end;
			/*  p1 is a point on this circle  */
			p1.x = csi->rj;
			p1.y = 0.0;
			p1.z = csi->gj;
			MatrixVec(&p1, csi->rot_c, &p1); /*quat_rotate(&p1, &p1, &csi->quat_c);*/
			idx_end = sarena->sphere_idx[csi->i*2+1];
			for (idx = sarena->sphere_idx[csi->i*2]; idx < idx_end; idx++) {
				/*  Is it inside any crossing circle on sphere csi->i?  */
				if (i == idx)
					continue;
				if (s_is_inside(&sarena->spheres[idx], &p1)) {
					csi->skip = 1;
				#if SP_DEBUG
					s_sp_debug(arena, "Circle (%d,%d) is buried within multiple circles on sphere %d", csi->i+1, csi->j+1, i+1);
				#endif
					break;
				}
			}
		}
	}
#if SP_DEBUG
	s_sp_debug(arena, "Number of intersection points: %d", sarena->nintersects);
	for (i = 0; i < sarena->nspheres; i++) {
		crossing_sphere_record *csi = &sarena->spheres[i];
		for (j = csi->int_start; j < csi->int_end; j++) {
			intersect_record *irj = &sarena->intersects[j];
			Vector p = sarena->points[irj->point_id].v;
		/*	s_calc_intersection_coord(&p, irj, arena); */
		/*	p.x = csi->rj * cos(irj->s);
			p.y = csi->rj * sin(irj->s);
			p.z = csi->gj;
			quat_rotate(&p, &p, &csi->quat); */
		/*	VecInc(p, arena->mol->atoms[csi->i].r); */
			s_sp_debug(arena, "Circle %d: entry %d, {i=%d,j=%d,k=%d,s=%f,id=%d,p={%f,%f,%f},idx=%d",
				i, j, irj->i+1,irj->j+1, irj->k+1, irj->s, irj->point_id, p.x, p.y, p.z,
				irj->idx);
		}
	}
#endif
}

/*  Re-calculate the coordinates and arc parameters of intersection points  */
/*  Note: the intersection points are sorted to maintain the ascending order of
 *  the arc parameters  */
/*  If calculation failed, then returns zero. In that case, the lists of the
 *  crossing spheres and intersection points must be rebuilt  */
static int
s_recalc_intersection_points(MDArena *arena)
{
	Int i;
	Atom *atoms = arena->mol->atoms;
	SPArena *sarena = arena->sp_arena;
	Int npoints = sarena->npoints;
	intersect_point_record *ipr;
	intersect_record *ir;
	crossing_sphere_record *csr;
	Vector p1, p2;

	/*  Recalculate the intersection points  */
	for (i = 0, ipr = sarena->points; i < npoints; i++, ipr++) {
		if (s_calc_intersection_points(&sarena->spheres[ipr->idx_j], &sarena->spheres[ipr->idx_k], &p1, &p2) == 0)
			return 0;  /*  Failure  */
		if (ipr->sign == 0) {
			ipr->v = p1;
			if (i < npoints - 1 && ipr[1].i == ipr->i && ipr[1].idx_j == ipr->idx_j && ipr[1].idx_k == ipr->idx_k) {
				ipr[1].v = p2;
				i++;
				ipr++;
			}
		} else {
			ipr->v = p2;
		}
	}
	
	/*  Recalculate the arc parameters  */
	for (i = 0, ir = sarena->intersects; i < sarena->nintersects; i++, ir++) {
		ipr = &sarena->points[ir->point_id];
		csr = &sarena->spheres[ir->idx];
		p1 = atoms[ir->i].r;
		VecSub(p2, ipr->v, p1);
		MatrixVec(&p2, csr->rot, &p2);
		ir->s = atan2(p2.y, p2.x);
	}
	
	/*  Sort the points by arc parameters  */
	for (i = 0, csr = sarena->spheres; i < sarena->nspheres; i++) {
		Int start = csr->int_start;
		Int end = csr->int_end;
		if (start - end >= 2) {
			qsort(&sarena->intersects[start], end - start, sizeof(intersect_record), s_compare_intersect_records_s_only);
		}
	}
	
	return 0;
}

#if SP_DEBUG
static void
s_print_internal_information(MDArena *arena)
{
	struct temp_record {
		int i, j, k, count;
		Vector p;
	};
	int i;
	Int ntemps = 0;
	struct temp_record *temps = NULL;
	Vector p, dp;
	int k, ii, jj, kk;
	struct temp_record *tp = NULL;

	if (arena->debug_result == NULL || arena->debug_output_level == 0)
		return;
	
	fprintf(arena->debug_result, "The spheres and the crossing points:\n");
	fprintf(arena->debug_result, "%d\n", arena->mol->natoms + arena->sp_arena->nintersects / 6);
	for (i = 0; i < arena->mol->natoms; i++) {
		Vector *rp = &arena->mol->atoms[i].r;
		fprintf(arena->debug_result, "XX %f %f %f %f\n", rp->x, rp->y, rp->z, arena->sp_arena->atom_rad[i]);
	}
	for (i = 0; i < arena->sp_arena->nintersects; i++) {
		intersect_record *ip = &arena->sp_arena->intersects[i];
		crossing_sphere_record *csi = &arena->sp_arena->spheres[ip->idx];
		p.x = csi->rj * cos(ip->s);
		p.y = csi->rj * sin(ip->s);
		p.z = csi->gj;
		MatrixVec(&p, csi->rot_c, &p); /*quat_rotate(&p, &p, &csi->quat_c);*/
		VecInc(p, arena->mol->atoms[csi->i].r);
		ii = ip->i;
		jj = ip->j;
		kk = ip->k;
		if (ii > jj) { k = ii; ii = jj; jj = k; }
		if (jj > kk) { k = jj; jj = kk; kk = k; }
		if (ii > jj) { k = ii; ii = jj; jj = k; }
		/*  Search for temps[], and register this point if not present  */
		for (k = 0; k < ntemps; k++) {
			tp = &temps[k];
			VecSub(dp, p, tp->p);
			if (VecLength2(dp) < 1e-6)
				break;
		}
		if (k < ntemps) {
			if (tp->i != ii || tp->j != jj || tp->k != kk) {
				fprintf(arena->debug_result, "!!! The point {%f,%f,%f} appeared twice with difference indices {%d,%d,%d} and {%d,%d,%d}\n",
					p.x, p.y, p.z, tp->i+1, tp->j+1, tp->k+1, ii+1, jj+1, kk+1);
			} else {
				tp->count++;
			}
		} else {
			tp = (struct temp_record *)AssignArray(&temps, &ntemps, sizeof(temps[0]), ntemps, NULL);
			tp->i = ii;
			tp->j = jj;
			tp->k = kk;
			tp->count = 1;
			tp->p = p;
		}
	}
	for (i = 0; i < ntemps; i++) {
		p = temps[i].p;
		fprintf(arena->debug_result, "XX %f %f %f %f", p.x, p.y, p.z, 0.500);
		if (temps[i].count != 6) {
			fprintf(arena->debug_result, "  !!! count = %d != 6 (i=%d,j=%d,k=%d)", temps[i].count, temps[i].i+1, temps[i].j+1, temps[i].k+1);
		}
		fprintf(arena->debug_result, "\n");
	}
	fflush(arena->debug_result);
	if (temps != NULL)
		free(temps);
}
#endif

/*  The gradient vectors for force calculation  */
/*  For force calculation we need three gradients, namely for vi, vj, vk 
 *  (these are the center of spheres i, j, k in global coordinates).
 *  Actually we calculate two gradients for xj and xk, where xj = vj - vi
 *  and xk = vk - vi. Gradients for vi, vj, vk can then be calculated 
 *  by grad_vi = -(grad_xj + grad_xk), grad_vj = grad_xj, grad_vk = grad_xk.
 */
/*  Notations:
 *  q : the position of the intersection point in the (i,j) local coordinates
 *  (i.e. the coordinates in which the center of sphere i is (0,0,0) and
 *  the center of sphere j is (0,0,csr1->gj) )
 *  gvj : (0,0,csr1->gj) = the center of circle j in the (i,j) local coordinates
 *  gvk : the center of circle k in the (i,j) local coordinates
 *  gvj and gvk are related to xj and xk by the following equations:
 *  gvj = quat_rotate((csr1->gj / csr1->dij) * csr1->vij, csr1->quat)
 *  gvk = quat_rotate((csr2->gj / csr2->dij) * csr2->vij, csr1->quat)
 *  where csr1 and csr2 are the crossing_sphere_records for (i,j) and (i,k)
 *  pairs, respectively. Note that the quaternion in the gvk equation is
 *  csr1->quat, not csr2->quat!
 */
/*  The following gradients are calculated here:
 *  qx_xj = grad_xj(qx) = (d(q.x)/d(xj.x), d(q.x)/d(xj.y), d(q.x)/d(xj.z))
 *  qy_xj = grad_xj(qy) = (d(q.y)/d(xj.x), d(q.y)/d(xj.y), d(q.y)/d(xj.z))
 *  gj_xj = grad_xj(gj) = (d(gj)/d(xj.x), d(gj)/d(xj.y), d(gj)/d(xj.z))
 *  s_xj = grad_xj(s) = (ds/d(xj.x), ds/d(xj.y), ds/d(xj.z))
 *  qx_xk = grad_xk(qx) = (d(q.x)/d(xk.x), d(q.x)/d(xk.y), d(q.x)/d(xk.z))
 *  qy_xk = grad_xk(qy) = (d(q.y)/d(xk.x), d(q.y)/d(xk.y), d(q.y)/d(xk.z))
 *  gk_xk = grad_xk(gk) = (d(gk)/d(xk.x), d(gk)/d(xk.y), d(gk)/d(xk.z))
 *  s_xk = grad_xk(s) = (ds/d(xk.x), ds/d(xk.y), ds/d(xk.z))
 *  (s is the arc parameter on circle j.)
 *  Note that grad_xk(gj) and grad_xj(gk) are zero.
 */
static void
s_calc_gradient_vectors(gradient_record *gp, const MDArena *arena, const crossing_sphere_record *csr1, const crossing_sphere_record *csr2, const Vector *qp)
{
	Vector q, gvj, gvk;
	Vector gvkx_xj, gvky_xj, gvkz_xj;
	Vector gvkx_xk, gvky_xk, gvkz_xk;
	Vector va1, va2;
	Double rho_i, rho_j, rho_k;
	Double w1, w2, w3 ,w4, w5, w6;
	const Double *atom_rad;

	/*  The crossing point in (i,j) local coordinates  */
	MatrixVec(&q, csr1->rot, qp); /*quat_rotate(&q, qp, &csr1->quat); *//* *qp is the crossing point in global coords */
	
	/*  The radii of the spheres  */
	atom_rad = arena->sp_arena->atom_rad;
	rho_i = atom_rad[csr1->i];
	rho_j = atom_rad[csr1->j];
	rho_k = atom_rad[csr2->j];
	
	/*  The center of sphere k in (i,j) local coordinates  */
	gvj.x = gvj.y = 0.0;
	gvj.z = csr1->gj;
	w1 = csr2->gj / csr2->dij;
	VecScale(gvk, csr2->vij, w1);
	MatrixVec(&gvk, csr1->rot, &gvk); /*quat_rotate(&gvk, &gvk, &csr1->quat);*/

	/*  Gradients by xj  */

	/*  Gradients of local-to-global matrix M  */
	/*  M converts the (i,j)-local coordinates to the global coordinates as:
	 *   x(global).x = M[0]*x(local).x + M[1]*x(local).y + M[2]*x(local).z, etc.
	 *  mx, my, mz are defined as follows:
	 *   mx[n] = d(M[n])/d(xj.x)
	 *   my[n] = d(M[n])/d(xj.y)
	 *   mz[n] = d(M[n])/d(xj.z)
	*/
	{
		Double wxx, wyy, wzz, wdd, w2sq;
		wxx = csr1->vij.x * csr1->vij.x;
		wyy = csr1->vij.y * csr1->vij.y;
		wzz = csr1->vij.z * csr1->vij.z;
		wdd = csr1->dij * csr1->dij;
		w2 = wxx + wyy;
		w2sq = sqrt(w2);
		w3 = 1.0 / (wdd * csr1->dij);
		w4 = 1.0 / (w2 * w2sq);
		w5 = w3 * w4;
		gp->mx[0] = csr1->vij.z * w5 * (-wxx * w2 + wyy * wdd);
		gp->mx[1] = -csr1->vij.x * csr1->vij.y * csr1->vij.z * w5 * (w2 + wdd);
		w6 = csr1->vij.z * csr1->vij.z * w2 * w5;
		gp->mx[2] = -csr1->vij.x * w6;
		gp->my[0] = gp->mx[1];
		gp->my[1] = csr1->vij.z * w5 * (-wyy * w2 + wxx * wdd);
		gp->my[2] = -csr1->vij.y * w6;
		w6 = w2sq * w3;
		gp->mz[0] = csr1->vij.x * w6;
		gp->mz[1] = csr1->vij.y * w6;
		gp->mz[2] = csr1->vij.z * w6;
		gp->mx[3] = csr1->vij.x * csr1->vij.y * w4;
		gp->mx[4] = wyy * w4;
		gp->my[3] = -wxx * w4;
		gp->my[4] = -gp->mx[3];
		gp->mx[5] = gp->my[5] = 0.0;
		gp->mz[3] = gp->mz[4] = gp->mz[5] = 0.0;
		gp->mx[6] = (wyy + wzz) * w3;
		gp->mx[7] = -csr1->vij.x * csr1->vij.y * w3;
		gp->mx[8] = -csr1->vij.x * csr1->vij.z * w3;
		gp->my[6] = gp->mx[7];
		gp->my[7] = (wxx + wzz) * w3;
		gp->my[8] = -csr1->vij.y * csr1->vij.z * w3;
		gp->mz[6] = gp->mx[8];
		gp->mz[7] = gp->my[8];
		gp->mz[8] = (wxx + wyy) * w3;
	}
	
	/*  gvkx_xj = grad_xj(gvk.x) = gk/vk*grad_xj((M^-1)*xk)
	 *   = gk/vk*grad_xj((M[0],M[3],M[6])*xk)
	 *   = gk/vk*{(mx[0],mx[3],mx[6])*xk, (my[0],my[3],my[6])*xk, (mz[0],mz[3],mz[6])*xk}
	 *  Also for gvky_xj, gvkz_xj  */
	va1 = csr2->vij;
	w1 = csr2->gj / csr2->dij;
	gvkx_xj.x = w1 * (gp->mx[0] * va1.x + gp->mx[1] * va1.y + gp->mx[2] * va1.z);
	gvkx_xj.y = w1 * (gp->my[0] * va1.x + gp->my[1] * va1.y + gp->my[2] * va1.z);
	gvkx_xj.z = w1 * (gp->mz[0] * va1.x + gp->mz[1] * va1.y + gp->mz[2] * va1.z);
	gvky_xj.x = w1 * (gp->mx[3] * va1.x + gp->mx[4] * va1.y + gp->mx[5] * va1.z);
	gvky_xj.y = w1 * (gp->my[3] * va1.x + gp->my[4] * va1.y);
	gvky_xj.z = 0.0;
	gvkz_xj.x = w1 * (gp->mx[6] * va1.x + gp->mx[7] * va1.y + gp->mx[8] * va1.z);
	gvkz_xj.y = w1 * (gp->my[6] * va1.x + gp->my[7] * va1.y + gp->my[8] * va1.z);
	gvkz_xj.z = w1 * (gp->mz[6] * va1.x + gp->mz[7] * va1.y + gp->mz[8] * va1.z);

#if 0
	/*  Intermediate vectors: gradients of gvk.x, gvk.y, gvk.y by xj  */
	/*  gvkx_xj = grad_xj(gvk.x)  */
	w2 = csr1->vij.x * csr1->vij.x + csr1->vij.y * csr1->vij.y;
	w2sq = sqrt(w2);
	w3 = w1 / (csr1->dij * csr1->dij * csr1->dij * w2 * w2sq);
	w4 = csr1->vij.z * (-csr1->vij.x * csr1->vij.x * w2 + csr1->vij.y * csr1->vij.y * csr1->dij * csr1->dij);
	w5 = -csr1->vij.x * csr1->vij.y * csr1->vij.z * (w1 + csr1->dij * csr1->dij);
	w6 = -csr1->vij.z * csr1->vij.z * w2;
	gvkx_xj.x = w3 * (csr2->vij.x * w4 + csr2->vij.y * w5 + csr2->vij.z * csr1->vij.x * w6);
	w4 = csr1->vij.z * (-csr1->vij.y * csr1->vij.y * w2 + csr1->vij.x * csr1->vij.x * csr1->dij * csr1->dij);
	gvkx_xj.y = w3 * (csr2->vij.x * w5 + csr2->vij.y * w4 + csr2->vij.z * csr1->vij.y * w6);
	gvkx_xj.z = w3 * w2 * w2 * VecDot(csr2->vij, csr1->vij);
	/*  gvky_xj = grad_xj(gvk.y)  */
	w3 = w1 / (w2 * w2sq) * (csr2->vij.x * csr1->vij.x + csr2->vij.y * csr1->vij.y);
	gvky_xj.x = csr1->vij.y * w3;
	gvky_xj.y = -csr1->vij.x * w3;
	gvky_xj.z = 0.0;
	/*  gvkz_xj = grad_xj(gvk.z)  */
	w3 = w1 / (csr1->dij * csr1->dij * csr1->dij);
	gvkz_xj.x = w3 * csr2->vij.x * (csr1->dij * csr1->dij - csr1->vij.x * csr1->vij.x);
	gvkz_xj.y = w3 * csr2->vij.y * (csr1->dij * csr1->dij - csr1->vij.y * csr1->vij.y);
	gvkz_xj.z = w3 * csr2->vij.z * (csr1->dij * csr1->dij - csr1->vij.z * csr1->vij.z);
#endif

	/*  Intermediate vectors:
	 *   va1 = gvk.x * qx_xj + gvk.y * qy_xj
	 *   va2 = q.x * qx_xj + qy * qy_xj  */
	w1 = csr1->dij;
	w1 = (rho_j * rho_j - rho_i * rho_i + w1 * w1) / (2.0 * w1 * w1 * w1);
	w2 = -gvk.z * w1;
	VecScale(va1, csr1->vij, w2);
	VecScaleSelf(gvkx_xj, q.x);
	VecScaleSelf(gvky_xj, q.y);
	VecScaleSelf(gvkz_xj, csr1->gj);
	va1.x -= gvkx_xj.x + gvky_xj.x + gvkz_xj.x;
	va1.y -= gvkx_xj.y + gvky_xj.y + gvkz_xj.y;
	va1.z -= gvkx_xj.z + gvky_xj.z + gvkz_xj.z;
	w2 = -csr1->gj * w1;
	VecScale(va2, csr1->vij, w2);
	w2 = 1.0 / (gvk.x * q.y - gvk.y * q.x);
	w3 = q.y * w2;
	w4 = -gvk.y * w2;
	w5 = -q.x * w2;
	w6 = gvk.x * w2;
	VecScale(gp->qx_xj, va1, w3);
	VecScaleInc(gp->qx_xj, va2, w4);
	VecScale(gp->qy_xj, va1, w5);
	VecScaleInc(gp->qy_xj, va2, w6);
	VecScale(gp->gj_xj, csr1->vij, w1);

/*	w2 = 1.0 / (q.x * gvk.y - q.y * gvk.x);
	w3 = w1 * w2;
	w4 = (-gvk.y * csr1->gj + gvk.z * q.y) * w3;
	VecScale(gp->qx_xj, csr1->vij, w4);
	w4 = (-gvk.x * csr1->gj - gvk.z * q.x) * w3;
	VecScale(gp->qy_xj, csr1->vij, w4);
	VecScale(gp->gj_xj, csr1->vij, w1); */
/*	gp->qx_xj.x = gp->qx_xj.y = 0.0;
	gp->qx_xj.z = (-gvk.y * csr1->gj + gvk.z * q.y) * w3;
	gp->qy_xj.x = gp->qy_xj.y = 0.0;
	gp->qy_xj.z = (gvk.x * csr1->gj - gvk.z * q.x) * w3;
	gp->gj_xj.x = gp->gj_xj.y = 0.0;
	gp->gj_xj.z = w1;
	quat_rotate(&gp->qx_xj, &gp->qx_xj, &csr1->quat);
	quat_rotate(&gp->qy_xj, &gp->qy_xj, &csr1->quat);
	quat_rotate(&gp->gj_xj, &gp->gj_xj, &csr1->quat); */

	/*  Gradients by xk  */

	/*  Intermediate vectors: gradients of gvk.x, gvk.y, gvk.y by xk  */
	/*  gvkx_xk = grad_xk(gvk.x)  */
	w1 = csr2->dij * csr2->dij;
	w2 = (rho_i * rho_i - rho_k * rho_k + w1) / (2.0 * w1);
	w3 = -(rho_i * rho_i - rho_k * rho_k) / (w1 * w1);
	va1.x = 1.0;
	va1.y = va1.z = 0.0;
	MatrixVec(&va1, csr1->rot_c, &va1); /*quat_rotate(&va1, &va1, &csr1->quat_c);*/
	w4 = VecDot(va1, csr2->vij) * w3;
	VecScale(gvkx_xk, va1, w2);
	VecScaleInc(gvkx_xk, csr2->vij, w4);
	/*  gvky_xk = grad_xk(gvk.y)  */
	va1.y = 1.0;
	va1.x = va1.z = 0.0;
	MatrixVec(&va1, csr1->rot_c, &va1); /*quat_rotate(&va1, &va1, &csr1->quat_c);*/
	w4 = VecDot(va1, csr2->vij) * w3;
	VecScale(gvky_xk, va1, w2);
	VecScaleInc(gvky_xk, csr2->vij, w4);
	/*  gvkz_xk = grad_xk(gvk.z)  */
	va1.z = 1.0;
	va1.x = va1.y = 0.0;
	MatrixVec(&va1, csr1->rot_c, &va1); /*quat_rotate(&va1, &va1, &csr1->quat_c);*/
	w4 = VecDot(va1, csr2->vij) * w3;
	VecScale(gvkz_xk, va1, w2);
	VecScaleInc(gvkz_xk, csr2->vij, w4);

	/*   va1 = gvk.x * qx_xk + gvk.y * qy_xk  */
	w1 = csr2->dij * csr2->dij;
	w1 = (rho_k * rho_k - rho_i * rho_i + w1) / (2.0 * w1 * csr2->dij);
	w2 = 2.0 * csr2->gj * w1;
	VecScale(va1, csr2->vij, w2);
	VecScaleSelf(gvkx_xk, q.x);
	VecScaleSelf(gvky_xk, q.y);
	VecScaleSelf(gvkz_xk, csr1->gj);
	va1.x -= gvkx_xk.x + gvky_xk.x + gvkz_xk.x;
	va1.y -= gvkx_xk.y + gvky_xk.y + gvkz_xk.y;
	va1.z -= gvkx_xk.z + gvky_xk.z + gvkz_xk.z;
	
	w2 = 1.0 / (gvk.x * q.y - gvk.y * q.x);
	w3 = q.y * w2;
	w5 = -q.x * w2;
	VecScale(gp->qx_xk, va1, w3);
	VecScale(gp->qy_xk, va1, w5);
	VecScale(gp->gk_xk, csr2->vij, w1);

#if 0
	w1 = csr2->dij;
	w1 = (rho_k * rho_k - rho_i * rho_i + w1 * w1) / (2.0 * w1 * w1 * w1);
	w2 = 1.0 / (q.x * gvk.y - q.y * gvk.x);
	w3 = w1 * w2;
	gp->qx_xk.x = (2.0 * gvk.x - q.x) * w3;
	gp->qx_xk.y = (2.0 * gvk.y - q.y) * w3;
	gp->qx_xk.z = (2.0 * gvk.z - csr1->gj) * w3;
	gp->qy_xk = gp->qx_xk;
	VecScaleSelf(gp->qx_xk, q.y);
	VecScaleSelf(gp->qy_xk, q.x);
/*	w3 = w1 / csr2->dij;  */
	VecScale(gp->gk_xk, csr2->vij, w1);
	MatrixVec(&gp->qx_xk, csr1->rot, &gp->qx_xk); /*quat_rotate(&gp->qx_xk, &gp->qx_xk, &csr1->quat);*/
	MatrixVec(&gp->qy_xk, csr1->rot, &gp->qy_xk); /*quat_rotate(&gp->qy_xk, &gp->qy_xk, &csr1->quat);*/
/*	quat_rotate(&gp->gk_xk, &gp->gk_xk, &csr1->quat); */
#endif

	/*  Gradients of the arc parameter s  */
	w1 = 1.0 / (csr1->rj * csr1->rj);
	gp->s_xj.x = (-q.y * gp->qx_xj.x + q.x * gp->qy_xj.x) * w1;
	gp->s_xj.y = (-q.y * gp->qx_xj.y + q.x * gp->qy_xj.y) * w1;
	gp->s_xj.z = (-q.y * gp->qx_xj.z + q.x * gp->qy_xj.z) * w1;
	gp->s_xk.x = (-q.y * gp->qx_xk.x + q.x * gp->qy_xk.x) * w1;
	gp->s_xk.y = (-q.y * gp->qx_xk.y + q.x * gp->qy_xk.y) * w1;
	gp->s_xk.z = (-q.y * gp->qx_xk.z + q.x * gp->qy_xk.z) * w1;

#if SP_DEBUG
	if (arena->debug_result != NULL && arena->debug_output_level > 2) {
		Int i, j, k, c;
		FILE *fp = arena->debug_result;  /*  Just an alias for shorter typing */
		i = csr1->i + 1;
		j = csr1->j + 1;
		k = csr2->j + 1;
		c = (gp->sign ? '+' : '-');
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: qx    =  %10.6f\n", i, j, k, c, q.x);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: qx_xj = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->qx_xj.x, gp->qx_xj.y, gp->qx_xj.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: qx_xk = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->qx_xk.x, gp->qx_xk.y, gp->qx_xk.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: qy    =  %10.6f\n", i, j, k, c, q.y);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: qy_xj = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->qy_xj.x, gp->qy_xj.y, gp->qy_xj.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: qy_xk = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->qy_xk.x, gp->qy_xk.y, gp->qy_xk.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gj    =  %10.6f\n", i, j, k, c, csr1->gj);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gj_xj = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->gj_xj.x, gp->gj_xj.y, gp->gj_xj.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gk    =  %10.6f\n", i, j, k, c, csr2->gj);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gk_xk = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->gk_xk.x, gp->gk_xk.y, gp->gk_xk.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gvk   = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gvk.x, gvk.y, gvk.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gvkx_xk={%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gvkx_xk.x/q.x, gvkx_xk.y/q.x, gvkx_xk.z/q.x);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gvky_xk={%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gvky_xk.x/q.y, gvky_xk.y/q.y, gvky_xk.z/q.y);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: gvkz_xk={%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gvkz_xk.x/csr1->gj, gvkz_xk.y/csr1->gj, gvkz_xk.z/csr1->gj);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: s     =  %10.6f\n", i, j, k, c, atan2(q.y, q.x));
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: s_xj  = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->s_xj.x, gp->s_xj.y, gp->s_xj.z);
		fprintf(fp, "CALC_GRADIENT[%d,%d,%d,%c]: s_xk  = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, gp->s_xk.x, gp->s_xk.y, gp->s_xk.z);
	}
#endif
}			

/*  Find the intersection point with the same point_id  */
static int
s_find_partner_point(MDArena *arena, const SPArena *sarena, const intersect_record *irp)
{
	Int i = irp->i;
	Int start = sarena->sphere_idx[i * 2];
	Int end = sarena->sphere_idx[i * 2 + 1];
	Int j, n;
	const crossing_sphere_record *csj;
	for (j = start; j < end; j++) {
		if (j == irp->idx)
			continue;
		csj = &sarena->spheres[j];
		for (n = csj->int_start; n < csj->int_end; n++) {
			if (sarena->intersects[n].point_id == irp->point_id)
				break;
		}
		if (n < csj->int_end)
			return n;
	}
	md_warning(arena, "Internal inconsistency: the intersection point (%d,%d,%d,%d) appears only once", irp->i+1, irp->j+1, irp->k+1, irp->point_id);
#if DEBUG
	for (j = 0; j < sarena->nspheres; j++) {
		csj = &sarena->spheres[j];
		for (n = csj->int_start; n < csj->int_end; n++) {
			intersect_record *irn = &sarena->intersects[n];
			md_warning(arena, "DEBUG: point %d = (%d,%d,%d;%d;%d), sphere %d = (%d,%d)", n, irn->i+1, irn->j+1, irn->k+1, irn->point_id, irn->idx, j, csj->i+1, csj->j+1);
			if (irn->i != csj->i || irn->j != csj->j)
				md_warning(arena, "DEBUG: internal inconsistency, point %d does not belong to the correct crossing sphere record", n);
		}
	}
#endif			
	return -1;
}

/*  Print the debug information  */
#if SP_DEBUG
static void
s_print_surface_area_debug(const MDArena *arena, const gradient_record *gp, const intersect_record *ip, Double da, const char *tag)
{
	if (arena->debug_result != NULL && arena->debug_output_level > 2) {
		Int i, j, k, c;
		FILE *fp = arena->debug_result;  /*  Just an alias for shorter typing */
		i = ip->i + 1;
		j = ip->j + 1;
		k = ip->k + 1;
		c = (gp->sign ? '+' : '-');
		fprintf(fp, "SURF_AREA_GRAD[%d,%d,%d,%c]: %-6s: da   =  %10.6f\n", i, j, k, c, tag, da);
		fprintf(fp, "SURF_AREA_GRAD[%d,%d,%d,%c]: %-6s: a_xj = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, tag, gp->a_xj.x, gp->a_xj.y, gp->a_xj.z);
		fprintf(fp, "SURF_AREA_GRAD[%d,%d,%d,%c]: %-6s: a_xk = {%10.6f,%10.6f,%10.6f}\n", i, j, k, c, tag, gp->a_xk.x, gp->a_xk.y, gp->a_xk.z);
	}
}
#define SA_DEBUG(tag, ir) s_print_surface_area_debug(arena, &grad, ir, da, tag)
#else
#define SA_DEBUG(tag, ir)
#endif

/*  Calculate the surface area for each atom  */
static void
s_calc_surface_area(MDArena *arena)
{
	Int i;      /*  The atom index  */
	Int idx1, idx2;   /*  The circle (or crossing sphere) indices  */
	Int n1, n2; /*  The intersection point index  */
	Int count;  /*  The number of intersection points processed  */
	intersect_record *ir1, *ir2;
	Vector p1, p2;
	gradient_record grad;
	Double cosomega, sinomega, sinomega_inv, omega, ss;
	Int natoms = arena->mol->natoms;
	SPArena *sarena = arena->sp_arena;
	Double *energies = &arena->energies[kSurfaceIndex];
	Vector *forces = &arena->forces[kSurfaceIndex * arena->mol->natoms];

#if SP_DEBUG
	s_sp_debug(arena, "SURFACE AREA CALCULATION");
#endif
	memset(forces, 0, sizeof(Vector) * arena->mol->natoms);
	for (i = 0; i < natoms; i++) {
		Int start, end;
		Double area = 0.0;  /*  Accessible area  */
		Double rho = sarena->atom_rad[i];  /*  The radius of atom i  */
		Double da;
		crossing_sphere_record *csr1, *csr2;
		start = sarena->sphere_idx[i * 2];
		end = sarena->sphere_idx[i * 2 + 1];
		/*  Treat the special cases  */
		if (sarena->atom_buried[i]) {
			sarena->atom_area[i] = 0.0;
		#if SP_DEBUG
			s_sp_debug(arena, "atom_area[%d] = 0.0, because atom[%d] is buried", i+1, i+1);
		#endif
			continue;
		}
		if (start == end) {
			/*  No crossing spheres (isolated sphere)  */
			sarena->atom_area[i] = 4 * PI * rho * rho;
		#if SP_DEBUG
			s_sp_debug(arena, "atom_area[%d] = 4*PI*rho**2, because atom[%d] is isolated", i+1, i+1);
		#endif
			continue;
		}
		/*  Find the circles without any intersection points */
		for (idx1 = start; idx1 < end; idx1++) {
			csr1 = &sarena->spheres[idx1];
			if (csr1->skip)
				continue;
			if (csr1->int_end == csr1->int_start) {
				/*  Is this circle buried? */
			/*	p1.x = csr1->rj;
				p1.y = 0.0;
				p1.z = csr1->gj;
				quat_rotate(&p1, &p1, &csr1->quat_c);
				for (idx2 = start; idx2 < end; idx2++) {
					if (idx1 == idx2)
						continue;
					if (s_is_inside(sarena, &sarena->spheres[idx2], &p1))
						break;
				}
				if (idx2 < end) {
				} else {
				#if SP_DEBUG
					s_sp_debug(arena, "Circle (%d,%d) is buried", i+1, csr1->j+1);
				#endif
				}
			*/
				Double w1 = sarena->atom_pot[i];
			/*	Double w2 = sarena->atom_rad[csr1->j];  */
				da = 2 * PI * (1.0 + csr1->gj / rho);
				area += da;
			/*	w1 = w1 * (csr1->dij * csr1->dij + w2 * w2 - rho * rho) / (2.0 * csr1->dij * csr1->dij * csr1->dij); */
				w1 = w1 * PI * rho / csr1->dij;
				VecScale(p1, csr1->vij, w1);
				VecDec(forces[csr1->j], p1);
				VecInc(forces[i], p1);
			#if SP_DEBUG
				s_sp_debug(arena, "Circle (%d,%d) forms a boundary itself, dA/(rho^2) = %f", i+1, csr1->j+1, da);
			#endif
			}
		}
		/*  Traverse the set of arcs divided by the intersection points  */
		n1 = -1;
		ir1 = NULL;
		csr1 = NULL;
		count = 0;
		while (1) {
			Vector q1, pw1, pw2;
			Double w1, w2, w3, w4;
			int do_grad = 1;
			if (n1 < 0) {
				/*  Find the first unprocessed intersection point  */
				ir1 = NULL;
			#if SP_DEBUG
				s_sp_debug(arena, "Trying to find a new intersection point...");
			#endif
				for (idx1 = start; idx1 < end; idx1++) {
					csr1 = &sarena->spheres[idx1];
					for (n1 = csr1->int_start; n1 < csr1->int_end; n1++) {
						if (sarena->intersects[n1].flag == 0) {
							ir1 = &sarena->intersects[n1];
							break;
						}
					}
					if (ir1 != NULL)
						break;
				}
				if (ir1 == NULL) {
				#if SP_DEBUG
					s_sp_debug(arena, "All points have been processed.");
				#endif
					break;  /*  All points have been processed  */
				}
			#if SP_DEBUG
				s_sp_debug(arena, "A new intersection point (%d,%d,%d;%d) was found.", ir1->i+1, ir1->j+1, ir1->k+1, ir1->point_id);
			#endif
				/*  Determine whether this point is starting point or ending point */
				/*  ss is the point between this point and the next point  */
				if (n1 == csr1->int_end - 1) {
					ss = 0.5 * (ir1->s + sarena->intersects[csr1->int_start].s + 2 * PI);
				} else {
					ss = 0.5 * (ir1->s + (ir1 + 1)->s);
				}
				p1.x = csr1->rj * cos(ss);
				p1.y = csr1->rj * sin(ss);
				p1.z = csr1->gj;
				MatrixVec(&p1, csr1->rot_c, &p1); /*quat_rotate(&p1, &p1, &csr1->quat_c);*/
			#if SP_DEBUG
				s_sp_debug(arena, "Examining the point for s=%f, {%f,%f,%f}, whether it is inside any circle", ss, p1.x, p1.y, p1.z);
			#endif
				/*  Is this point in any of the circle on the same sphere? */
				for (idx2 = start; idx2 < end; idx2++) {
					if (idx2 == idx1) /*|| sarena->spheres[idx2].j == ir1->k)*/
						continue;
					if (s_is_inside(&sarena->spheres[idx2], &p1)) {
					#if SP_DEBUG
						s_sp_debug(arena, "The point is within the circle (%d,%d)", sarena->spheres[idx2].i+1, sarena->spheres[idx2].j+1);
					#endif
						break;
					}
				}
				/*  If yes, then this is the end point. Otherwise, it is the start point. */
				/*  We start traversing from the end point; end point->next arc's start
				 *  point->traverse the arc->end point->next arc's start point->... */
				if (idx2 >= end) {
					/*  It is the start point, so start from the partner point */
					n1 = s_find_partner_point(arena, sarena, ir1);
					if (n1 >= 0) {
						ir1 = &sarena->intersects[n1];
						idx1 = ir1->idx;
						csr1 = &sarena->spheres[idx1];
					}
				}
			#if SP_DEBUG
				s_sp_debug(arena, "Start traversing from point (%d,%d,%d;%d)", ir1->i+1, ir1->j+1, ir1->k+1, ir1->point_id);
			#endif
			} /* end if (n1 < 0) */
			if (n1 < 0)
				break; /* Abnormal end */

			/*  Now we are at the 'end point' of a particular arc. */

			/*  Find the partner point, which is the 'start point' of the next arc. */
			n2 = s_find_partner_point(arena, sarena, ir1);
			if (n2 < 0)
				break;  /* Abnormal end */
			ir2 = &sarena->intersects[n2];
			idx2 = ir2->idx;
			csr2 = &sarena->spheres[idx2];
		#if SP_DEBUG
			s_sp_debug(arena, "The partner point is (%d,%d,%d;%d).", ir2->i+1, ir2->j+1, ir2->k+1, ir2->point_id);
		#endif

			/*  The tangent vector (p1) and position vector (q1) at this point  */
			/*  (in global coordinates)  */
			p1.x = -sin(ir1->s);
			p1.y = cos(ir1->s);
			p1.z = 0.0;
			q1.y = -csr1->rj * p1.x;
			q1.x = csr1->rj * p1.y;
			q1.z = csr1->gj;
			MatrixVec(&p1, csr1->rot_c, &p1); /*quat_rotate(&p1, &p1, &csr1->quat_c);*/
			MatrixVec(&q1, csr1->rot_c, &q1); /*quat_rotate(&q1, &q1, &csr1->quat_c);*/

			/*  The tangent vector at the partner point  */
			p2.x = -sin(ir2->s);
			p2.y = cos(ir2->s);
			p2.z = 0.0;
			MatrixVec(&p2, csr2->rot_c, &p2); /*quat_rotate(&p2, &p2, &csr2->quat_c);*/

			/*  The external angle (omega)  */
			cosomega = VecDot(p1, p2)/(VecLength(p1) * VecLength(p2));
			sinomega = sqrt(1.0 - cosomega * cosomega);
			omega = atan2(sinomega, cosomega);

			/*  If omega is close to 0 or PI, then avoid calculation of the gradient
			 *  because the surface area function becomes singular  */
			if (cosomega > 0.95 || cosomega < -0.95) {
				VecZero(grad.a_xj);
				VecZero(grad.a_xk);
				do_grad = 0;
			} else do_grad = 1;

			/*  Gradients of s (arc parameter) at this point */
		#if SP_DEBUG
			grad.sign = 0;
		#endif
			if (do_grad)
				s_calc_gradient_vectors(&grad, arena, csr1, csr2, &q1);

			/*  The arc parameter term of the surface area */
			/*  (The arc parameter term is added at each 'end point', and
			 *  subtracted at each 'start point'.)  */
			da = ir1->s * csr1->gj / rho;
		#if SP_DEBUG
			s_sp_debug(arena, "The arc term [end point] = %f (s = %f)", ir1->s * csr1->gj / rho, ir1->s);
		#endif

			/*  Gradients of the arc parameter term  */
			if (do_grad) {
				grad.a_xj.x = (csr1->gj * grad.s_xj.x + ir1->s * grad.gj_xj.x) / rho;
				grad.a_xj.y = (csr1->gj * grad.s_xj.y + ir1->s * grad.gj_xj.y) / rho;
				grad.a_xj.z = (csr1->gj * grad.s_xj.z + ir1->s * grad.gj_xj.z) / rho;
				grad.a_xk.x = (csr1->gj * grad.s_xk.x) / rho;
				grad.a_xk.y = (csr1->gj * grad.s_xk.y) / rho;
				grad.a_xk.z = (csr1->gj * grad.s_xk.z) / rho;
				SA_DEBUG("ARC1", ir1);
			}

			/*  The omega (external angle) term of the surface area  */
			VecCross(pw2, p1, p2);
			if (VecDot(q1, pw2) < 0) {
				omega = -omega;
				sinomega = -sinomega;
			}
			da += omega;
		#if SP_DEBUG
			s_sp_debug(arena, "The external angle is %f (%f deg); the tangent vectors are {%f,%f,%f}, {%f,%f,%f}", omega, omega*180/PI, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
		#endif

			sinomega_inv = 0.0;  /*  Not significant; just for suppress compiler warnings */

			/*  pw1(pw2): vector from the center of the circle j(k) to q1 */
			w1 = -csr1->gj / csr1->dij;
			pw1 = q1;
			VecScaleInc(pw1, csr1->vij, w1);
			w1 = -csr2->gj / csr2->dij;
			pw2 = q1;
			VecScaleInc(pw2, csr2->vij, w1);

			/*  Gradients of the omega term, the first half  */
			if (do_grad) {
				sinomega_inv = (fabs(sinomega) < 1e-5 ? 0.0 : 1.0 / sinomega);
				w2 = sinomega_inv * VecDot(pw1, p2) / csr1->rj;
				w3 = -sin(ir1->s) * sinomega_inv;
				w4 = cos(ir1->s) * sinomega_inv;
				grad.a_xj.x += w2 * grad.s_xj.x
					- w3 * (grad.mx[0]*p2.x + grad.mx[1]*p2.y + grad.mx[2]*p2.z)
					- w4 * (grad.mx[3]*p2.x + grad.mx[4]*p2.y + grad.mx[5]*p2.z);
				grad.a_xj.y += w2 * grad.s_xj.y
					- w3 * (grad.my[0]*p2.x + grad.my[1]*p2.y + grad.my[2]*p2.z)
					- w4 * (grad.my[3]*p2.x + grad.my[4]*p2.y + grad.my[5]*p2.z);
				grad.a_xj.z += w2 * grad.s_xj.z
					- w3 * (grad.mz[0]*p2.x + grad.mz[1]*p2.y + grad.mz[2]*p2.z)
					- w4 * (grad.mz[3]*p2.x + grad.mz[4]*p2.y + grad.mz[5]*p2.z);
				VecScaleInc(grad.a_xk, grad.s_xk, w2);
			}

		#if 0
		{
			s_sp_debug(arena, "[%d,%d,%d] p1.x = %10.6f, d(p1.x)/dxj = {%10.6f,%10.6f,%10.6f} d(p1.x)/dxk = {%10.6f,%10.6f,%10.6f}", ir1->i+1, ir1->j+1, ir1->k+1, p1.x,
				(grad.mx[0]*w3+grad.mx[3]*w4)*sinomega - pw1.x/csr1->rj*grad.s_xj.x,
				(grad.my[0]*w3+grad.my[3]*w4)*sinomega - pw1.x/csr1->rj*grad.s_xj.y,
				(grad.mz[0]*w3+grad.mz[3]*w4)*sinomega - pw1.x/csr1->rj*grad.s_xj.z,
				-pw1.x/csr1->rj*grad.s_xk.x,
				-pw1.x/csr1->rj*grad.s_xk.y,
				-pw1.x/csr1->rj*grad.s_xk.z);
			s_sp_debug(arena, "[%d,%d,%d] p1.y = %10.6f, d(p1.y)/dxj = {%10.6f,%10.6f,%10.6f} d(p1.y)/dxk = {%10.6f,%10.6f,%10.6f}", ir1->i+1, ir1->j+1, ir1->k+1, p1.y,
				(grad.mx[1]*w3+grad.mx[4]*w4)*sinomega - pw1.y/csr1->rj*grad.s_xj.x,
				(grad.my[1]*w3+grad.my[4]*w4)*sinomega - pw1.y/csr1->rj*grad.s_xj.y,
				(grad.mz[1]*w3+grad.mz[4]*w4)*sinomega - pw1.y/csr1->rj*grad.s_xj.z,
				-pw1.y/csr1->rj*grad.s_xk.x,
				-pw1.y/csr1->rj*grad.s_xk.y,
				-pw1.y/csr1->rj*grad.s_xk.z);
			s_sp_debug(arena, "[%d,%d,%d] p1.z = %10.6f, d(p1.z)/dxj = {%10.6f,%10.6f,%10.6f} d(p1.z)/dxk = {%10.6f,%10.6f,%10.6f}", ir1->i+1, ir1->j+1, ir1->k+1, p1.y,
				(grad.mx[2]*w3+grad.mx[5]*w4)*sinomega - pw1.z/csr1->rj*grad.s_xj.x,
				(grad.my[2]*w3+grad.my[5]*w4)*sinomega - pw1.z/csr1->rj*grad.s_xj.y,
				(grad.mz[2]*w3+grad.mz[5]*w4)*sinomega - pw1.z/csr1->rj*grad.s_xj.z,
				-pw1.z/csr1->rj*grad.s_xk.x,
				-pw1.z/csr1->rj*grad.s_xk.y,
				-pw1.z/csr1->rj*grad.s_xk.z);
		}
		#endif

			/*  Calculate the gradients for the partner point */
		#if SP_DEBUG
			grad.sign = 1;
		#endif
			if (do_grad)
				s_calc_gradient_vectors(&grad, arena, csr2, csr1, &q1);

			/*  Gradients of the omega term, the second half  */
			if (do_grad) {
				w2 = sinomega_inv * VecDot(pw2, p1) / csr2->rj;
				w3 = -sin(ir2->s) * sinomega_inv;
				w4 = cos(ir2->s) * sinomega_inv;
				grad.a_xk.x += w2 * grad.s_xj.x
					- w3 * (grad.mx[0]*p1.x + grad.mx[1]*p1.y + grad.mx[2]*p1.z)
					- w4 * (grad.mx[3]*p1.x + grad.mx[4]*p1.y + grad.mx[5]*p1.z);
				grad.a_xk.y += w2 * grad.s_xj.y
					- w3 * (grad.my[0]*p1.x + grad.my[1]*p1.y + grad.my[2]*p1.z)
					- w4 * (grad.my[3]*p1.x + grad.my[4]*p1.y + grad.my[5]*p1.z);
				grad.a_xk.z += w2 * grad.s_xj.z
					- w3 * (grad.mz[0]*p1.x + grad.mz[1]*p1.y + grad.mz[2]*p1.z)
					- w4 * (grad.mz[3]*p1.x + grad.mz[4]*p1.y + grad.mz[5]*p1.z);
				VecScaleInc(grad.a_xj, grad.s_xk, w2);
				SA_DEBUG("OMEGA", ir1);
			}

		#if 0
		{
			s_sp_debug(arena, "[%d,%d,%d] p2.x = %10.6f, d(p2.x)/dxj = {%10.6f,%10.6f,%10.6f} d(p2.x)/dxk = {%10.6f,%10.6f,%10.6f}", ir1->i+1, ir1->j+1, ir1->k+1, p2.x,
				-pw2.x/csr2->rj*grad.s_xk.x,
				-pw2.x/csr2->rj*grad.s_xk.y,
				-pw2.x/csr2->rj*grad.s_xk.z,
				(grad.mx[0]*w3+grad.mx[3]*w4)*sinomega - pw2.x/csr2->rj*grad.s_xj.x,
				(grad.my[0]*w3+grad.my[3]*w4)*sinomega - pw2.x/csr2->rj*grad.s_xj.y,
				(grad.mz[0]*w3+grad.mz[3]*w4)*sinomega - pw2.x/csr2->rj*grad.s_xj.z);
			s_sp_debug(arena, "[%d,%d,%d] p2.y = %10.6f, d(p2.y)/dxj = {%10.6f,%10.6f,%10.6f} d(p2.y)/dxk = {%10.6f,%10.6f,%10.6f}", ir1->i+1, ir1->j+1, ir1->k+1, p2.y,
				-pw2.y/csr2->rj*grad.s_xk.x,
				-pw2.y/csr2->rj*grad.s_xk.y,
				-pw2.y/csr2->rj*grad.s_xk.z,
				(grad.mx[1]*w3+grad.mx[4]*w4)*sinomega - pw2.y/csr2->rj*grad.s_xj.x,
				(grad.my[1]*w3+grad.my[4]*w4)*sinomega - pw2.y/csr2->rj*grad.s_xj.y,
				(grad.mz[1]*w3+grad.mz[4]*w4)*sinomega - pw2.y/csr2->rj*grad.s_xj.z);
			s_sp_debug(arena, "[%d,%d,%d] p2.z = %10.6f, d(p2.z)/dxj = {%10.6f,%10.6f,%10.6f} d(p2.z)/dxk = {%10.6f,%10.6f,%10.6f}", ir1->i+1, ir1->j+1, ir1->k+1, p2.z,
				-pw2.z/csr2->rj*grad.s_xk.x,
				-pw2.z/csr2->rj*grad.s_xk.y,
				-pw2.z/csr2->rj*grad.s_xk.z,
				(grad.mx[2]*w3+grad.mx[5]*w4)*sinomega - pw2.z/csr2->rj*grad.s_xj.x,
				(grad.my[2]*w3+grad.my[5]*w4)*sinomega - pw2.z/csr2->rj*grad.s_xj.y,
				(grad.mz[2]*w3+grad.mz[5]*w4)*sinomega - pw2.z/csr2->rj*grad.s_xj.z);
		}
		#endif

			/*  The arc parameter term for the partner point  */
			/*  Substract 2*PI if this is the last point registered for circle k,
			 *  because the s parameter of next point will cross the [-PI,PI] boundary. */
			ss = ir2->s;
			if (n2 == csr2->int_end - 1)
				ss -= 2 * PI;
			da -= ss * csr2->gj / rho;
		#if SP_DEBUG
			s_sp_debug(arena, "The arc term [start point] = %f (s = %f%s)", ss * csr2->gj / rho, ir2->s, (n2 == csr2->int_end - 1 ? " - 2*PI" : ""));
		#endif
			/*  Gradients of the arc parameter term  */
			if (do_grad) {
				grad.a_xk.x -= (csr2->gj * grad.s_xj.x + ss * grad.gj_xj.x) / rho;
				grad.a_xk.y -= (csr2->gj * grad.s_xj.y + ss * grad.gj_xj.y) / rho;
				grad.a_xk.z -= (csr2->gj * grad.s_xj.z + ss * grad.gj_xj.z) / rho;
				grad.a_xj.x -= (csr2->gj * grad.s_xk.x) / rho;
				grad.a_xj.y -= (csr2->gj * grad.s_xk.y) / rho;
				grad.a_xj.z -= (csr2->gj * grad.s_xk.z) / rho;
				SA_DEBUG("ARC2", ir1);
			}

			/*  Mark this point (two entries) as processed  */
			ir1->flag = ir2->flag = 1;

			/*  Add up the partial area and force  */
			w1 = sarena->atom_pot[ir1->i] * rho * rho;
			area += da;
			VecScale(pw1, grad.a_xj, w1);
			VecScale(pw2, grad.a_xk, w1);
			VecDec(forces[ir1->j], pw1);
			VecDec(forces[ir1->k], pw2);
			VecInc(pw1, pw2);
			VecInc(forces[ir1->i], pw1);

		#if DEBUG
			if (VecLength2(pw1) > 1.0 || VecLength2(pw2) > 1.0) {
				md_warning(arena, "Warning: The surface force becomes very large. Maybe something is wrong?");
				md_warning(arena, "Warning: STEP %d, (i,j,k,id)=(%d,%d,%d,%d), pw1={%f,%f,%f}, pw2={%f,%f,%f}", arena->step, ir1->i+1, ir1->j+1, ir1->k+1, ir1->point_id, pw1.x, pw1.y, pw1.z, pw2.x, pw2.y, pw2.z);
			}
		#endif
		
			count++;

			/*  Find the next point on circle k */
			n1 = n2 + 1;
			if (n1 >= csr2->int_end)
				n1 = csr2->int_start;
			ir1 = &sarena->intersects[n1];
		#if SP_DEBUG
			s_sp_debug(arena, "The next point on the arc is (%d,%d,%d;%d)%s", ir1->i+1, ir1->j+1, ir1->k+1, ir1->point_id, (ir1->flag ? ", which is already processed" : "."));
		#endif
			if (ir1->flag) {
				/*  This point is already processed; then close this boundary  */
				n1 = -1;
			} else {
				idx1 = idx2;
				csr1 = csr2;
			}
			

		} /* end while (1) */
		if (count > 0) { /* Some intersection points were processed */
			area += 2 * PI;
		}
		area -= floor(area / (4 * PI)) * (4 * PI);
		sarena->atom_area[i] = area * rho * rho;
		*energies += sarena->atom_area[i] * sarena->atom_pot[i];
	#if SP_DEBUG
		s_sp_debug(arena, "The total accessible area of atom %d is %f", i+1, area*rho*rho);
	#endif
	}
}

/*  Calculate the surface potential/force of a molecule */
void
calc_surface_force(MDArena *arena)
{
	if (arena->sp_arena == NULL)
		s_allocate_sp_arena(arena);
	arena->sp_arena->nspheres = 0;
	arena->sp_arena->nintersects = 0;
	arena->sp_arena->npoints = 0;

	s_list_crossing_spheres(arena);
	s_list_intersection_points(arena);
#if SP_DEBUG
	s_print_internal_information(arena);
#endif
	s_calc_surface_area(arena);

#if DEBUG
	if (arena->debug_result && arena->debug_output_level > 1) {
		int i;
		for (i = 0; i < arena->mol->natoms; i++) {
			Vector *fp = &arena->forces[kSurfaceIndex * arena->mol->natoms + i];
			fprintf(arena->debug_result, "surface of atom %d: area=%f, pot=%f, {%f %f %f}\n", i+1, arena->sp_arena->atom_area[i], arena->sp_arena->atom_area[i]*arena->sp_arena->atom_pot[i], fp->x * INTERNAL2KCAL, fp->y * INTERNAL2KCAL, fp->z * INTERNAL2KCAL);
		}
	}
#endif
}

/*  Clear the surface potential record  */
void
clear_sp_arena(MDArena *arena)
{
	if (arena->sp_arena != NULL) {
		SPArena *sarena = arena->sp_arena;
		if (sarena->spheres != NULL)
			free(sarena->spheres);
		if (sarena->intersects != NULL)
			free(sarena->intersects);
		if (sarena->points != NULL)
			free(sarena->points);
		if (sarena->atom_rad != NULL)
			free(sarena->atom_rad);
		if (sarena->atom_pot != NULL)
			free(sarena->atom_pot);
		if (sarena->atom_area != NULL)
			free(sarena->atom_area);
		if (sarena->atom_buried != NULL)
			free(sarena->atom_buried);
		if (sarena->sphere_idx != NULL)
			free(sarena->sphere_idx);
		free(sarena);
		arena->sp_arena = NULL;
	}
}

/*  Print surface area information  */
void
print_surface_area(MDArena *arena)
{
	Int i, natoms;
	if (arena->sp_arena == NULL)
		return;
	natoms = arena->mol->natoms;
	printf("Atom     area    energy force\n");
	for (i = 0; i < natoms; i++) {
		Vector f = arena->forces[kSurfaceIndex * natoms + i];
		Double area = arena->sp_arena->atom_area[i];
		Double energy = area * arena->sp_arena->atom_pot[i];
		VecScaleSelf(f, INTERNAL2KCAL);
		energy *= INTERNAL2KCAL;
		printf("%5d %11.5f %11.5f %11.5f %11.5f %11.5f\n",
			i+1, area, energy, f.x, f.y, f.z);
	}
}
