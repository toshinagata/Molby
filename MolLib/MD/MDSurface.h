/*
 *  MDSurface.h
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

#ifndef __MDSURFACE_H__
#define __MDSURFACE_H__

#include "MDCore.h"

#ifdef __cplusplus
extern "C" {
#endif
	
/*  Internal information for a pair of intercrossing spheres  */
typedef struct crossing_sphere_record {
	Int    i;        /*  The index of the target sphere  */
	Int    j;        /*  The index of the crossing sphere  */
	Vector vij;      /*  Center-to-center vector  */
	Double  dij;      /*  The length of vij  */
	Double  rj;       /*  The radius of the intersection circle  */
	Double  gj;       /*  sqrt(atom_rad[i]**2 - rj**2)  */
	Mat33  rot;      /*  The rotation matrix that convert vij->(0,0,1) and qij->(1,0,0) 
	                     where qij = vij x (1,0,0) {if vij is not parallel with (1,0,0)} 
	                              or vij x (0,1,0) {otherwise}  */
	Mat33  rot_c;    /*  The inverse rotation  */
	Int    skip;     /*  This entry should be skipped (i.e. the intersection circle is included within another circle) */
	Int    int_start, int_end;  /* The start and end indices of the intersect_record array (int_end is not included, so int_end - int_start is the number of records) */
} crossing_sphere_record;

/*  Intersection points on one particular crossing circle  */
typedef struct intersect_record {
	Int    i;        /*  The index of the target sphere  */
	Int    j;        /*  The index of the sphere to define the crossing circle  */
	Int    k;        /*  The index of the sphere that causes this intersection point */
	Int    idx;      /*  The index to the spheres[] (to allow quick access) */
	unsigned int flag : 1;      /*  The flag to show this point is already processed  */
	unsigned int point_id : 31; /*  The id number of this point = index to points[] */
	Double  s;        /*  The arc-length parameter (radian in [-pi,pi]) */
} intersect_record;

/*  One particular intersection point  */
typedef struct intersect_point_record {
	unsigned int sign : 1;  /*  The flag to distinguish the two crossing points on the same (i,j,k) sphere combination  */
	unsigned int i : 31;    /*  The index of the target sphere  */
	Int    idx_j, idx_k; /*  The crossing circles; indices to spheres[]  */
	Vector v;        /*  The global coordinate  */
} intersect_point_record;

/*  A record for working on gradients  */
typedef struct gradient_record {
	Vector qx_xj, qy_xj, gj_xj, s_xj;
	Vector qx_xk, qy_xk, gk_xk, s_xk;
	Mat33  mx, my, mz;  /*  Differentials of (i,j)-local-to-global matrix M by xj.x, xj.y, xj.z  */
	Vector a_xj, a_xk;
#if DEBUG
	/*  Not significant for calculation; only used for printing debug information. */
	Byte sign;  /* 1 if this point is 'startint point', 0 otherwise. */
#endif
} gradient_record;

typedef struct SPArena {
	Int    nspheres;
	Int    maxspheres;
	crossing_sphere_record *spheres;
	Int    nintersects;
	Int    maxintersects;
	intersect_record *intersects;
	Int    npoints;
	Int    maxpoints;
	intersect_point_record *points;    /*  The global coordinates of the intersection points  */
	Double  *atom_rad;  /*  The radius of the atoms (sum of the vdw and probe radii)  */
	Double  *atom_pot;  /*  The contribution of each atom to the total surface potential */
	Double  *atom_area; /*  The accessible surface of each atom  */
	Byte   *atom_buried;  /*  1 if this atom is buried in another atom  */
	Int    *sphere_idx; /*  The start/end indices of the crossing_sphere_record array (the size of the array is natoms*2) */
} SPArena;

void clear_sp_arena(MDArena *arena);
void calc_surface_force(MDArena *arena);
void print_surface_area(MDArena *arena);

#ifdef __cplusplus
}
#endif
		
#endif /* __MDSURFACE_H__ */

