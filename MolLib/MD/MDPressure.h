/*
 *  MDPressure.h
 *  Molby
 *
 *  Created by Toshi Nagata on 07/10/12.
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

#ifndef __MDPRESSURE_H__
#define __MDPRESSURE_H__

#include "MDCore.h"

#ifdef __cplusplus
extern "C" {
#endif
	
struct MDPressureArena {
	Int    disabled;
	Int    freq;
	Double  coupling;
	Double  trial_width;    /* Move width for MC trials. Default 0.01.  */
	Mat33  apply;          /*  A tensor. The unit is bar (NOT internal unit)  */
	Int    mc_accept, mc_reject;  /*  Monte Carlo statistics  */
	Double  cell_flexibility[8];  /*  How the cell vectors are modified in each trial step  */
	/*  cell_flexibility[0..2] determines whether the cell vector is modified along its
	 axis. If the value is 0, the cell vector is not changed in that direction. If the
	 value is positive, the cell vector is changed. If the value is negative, the cell 
	 vector is changed, but the amount of change is the same as the last axis which has
	 a positive flexibility value.
	 cell_flexibility[3..5] determines whether two axes interact each other. The
	 index 3, 4, 5 corresponds to b/c, c/a, a/b pair, respectively. If the value is
	 0, then the two axes do not interact. If the value is positive, then the two axes are
	 modified along the other axis; for b/c pair, they will be modifed as b->b + xc,
	 c->c + xb. To avoid unphysical rotation of the axes system, the same modifier
	 coefficient x is used for both transformation.
	 cell_flexibility[6] allows (if non-zero) the cell origin to fluctuate.
	 pressure_fluctuate_cell_origin can be set to control the amount of fluctuation.
	 cell_flexibility[7] allows the cell orientation to fluctuate. 1: rotation
	 around x axis, 2; rotation around y axis, 3: rotation around z axis, 4: free
	 rotation. The amount of fluctuation is controlled by pressure_fluctuate_cell_orientation.
	 The absolute value of cell_flexibility[0..5] are multiplied to trial_width, so that
	 different trial_width for the cell parameters are allowed. */
	Int    control_algorithm;
	Double  *temporary_energies;  /*  Double[natoms_uniq]  */
	Vector *temporary_velocities; /*  Vector[natoms_uniq]  */
	Int    use_atomic_move;  /*  If non-zero, all atoms positions are scaled. Otherwise, only the center of mass of each molecule is scaled.  */
	Double  fluctuate_cell_origin;  /*  Default 0.01. */
	Double  fluctuate_cell_orientation; /*  Default 0.01.  */
	Double  mc_delta_potential;
	Double  mc_delta_kinetic;
	Double  mc_stress;
	Double  mc_delta_volume;
};

typedef struct MDPressureArena MDPressureArena;

MDPressureArena *pressure_new(void);
void pressure_prepare(MDArena *arena);
void pressure_release(MDPressureArena *pressure);
void pressure_control(MDArena *arena);

#ifdef __cplusplus
}
#endif
		
#endif /* __MDPRESSURE_H__ */
