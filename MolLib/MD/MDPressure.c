/*
 *  MDPressure.c
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

#include "MDPressure.h"
#include "MDForce.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

/*  Calculate the lattice deformation stress  */
static Double
s_lattice_deformation_stress(const Mat33 apply, const Transform celltr, const Transform old_celltr)
{
	Mat33 grad, grad_tr, grad_inv, grad_inv_tr;
	Mat33 strain, tension;
	Double eng, volume;
	Int i;

	/*  Some of the calculations duplicate with md_scale_cell(), but it will not cause
		a serious slowdown  */
	volume = fabs(MatrixDeterminant(celltr));

	/*  Deformation gradient (temp) = celltr * old_rcelltr  */
	MatrixInvert(grad, old_celltr);
	MatrixMul(grad, celltr, grad);
	MatrixTranspose(grad_tr, grad);
	MatrixInvert(grad_inv, grad);
	MatrixTranspose(grad_inv_tr, grad_inv);
	MatrixMul(strain, grad_tr, grad);
	strain[0] -= 1.0;
	strain[4] -= 1.0;
	strain[8] -= 1.0;
	MatrixScale(strain, strain, 0.5);
	MatrixMul(tension, apply, grad_inv_tr);
	MatrixMul(tension, grad_inv, tension);
	MatrixScale(tension, tension, MatrixDeterminant(grad));
	eng = 0.0;
	for (i = 0; i < 9; i++)
		eng += strain[i] * tension[i];
	eng *= volume * BAR2INTERNALPRESSURE;
	return eng;
}

MDPressureArena *
pressure_new(void)
{
	MDPressureArena *pressure = (MDPressureArena *)calloc(sizeof(MDPressureArena), 1);
	if (pressure == NULL)
		return NULL;
	pressure->freq = 20;
	pressure->coupling = 0.3;
	pressure->trial_width = 0.001;
	pressure->control_algorithm = 0;
	pressure->fluctuate_cell_origin = 0.01;
	pressure->fluctuate_cell_orientation = 0.01;
	pressure->apply[0] = 1;
	pressure->apply[1] = 0;
	pressure->apply[2] = 0;
	pressure->apply[3] = 0;
	pressure->apply[4] = 1;
	pressure->apply[5] = 0;
	pressure->apply[6] = 0;
	pressure->apply[7] = 0;
	pressure->apply[8] = 1;
	pressure->cell_flexibility[0] = -1;
	pressure->cell_flexibility[1] = -1;
	pressure->cell_flexibility[2] = -1;
	pressure->cell_flexibility[3] = 0;
	pressure->cell_flexibility[4] = 0;
	pressure->cell_flexibility[5] = 0;
	pressure->cell_flexibility[6] = 0;
	pressure->cell_flexibility[7] = 0;
	return pressure;
}

void
pressure_release(MDPressureArena *pressure)
{
	if (pressure == NULL)
		return;
	free(pressure);
}

void
pressure_prepare(MDArena *arena)
{
	MDPressureArena *pressure = arena->pressure;
	if (pressure == NULL)
		return;
	pressure->mc_accept = pressure->mc_reject = 0;
	if (pressure->temporary_energies != NULL)
		free(pressure->temporary_energies);
	pressure->temporary_energies = (Double *)calloc(sizeof(Double), arena->natoms_uniq);
	if (pressure->temporary_energies == NULL)
		md_panic(arena, "Low memory");
	if (pressure->temporary_velocities != NULL)
		free(pressure->temporary_velocities);
	pressure->temporary_velocities = (Vector *)calloc(sizeof(Vector), arena->natoms_uniq);
	if (pressure->temporary_velocities == NULL)
		md_panic(arena, "Low memory");
}

/*  Pressure control  */
void
pressure_control(MDArena *arena)
{
	Vector cella_save, cellb_save, cellc_save, cello_save;
	Double de, total_energy_save, energies_save[kEndIndex];
	MDPressureArena *pressure = arena->pressure;
	Atom *ap;
	Transform tf;
	Vector v;
	Vector cello_new;
	Mat33 celltr_save;
	Double w, w0;
	int i;
	int needs_cell_recalculate = 0;
	XtalCell *cell = arena->mol->cell;
	
	if (pressure == NULL || (!arena->periodic_a && !arena->periodic_b && !arena->periodic_c))
		return;
	if (pressure->freq == 0 || arena->step % pressure->freq != 0)
		return;

	while (md_rand() < pressure->coupling) {

		Double flex;
		Transform tf0;
	
		/*  Trial move of the periodic cell  */
		memset(tf, 0, sizeof(tf));
		tf[0] = tf[4] = tf[8] = 1;
		memmove(tf0, tf, sizeof(tf));
		w0 = (md_rand() - 0.5) * pressure->trial_width;
		for (i = 0; i < 3; i++) {
			static char idx[] = {4, 8, 7, 5, 0, 8, 6, 2, 0, 4, 3, 1};
			flex = pressure->cell_flexibility[i + 3];
			if (flex > 0) {
				w = w0 = (md_rand() + md_rand() - 1.0) * pressure->trial_width * 0.5;
			} else if (flex < 0) {
				w = w0;
			} else w = 0;
			w *= fabs(flex);
			tf[idx[i * 4]] = tf[idx[i * 4 + 1]] = sqrt(1 - w * w);
			tf[idx[i * 4 + 2]] = tf[idx[i * 4 + 3]] = w;
			TransformMul(tf0, tf0, tf);
			tf[idx[i * 4]] = tf[idx[i * 4 + 1]] = 1;
			tf[idx[i * 4 + 2]] = tf[idx[i * 4 + 3]] = 0;
		}
		memmove(tf, tf0, sizeof(tf));
		w0 = (md_rand() - 0.5) * pressure->trial_width;
		for (i = 0; i < 3; i++) {
			int j;
			flex = pressure->cell_flexibility[i];
			if (flex > 0) {
				w = w0 = (md_rand() - 0.5) * pressure->trial_width;
			} else if (flex < 0) {
				w = w0;
			} else w = 0;
			for (j = 0; j < 3; j++) {
				tf[i * 3 + j] *= 1 + w * fabs(flex);
			}
		}
		
		/*  The new cell vector (bij) = (tik)(akj), so B = A*T  */
		/*  The cell transform matrix R = A*T*(A^-1)  */
		MatrixMul(tf, cell->tr, tf);
		MatrixMul(tf, tf, cell->rtr);
	/*	MatrixInvert(mat, arena->celltr);
		MatrixMul(tf, tf, mat); */

		cella_save = cell->axes[0];
		cellb_save = cell->axes[1];
		cellc_save = cell->axes[2];
		cello_save = cell->origin;
		memmove(celltr_save, cell->tr, sizeof(celltr_save));

		/*  Save a snapshot  */
		md_snapshot(arena, 0);
		total_energy_save = arena->total_energy;
		memmove(energies_save, arena->energies, sizeof(Double) * kEndIndex);

		/*  Transform the cell  */
		md_scale_cell(arena, tf, (pressure->use_atomic_move ? 1 : 2));

		/*  Fluctuate cell origin and orientation if requested  */
		if (pressure->cell_flexibility[6] != 0) {
			cello_new = cell->origin;
			cello_new.x += (md_rand() - 0.5) * pressure->fluctuate_cell_origin;
			cello_new.y += (md_rand() - 0.5) * pressure->fluctuate_cell_origin;
			cello_new.z += (md_rand() - 0.5) * pressure->fluctuate_cell_origin;
			cell->origin = cello_new;
			needs_cell_recalculate = 1;
		}
		if (pressure->cell_flexibility[7] != 0) {
			Vector axis;
			Mat33 mat2;
			switch ((int)pressure->cell_flexibility[7]) {
				case 1: axis.x = 1.0; axis.y = axis.z = 0.0; break;
				case 2: axis.y = 1.0; axis.x = axis.z = 0.0; break;
				case 3: axis.z = 1.0; axis.x = axis.y = 0.0; break;
				case 4:
					axis.x = md_rand() - 0.5;
					axis.y = md_rand() - 0.5;
					axis.z = md_rand() - 0.5;
					w = 1.0 / VecLength(axis);
					VecScaleSelf(axis, w);
					break;
				default:
					axis.x = axis.y = axis.z = 0.0; /* Just to keep compiler happy */
					break;
			}
			MatrixRotation(mat2, &axis, (md_rand() - 0.5) * pressure->fluctuate_cell_orientation);
			MatrixVec(&cell->axes[0], mat2, &cell->axes[0]);
			MatrixVec(&cell->axes[1], mat2, &cell->axes[1]);
			MatrixVec(&cell->axes[2], mat2, &cell->axes[2]);
			needs_cell_recalculate = 1;
		}
		if (needs_cell_recalculate)
			md_update_cell(arena);

		/*  Recalc the energies  */
		md_amend_by_symmetry(arena);
		calc_force(arena);
		md_calc_kinetic_energy(arena);
		pressure->mc_delta_potential = arena->total_energy - total_energy_save;
		pressure->mc_delta_kinetic = arena->energies[kKineticIndex] - energies_save[kKineticIndex];
		pressure->mc_stress = -s_lattice_deformation_stress(pressure->apply, cell->tr, celltr_save);
		de = pressure->mc_delta_potential + pressure->mc_delta_kinetic + pressure->mc_stress;
		pressure->mc_delta_volume = fabs(MatrixDeterminant(cell->tr)) - fabs(MatrixDeterminant(celltr_save));
		
		/*  Metropolis scheme  */
		if (de > 0 && md_rand() > exp(-de / (BOLTZMANN * arena->temperature))) {
			/*  Reject  */
			cell->axes[0] = cella_save;
			cell->axes[1] = cellb_save;
			cell->axes[2] = cellc_save;
			cell->origin = cello_save;
			md_update_cell(arena);
			if (pressure->control_algorithm != 1) {
				md_restore(arena, 0);
				arena->total_energy = total_energy_save;
				memmove(arena->energies, energies_save, sizeof(Double) * kEndIndex);
			}
			pressure->mc_reject++;
		} else {
			pressure->mc_accept++;
			if (pressure->control_algorithm == 1) {
				Vector dr;
				int j, k;
				/*  Update positions and velocities  */
				for (i = 0, ap = arena->mol->atoms; i < arena->mol->natoms; i++, ap++) {
					if (ap->periodic_exclude)
						continue;
					if (i < arena->natoms_uniq) {
						k = arena->fragment_indices[i];
						v = pressure->temporary_velocities[i];
					} else if (ap->symop.alive && (j = ap->symbase) >= 0 && j < arena->natoms_uniq) {
						k = arena->fragment_indices[j];
						v = pressure->temporary_velocities[j];
					} else continue;
					dr = arena->fragment_info[k].pos;
					if (ap->symop.alive) {
						md_transform_vec_by_symmetry(arena, &dr, &dr, ap->symop, 1);
						md_transform_vec_by_symmetry(arena, &v, &v, ap->symop, 1);
					}
					VecInc(ap->r, dr);
					VecInc(arena->verlets_dr[i], dr);
					ap->v = v;
				}
			}
		}
	/*	if (pressure->mc_callback != NULL) {
			if (Tcl_EvalObjEx((Tcl_Interp *)(arena->tcl_interp), pressure->mc_callback, 0) != TCL_OK) {
				arena->request_abort = 1;
				return;
			}
		} */
	} /* end while (mrand() < pressure->coupling) */
}

