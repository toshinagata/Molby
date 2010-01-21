/*
 *  Dcd.h
 *  Molby
 *
 *  Created by Toshi Nagata on 09/01/20.
 *  Copyright 2009 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __DCD_H__
#define __DCD_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
		
typedef int Int32;
typedef float SFloat32;

typedef struct DcdRecord {
    int fd;                 /*  File descripter  */	
	int reverse_endian;     /*  Need to reverse endian?  */
    Int32 natoms;	        /*  Number of atoms  */
	Int32 nframes;	        /*  Number of frames  */
	Int32 nstart;           /*  Start timestep  */
	Int32 ninterval;        /*  Number of timesteps between frames  */
	Int32 nend;             /*  Last timestep  */
	Int32 with_unitcell;    /*  Has a unit cell information?  */
    SFloat32 delta;         /*  Step time  */
    SFloat32 globalcell[6]; /*  cell size and origin; used when with_unitcell == 0 */
    off_t header_size;      /*  Header size  */
} DcdRecord;

int DcdOpen(const char *name, DcdRecord *dr);
int DcdCreate(const char *name, DcdRecord *dr);
int DcdClose(DcdRecord *dr);
int DcdReadFrame(DcdRecord *dr, int index, SFloat32 *xp, SFloat32 *yp, SFloat32 *zp, SFloat32 *cellp);
int DcdWriteFrame(DcdRecord *dr, int index, const SFloat32 *xp, const SFloat32 *yp, const SFloat32 *zp, const SFloat32 *cellp);

#ifdef __cplusplus
}
#endif

#endif /* __DCD_H__ */
