/*
 *  Dcd.c
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

#include "Dcd.h"

static void
s_Swap4(char *cp)
{
	char w[4];
	w[0] = cp[0]; w[1] = cp[1]; w[2] = cp[2]; w[3] = cp[3];
	cp[0] = w[3]; cp[1] = w[2]; cp[2] = w[1]; cp[3] = w[0];
}

static void
s_Swap8(char *cp)
{
	char w[8];
	w[0] = cp[0]; w[1] = cp[1]; w[2] = cp[2]; w[3] = cp[3]; w[4] = cp[4]; w[5] = cp[5]; w[6] = cp[6]; w[7] = cp[7];
	cp[0] = w[7]; cp[1] = w[6]; cp[2] = w[5]; cp[3] = w[4]; cp[4] = w[3]; cp[5] = w[2]; cp[6] = w[1]; cp[7] = w[0];
}

#define s_SwapInt32(ip)  s_Swap4((char *)(ip))
#define s_SwapSFloat32(fp) s_Swap4((char *)(fp))

#define s_ReadInt32(dcd, ip) \
	(read(dcd->fd, ip, 4) == 4 ? \
		((dcd->reverse_endian ? s_SwapInt32(ip) : 0), 1) : \
		0)
#define s_ReadSFloat32(dcd, fp) \
	(read(dcd->fd, fp, 4) == 4 ? \
		((dcd->reverse_endian ? s_SwapSFloat32(fp) : 0), 1) : \
		0)

static int
s_Write4(int fd, const char *cp, int swap)
{
	if (swap) {
		char w[4];
		w[3] = cp[0]; w[2] = cp[1]; w[1] = cp[2]; w[0] = cp[3];
		return write(fd, w, 4);
	} else return write(fd, cp, 4);
}

static int
s_Write8(int fd, const char *cp, int swap)
{
	if (swap) {
		char w[8];
		w[7] = cp[0]; w[6] = cp[1]; w[5] = cp[2]; w[4] = cp[3]; w[3] = cp[4]; w[2] = cp[5]; w[1] = cp[6]; w[0] = cp[7];
		return write(fd, w, 8);
	} else return write(fd, cp, 8);
}

static int
s_WriteInt32(DcdRecord *dcd, Int32 i)
{
	return (s_Write4(dcd->fd, (const char *)(&i), dcd->reverse_endian) == 4);
}

static int
s_WriteSFloat32(DcdRecord *dcd, SFloat32 f)
{
	return (s_Write4(dcd->fd, (const char *)(&f), dcd->reverse_endian) == 4);
}

static int
s_WriteSFloat64(DcdRecord *dcd, SFloat64 d)
{
	return (s_Write8(dcd->fd, (const char *)(&d), dcd->reverse_endian) == 8);
}

static int
s_WriteZeros(int fd, int count)
{
	char buf[128];
	memset(buf, 0, sizeof(buf));
	while (count > 0) {
		int n = (count > 128 ? 128 : count);
		if (write(fd, buf, n) != n)
			return 0;
		count -= n;
	}
	return 1;
}

int
DcdOpen(const char *name, DcdRecord *dr)
{
    Int32 nn, nnn;
    char buf[5];
	char delta_buf[8];

	if (dr == NULL)
		return -1;  /*  Internal error  */
	memset(dr, 0, sizeof(DcdRecord));
	dr->fd = open(name, O_RDONLY, 0);
	if (dr->fd < 0)
		return -2;  /*  Cannot open file  */
	
	/*  Section 1: 'CORD', NFILE, NPRIV, NSAVC, NSTEP, 4*ZERO, NFIX,  */
	/* DELTA (4 bytes for charmm format, 8 bytes for x-plor format),  */
	/* NEXTRA (charmm only), N4DIM, 6*ZERO, NCHARM, 24 */
	if (!s_ReadInt32(dr, &nn))
		return 1;   /*  Bad format: premature EOF  */
	if (nn != 84) {
		s_SwapInt32(&nn);
		if (nn == 84)
			dr->reverse_endian = 1;
		else return 2;  /*  Bad format: bad block length of the first section  */
	}
    if (read(dr->fd, buf, 4) != 4 || strncmp(buf, "CORD", 4) != 0)
        return 3;   /*  Bad format: missing "CORD" signature  */
    s_ReadInt32(dr, &(dr->nframes));
	s_ReadInt32(dr, &(dr->nstart));
	s_ReadInt32(dr, &(dr->ninterval));
	s_ReadInt32(dr, &(dr->nend));
    lseek(dr->fd, 16, SEEK_CUR);  /*  Skip 4 zeros  */
	s_ReadInt32(dr, &(dr->nfix));
	read(dr->fd, delta_buf, 8);   /*  If charmm format, then this is { float; int32; } otherwise this is a double */
	s_ReadInt32(dr, &(dr->n4dim));
	lseek(dr->fd, 28, SEEK_CUR);  /*  Skip 7 zeros  */
	s_ReadInt32(dr, &(dr->ncharmver));
	s_ReadInt32(dr, &nn);  /*  This should be 84  */
	if (nn != 84)
		return 4;
	if (dr->ncharmver == 0) {
		if (dr->reverse_endian)
			s_Swap8(delta_buf);
		dr->delta = *((SFloat64 *)delta_buf);
		dr->nextra = 0;
	} else {
		if (dr->reverse_endian) {
			s_Swap4(delta_buf);
			s_Swap4(delta_buf + 4);
		}
		dr->delta = *((SFloat32 *)delta_buf);
		dr->nextra = *((Int32 *)(delta_buf + 4));
	}
	
    /*  Section 2: Title lines  */
    if (!s_ReadInt32(dr, &nn) || lseek(dr->fd, nn, SEEK_CUR) < 0 || !s_ReadInt32(dr, &nnn) || nn != nnn)
		return 5;   /*  Bad format: the second section looks strange  */
    
    /*  Section 3: Number of atoms  */
	if (!s_ReadInt32(dr, &nn) || nn != 4)
		return 6;   /*  Bad format: the third section is not correct  */
	s_ReadInt32(dr, &(dr->natoms));
	if (!s_ReadInt32(dr, &nn) || nn != 4)
		return 6;

	dr->header_size = lseek(dr->fd, 0, SEEK_CUR);
	return 0;
}

int
DcdCreate(const char *name, DcdRecord *dr)
{
	char buf[80];
	
	if (dr == NULL)
		return -1;  /*  Internal error  */
	dr->fd = open(name, O_WRONLY|O_CREAT|O_TRUNC, 0666);
	if (dr->fd < 0)
		return -2;  /*  Cannot create file  */
	memset(buf, ' ', 80);
	
	/*  Section 1: 'CORD', NFILE, NPRIV, NSAVC, NSTEP, 4*ZERO, NFIX,  */
	/* DELTA (4 bytes for charmm format, 8 bytes for x-plor format),  */
	/* NEXTRA (charmm only), N4DIM, 6*ZERO, NCHARM, 24 */
	if (!s_WriteInt32(dr, 84) ||
		write(dr->fd, "CORD", 4) != 4 ||
		!s_WriteInt32(dr, dr->nframes) ||
		!s_WriteInt32(dr, dr->nstart) ||
		!s_WriteInt32(dr, dr->ninterval) || 
		!s_WriteInt32(dr, dr->nend) ||
		!s_WriteZeros(dr->fd, 16) ||
		!s_WriteInt32(dr, dr->nfix) ||
		!s_WriteSFloat32(dr, dr->delta) ||
		!s_WriteInt32(dr, dr->nextra) ||
		!s_WriteZeros(dr->fd, 32) ||
		!s_WriteInt32(dr, dr->ncharmver) ||
		!s_WriteInt32(dr, 84))
		return 1;   /*  Cannot write  */
	
    /*  Section 2: Title lines (92 bytes)  */
	if (!s_WriteInt32(dr, 84) ||
		!s_WriteInt32(dr, 1) ||
		write(dr->fd, buf, 80) != 80 ||
		!s_WriteInt32(dr, 84))
		return 1;   /*  Cannot write  */
    
    /*  Section 3: Number of atoms (12 bytes)  */
	if (!s_WriteInt32(dr, 4) ||
		!s_WriteInt32(dr, dr->natoms) ||
		!s_WriteInt32(dr, 4))
		return 1;
	
	/*  The header size is 196 bytes  */
	dr->header_size = lseek(dr->fd, 0, SEEK_CUR);
    return 0;
}

int
DcdClose(DcdRecord *dr)
{
	if (dr->fd >= 0)
		close(dr->fd);
	dr->fd = -1;
	return 0;
}

/*  Read one frame. Xp, yp, zp should be larger enough for sizeof(SFloat32)*dr->natoms.
    Cellp may be NULL, but if non-NULL then it should be large enough for sizeof(SFloat32)*6.  */
int
DcdReadFrame(DcdRecord *dr, int index, SFloat32 *xp, SFloat32 *yp, SFloat32 *zp, SFloat32 *cellp)
{
    Int32 nn, nnn, i;
    off_t block_size = (index > 0 ? dr->block_size : 24 + (off_t)(dr->natoms * 12) + (dr->nextra ? 56 : 0));
    lseek(dr->fd, dr->header_size + block_size * index, SEEK_SET);
	block_size = 0;
	if (!s_ReadInt32(dr, &nn))
		goto error;
	if (nn == 48) {
		if (dr->nextra) {
			SFloat64 mycell[6];
			if (nn != 48)
				goto error;
			read(dr->fd, mycell, 48);
			if (cellp != NULL) {
				for (i = 0; i < 6; i++) {
					if (dr->reverse_endian)
						s_Swap8((char *)(mycell + i));
					cellp[i] = mycell[i];  /*  double -> float  */
				}
			}
			if (!s_ReadInt32(dr, &nn) || nn != 48)
				goto error;
		} else {
			lseek(dr->fd, 52, SEEK_CUR);
		}
		if (!s_ReadInt32(dr, &nn))
			goto error;
		block_size = 56;
	}
    if (nn != dr->natoms * 4 ||
		read(dr->fd, xp, nn) < nn ||
		!s_ReadInt32(dr, &nnn) || nn != nnn)
        goto error;
	if (dr->reverse_endian) {
		for (i = 0; i < dr->natoms; i++)
			s_SwapSFloat32(xp + i);
	}
    if (!s_ReadInt32(dr, &nn) || nn != dr->natoms * 4 || 
		read(dr->fd, yp, nn) < nn || 
		!s_ReadInt32(dr, &nnn) || nn != nnn)
        goto error;
	if (dr->reverse_endian) {
		for (i = 0; i < dr->natoms; i++)
			s_SwapSFloat32(yp + i);
	}
    if (!s_ReadInt32(dr, &nn) || nn != dr->natoms * 4 || 
		read(dr->fd, zp, nn) < nn || 
		!s_ReadInt32(dr, &nnn) || nn != nnn)
        goto error;
	if (dr->reverse_endian) {
		for (i = 0; i < dr->natoms; i++)
			s_SwapSFloat32(zp + i);
	}
	block_size += 24 + (off_t)(dr->natoms * 12);
	if (index == 0)
		dr->block_size = block_size;
    return 0;
error:
    return 1;
}

/*  Write one frame. Xp, yp, zp should be an array of SFloat32 * dr->natoms. Cellp may be
    NULL, but if non-NULL it should be an array of SFloat32 * 6.
    This function does _not_ do lseek, so it should be called sequentially from index 0.  */
int
DcdWriteFrame(DcdRecord *dr, int index, const SFloat32 *xp, const SFloat32 *yp, const SFloat32 *zp, const SFloat32 *cellp)
{
	Int32 i;
    if (dr->nextra) {
        if (!s_WriteInt32(dr, 48))
            goto error;
		for (i = 0; i < 6; i++)
			s_WriteSFloat64(dr, (cellp != NULL ? cellp[i] : 0.0));
		if (!s_WriteInt32(dr, 48))
            goto error;
    }
	
	/*  X coordinates  */
	if (!s_WriteInt32(dr, dr->natoms * 4))
		goto error;
	for (i = 0; i < dr->natoms; i++) {
		if (!s_WriteSFloat32(dr, xp[i]))
			goto error;
	}
	if (!s_WriteInt32(dr, dr->natoms * 4))
		goto error;

	/*  Y coordinates  */
	if (!s_WriteInt32(dr, dr->natoms * 4))
		goto error;
	for (i = 0; i < dr->natoms; i++) {
		if (!s_WriteSFloat32(dr, yp[i]))
			goto error;
	}
	if (!s_WriteInt32(dr, dr->natoms * 4))
		goto error;
	
	/*  Z coordinates  */
	if (!s_WriteInt32(dr, dr->natoms * 4))
		goto error;
	for (i = 0; i < dr->natoms; i++) {
		if (!s_WriteSFloat32(dr, zp[i]))
			goto error;
	}
	if (!s_WriteInt32(dr, dr->natoms * 4))
		goto error;

	return 0;

error:
	return 1;
}
