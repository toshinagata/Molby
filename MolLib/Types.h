/*
 *  Types.h
 *
 *  Created by Toshi Nagata on 06/03/11.
 *  Copyright 2006-2008 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __Types_h__
#define __Types_h__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if defined(__WXMAC__) || defined(__CMDMAC__)
/*  On Mac OS X, CLAPACK is in Accelerate.framework  */
#include <vecLib/clapack.h>
#else
#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
	
#define STUB extern

#ifndef kRad2Deg
#define kRad2Deg  (180./3.14159265358979)
#endif
#ifndef kDeg2Rad
#define kDeg2Rad  (3.14159265358979 / 180.)
#endif
#ifndef kLog10
#define kLog10 2.3025851
#endif
#ifndef kArcCot15
#define kArcCot15 3.73205075
#endif

#undef kBohr2Angstrom
#define kBohr2Angstrom 0.529177249

#undef kAngstrom2Bohr
#define kAngstrom2Bohr 1.8897259885

typedef double Double;
typedef int Int;
typedef unsigned char Byte;
typedef unsigned int UInt;

#define kInvalidIndex -99999999  /*  Used for terminating integer array  */

typedef struct Vector { Double x, y, z; } Vector;
typedef struct Quat { Double x, y, z, w; } Quat;   /*  A quaternion  */
typedef Double Mat33[9];

/*  Mat33 (rot) and Vector (translation)  */
/*  Transform[0..8] are row-first, and Transform[9..11] are the fourth column. So that, the
    matrix elements are in the following order;
    {a11, a12, a13, a21, a22, a23, a31, a32, a33, a14, a24, a34}
   Sorry for confusion!  */
typedef Double Transform[12];
/*typedef Double Matrix[4][4]; */

#define VecAdd(v3, v1, v2) ((v3).x=(v1).x+(v2).x, (v3).y=(v1).y+(v2).y, (v3).z=(v1).z+(v2).z)
#define VecInc(v1, v2) ((v1).x+=(v2).x, (v1).y+=(v2).y, (v1).z+=(v2).z)
#define VecSub(v3, v1, v2) ((v3).x=(v1).x-(v2).x, (v3).y=(v1).y-(v2).y, (v3).z=(v1).z-(v2).z)
#define VecDec(v1, v2) ((v1).x-=(v2).x, (v1).y-=(v2).y, (v1).z-=(v2).z)
#define VecScale(v2, v1, r) ((v2).x=(v1).x*(r), (v2).y=(v1).y*(r), (v2).z=(v1).z*(r))
#define VecScaleSelf(v1, r) ((v1).x*=(r), (v1).y*=(r), (v1).z*=(r))
#define VecScaleInc(v2, v1, r) ((v2).x+=(v1).x*(r), (v2).y+=(v1).y*(r), (v2).z+=(v1).z*(r))
#define VecLength2(v1) ((v1).x*(v1).x+(v1).y*(v1).y+(v1).z*(v1).z)
#define VecLength(v1) (sqrt(VecLength2(v1)))
#define VecDot(v1, v2) ((v1).x*(v2).x+(v1).y*(v2).y+(v1).z*(v2).z)
#define VecCross(v3, v1, v2) ((v3).x=(v1).y*(v2).z-(v1).z*(v2).y, (v3).y=(v1).z*(v2).x-(v1).x*(v2).z, (v3).z=(v1).x*(v2).y-(v1).y*(v2).x)
#define VecZero(v1) ((v1).x=(v1).y=(v1).z=0.0)
#define VecIndex(vp, i) (((Double *)(vp))[i])
	
/*  Vector utility functions  */
int NormalizeVec(Vector *vdst, const Vector *vsrc);  /*  Returns non-zero when zero-vector */

void MatrixRotation(Mat33 dst, const Vector *axis, Double angle);
void MatrixVec(Vector *dst, const Mat33 mat, const Vector *src);
void MatrixMul(Mat33 dst, const Mat33 src1, const Mat33 src2);
void MatrixTranspose(Mat33 dst, const Mat33 src);
void MatrixScale(Mat33 dst, const Mat33 src, Double factor);
Double MatrixDeterminant(const Mat33 src);
int MatrixInvert(Mat33 dst, const Mat33 src);  /*  Return non-zero when determinant is zero */
int MatrixSymDiagonalize(Mat33 mat, Double *out_values, Vector *out_vectors);
void MatrixGeneralRotation(Mat33 dst, const Vector *v1, const Vector *v2, const Vector *v3);

void TransformVec(Vector *dst, const Transform tf, const Vector *src);
void TransformMul(Transform dst, const Transform src1, const Transform src2);
Double TransformDeterminant(const Transform src);
int TransformInvert(Transform dst, const Transform src);  /*  Return non-zero when determinant is zero */
void TransformTranspose(Transform dst, const Transform src);

void TransformForInversion(Transform dst, const Vector *center);
int  TransformForReflection(Transform dst, const Vector *axis, const Vector *center);
int  TransformForRotation(Transform dst, const Vector *axis, Double angle, const Vector *center);

/*  Utility functions  */
void SetPanicFunc(void (*func)(const char *, ...));
void SetWarningFunc(void (*func)(const char *, ...));
extern void (*gPanicFunc)(const char *, ...);
extern void (*gWarningFunc)(const char *, ...);
#define Panic (*gPanicFunc)
#define Warning (*gWarningFunc)

void PanicByOutOfMemory(const char *msg);
void PanicByInternalError(const char *msg, const char *file, int line);

#define MALLOC_CHECK(p, s) ((p) == NULL ? PanicByOutOfMemory(s) : (void)0)
#define NULL_CHECK(p, s) ((p) == NULL ? PanicByInternalError((s), __FILE__, __LINE__) : (void)0)

void *AssignArray(void *base, Int *count, int item_size, int idx, const void *value);
void *NewArray(void *base, Int *count, int item_size, int nitems);
void *InsertArray(void *base, Int *count, int item_size, int idx, int nitems, const void *value);
void *DeleteArray(void *base, Int *count, int item_size, int idx, int nitems, void *outValue);
int ReadLine(char *buf, int size, FILE *stream, int *lineNumber);
int ReadFormat(const char *str, const char *fmt,...);

#ifdef __cplusplus
}
#endif
		
#endif /* __Types_h__ */
