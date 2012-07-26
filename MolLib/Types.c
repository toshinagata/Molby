/*
 *  Types.c
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

#include "Types.h"
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

static void
defaultPanic(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	abort();
}

static void
defaultWarning(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
}

void (*gPanicFunc)(const char *, ...) = defaultPanic;
void (*gWarningFunc)(const char *, ...) = defaultWarning;

void
SetPanicFunc(void (*func)(const char *, ...))
{
	if (func == NULL)
		gPanicFunc = defaultPanic;
	else gPanicFunc = func;
}

void
SetWarningFunc(void (*func)(const char *, ...))
{
	if (func == NULL)
		gWarningFunc = defaultWarning;
	else gWarningFunc = func;
}

void
PanicByOutOfMemory(const char *msg)
{
	Panic("Out of memory while %s", msg);
}

void
PanicByInternalError(const char *msg, const char *file, int line)
{
	Panic("Internal Error in %s; File %s line %d", msg, file, line);
}


#pragma mark ==== Vector and Matrix ====

int
NormalizeVec(Vector *vdst, const Vector *vsrc)
{
	double dval;
	dval = 1.0 / VecLength(*vsrc);
	VecScale(*vdst, *vsrc, dval);
	if (!isfinite(vdst->x) || !isfinite(vdst->y) || !isfinite(vdst->z)) {
		return 1;
	}
	return 0;
}

void
MatrixRotation(Mat33 dst, const Vector *axis, Double angle)
{
	/*  Taken from "Matrix and Quaternion FAQ"
		http://www.j3d.org/matrix_faq/matrfaq_latest.html  */
	double rcos = cos((double)angle);
	double rsin = sin((double)angle);
	dst[0] =            rcos + axis->x*axis->x*(1-rcos);
	dst[1] = -axis->z * rsin + axis->x*axis->y*(1-rcos);
	dst[2] =  axis->y * rsin + axis->x*axis->z*(1-rcos);
	dst[3] =  axis->z * rsin + axis->y*axis->x*(1-rcos);
	dst[4] =            rcos + axis->y*axis->y*(1-rcos);
	dst[5] = -axis->x * rsin + axis->y*axis->z*(1-rcos);
	dst[6] = -axis->y * rsin + axis->z*axis->x*(1-rcos);
	dst[7] =  axis->x * rsin + axis->z*axis->y*(1-rcos);
	dst[8] =            rcos + axis->z*axis->z*(1-rcos);
}

void
MatrixVec(Vector *dst, const Mat33 mat, const Vector *src)
{
	Vector temp;
	temp.x = mat[0] * src->x + mat[3] * src->y + mat[6] * src->z;
	temp.y = mat[1] * src->x + mat[4] * src->y + mat[7] * src->z;
	temp.z = mat[2] * src->x + mat[5] * src->y + mat[8] * src->z;
	*dst = temp;
}

void
MatrixMul(Mat33 dst, const Mat33 src1, const Mat33 src2)
{
	Mat33 temp;
	temp[0] = src1[0] * src2[0] + src1[3] * src2[1] + src1[6] * src2[2];
	temp[1] = src1[1] * src2[0] + src1[4] * src2[1] + src1[7] * src2[2];
	temp[2] = src1[2] * src2[0] + src1[5] * src2[1] + src1[8] * src2[2];
	temp[3] = src1[0] * src2[3] + src1[3] * src2[4] + src1[6] * src2[5];
	temp[4] = src1[1] * src2[3] + src1[4] * src2[4] + src1[7] * src2[5];
	temp[5] = src1[2] * src2[3] + src1[5] * src2[4] + src1[8] * src2[5];
	temp[6] = src1[0] * src2[6] + src1[3] * src2[7] + src1[6] * src2[8];
	temp[7] = src1[1] * src2[6] + src1[4] * src2[7] + src1[7] * src2[8];
	temp[8] = src1[2] * src2[6] + src1[5] * src2[7] + src1[8] * src2[8];
	memmove(dst, temp, sizeof(Mat33));
}

void
MatrixTranspose(Mat33 dst, const Mat33 src)
{
	Double w;
	dst[0] = src[0];
	dst[4] = src[4];
	dst[8] = src[8];
	w = src[3]; dst[3] = src[1]; dst[1] = w;
	w = src[6]; dst[6] = src[2]; dst[2] = w;
	w = src[7]; dst[7] = src[5]; dst[5] = w;
}

Double
MatrixDeterminant(const Mat33 src)
{
	return src[0] * src[4] * src[8] + src[1] * src[5] * src[6] + src[2] * src[3] * src[7]
		- src[0] * src[5] * src[7] - src[1] * src[3] * src[8] - src[2] * src[4] * src[6];
}

int
MatrixInvert(Mat33 dst, const Mat33 src)
{
	Mat33 temp;
	Double d = MatrixDeterminant(src);
	if (d == 0.0)
		return 1;
	d = 1.0 / d;
	temp[0] = d * ( src[4] * src[8] - src[5] * src[7]);
	temp[1] = d * (-src[1] * src[8] + src[2] * src[7]);
	temp[2] = d * ( src[1] * src[5] - src[2] * src[4]);
	temp[3] = d * (-src[3] * src[8] + src[5] * src[6]);
	temp[4] = d * ( src[0] * src[8] - src[2] * src[6]);
	temp[5] = d * (-src[0] * src[5] + src[2] * src[3]);
	temp[6] = d * ( src[3] * src[7] - src[4] * src[6]);
	temp[7] = d * (-src[0] * src[7] + src[1] * src[6]);
	temp[8] = d * ( src[0] * src[4] - src[1] * src[3]);
	memmove(dst, temp, sizeof(Mat33));
	return 0;
}

void
MatrixScale(Mat33 dst, const Mat33 src, Double factor)
{
	dst[0] = src[0] * factor;
	dst[1] = src[1] * factor;
	dst[2] = src[2] * factor;
	dst[3] = src[3] * factor;
	dst[4] = src[4] * factor;
	dst[5] = src[5] * factor;
	dst[6] = src[6] * factor;
	dst[7] = src[7] * factor;
	dst[8] = src[8] * factor;
}

/*  Get the matrix to rotate (1,0,0)->v1, (0,1,0)->v2, (0,0,1)->v3  */
void
MatrixGeneralRotation(Mat33 dst, const Vector *v1, const Vector *v2, const Vector *v3)
{
	dst[0] = v1->x; dst[3] = v2->x; dst[6] = v3->x;
	dst[1] = v1->y; dst[4] = v2->y; dst[7] = v3->y;
	dst[2] = v1->z; dst[5] = v2->z; dst[8] = v3->z;
}

int
MatrixSymDiagonalize(Mat33 mat, Double *out_values, Vector *out_vectors)
{
	__CLPK_integer n, lda, lwork, info;
	__CLPK_doublereal a[9], w[3], work[9];
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			a[i * 3 + j] = mat[i * 3 + j];
		}
	}
	n = lda = 3;
	lwork = 9;
	/*  For the meanings of the arguments, consult the LAPACK source; 
		http://www.netlib.org/lapack/double/dsyev.f */
	dsyev_("V", "U", &n, a, &lda, w, work, &lwork, &info);
	if (info == 0) {
		for (i = 0; i < 3; i++) {
			out_values[i] = w[i];
			out_vectors[i].x = a[i * 3];
			out_vectors[i].y = a[i * 3 + 1];
			out_vectors[i].z = a[i * 3 + 2];
		}
	}
	return info;
}

void
TransformVec(Vector *dst, const Transform tf, const Vector *src)
{
	Vector temp;
	temp.x = tf[0] * src->x + tf[3] * src->y + tf[6] * src->z + tf[9];
	temp.y = tf[1] * src->x + tf[4] * src->y + tf[7] * src->z + tf[10];
	temp.z = tf[2] * src->x + tf[5] * src->y + tf[8] * src->z + tf[11];
	*dst = temp;
}

void
TransformMul(Transform dst, const Transform src1, const Transform src2)
{
	Transform temp;
	temp[0] = src1[0] * src2[0] + src1[3] * src2[1] + src1[6] * src2[2];
	temp[1] = src1[1] * src2[0] + src1[4] * src2[1] + src1[7] * src2[2];
	temp[2] = src1[2] * src2[0] + src1[5] * src2[1] + src1[8] * src2[2];
	temp[3] = src1[0] * src2[3] + src1[3] * src2[4] + src1[6] * src2[5];
	temp[4] = src1[1] * src2[3] + src1[4] * src2[4] + src1[7] * src2[5];
	temp[5] = src1[2] * src2[3] + src1[5] * src2[4] + src1[8] * src2[5];
	temp[6] = src1[0] * src2[6] + src1[3] * src2[7] + src1[6] * src2[8];
	temp[7] = src1[1] * src2[6] + src1[4] * src2[7] + src1[7] * src2[8];
	temp[8] = src1[2] * src2[6] + src1[5] * src2[7] + src1[8] * src2[8];
	temp[9] = src1[0] * src2[9] + src1[3] * src2[10] + src1[6] * src2[11] + src1[9];
	temp[10] = src1[1] * src2[9] + src1[4] * src2[10] + src1[7] * src2[11] + src1[10];
	temp[11] = src1[2] * src2[9] + src1[5] * src2[10] + src1[8] * src2[11] + src1[11];
	memmove(dst, temp, sizeof(Transform));
}

Double
TransformDeterminant(const Transform src)
{
	return MatrixDeterminant(src);
}

void
TransformTranspose(Transform dst, const Transform src)
{
	Double w;
	dst[0] = src[0];
	dst[4] = src[4];
	dst[8] = src[8];
	w = src[3]; dst[3] = src[1]; dst[1] = w;
	w = src[6]; dst[6] = src[2]; dst[2] = w;
	w = src[7]; dst[7] = src[5]; dst[5] = w;
	dst[9] = src[9];
	dst[10] = src[10];
	dst[11] = src[11];
}

int
TransformInvert(Transform dst, const Transform src)
{
	Transform temp;
	int n = MatrixInvert(temp, src);
	if (n == 0) {
		temp[9] = -temp[0] * src[9] - temp[3] * src[10] - temp[6] * src[11];
		temp[10] = -temp[1] * src[9] - temp[4] * src[10] - temp[7] * src[11];
		temp[11] = -temp[2] * src[9] - temp[5] * src[10] - temp[8] * src[11];
		memmove(dst, temp, sizeof(Transform));
		return 0;
	} else return n;
}

void
TransformForInversion(Transform dst, const Vector *center)
{
	dst[0] = dst[4] = dst[8] = -1.0;
	dst[1] = dst[2] = dst[3] = dst[5] = dst[6] = dst[7] = 0.0;
	dst[9] = center->x * 2.0;
	dst[10] = center->y * 2.0;
	dst[11] = center->z * 2.0;
}

int
TransformForReflection(Transform dst, const Vector *axis, const Vector *center)
{
	Vector av = *axis;
	Double w;
	if ((w = VecLength2(av)) < 1e-15)
		return 1;
	w = 1.0 / sqrt(w);
	VecScaleSelf(av, w);
	/*  r' = r - 2 * VecDot(r-c, a) * a; a should be a unit vector */
	/*  (x',y',z') = (x,y,z) - 2*((x-cx)*ax + (y-cy)*ay + (z-cz)*az)*(ax,ay,az) */
	/*  = (1-2*ax*ax, -2*ax*ay, -2*ax*az)x + (-2*ax*ay, 1-2*ay*ay, -2*ay*az)y
	 + (-2*ax*az, -2*ay*az, 1-2*az*az)z + (ax*C, ay*C, az*C); C = 2*(ax*cx+ay*cy+az*cz)  */
	dst[0] = 1.0 - 2.0 * av.x * av.x;
	dst[1] = dst[3] = -2.0 * av.x * av.y;
	dst[2] = dst[6] = -2.0 * av.x * av.z;
	dst[4] = 1.0 - 2.0 * av.y * av.y;
	dst[5] = dst[7] = -2.0 * av.y * av.z;
	dst[8] = 1.0 - 2.0 * av.z * av.z;
	w = 2.0 * VecDot(av, *center);
	dst[9] = av.x * w;
	dst[10] = av.y * w;
	dst[11] = av.z * w;
	return 0;
}

int
TransformForRotation(Transform dst, const Vector *axis, Double angle, const Vector *center)
{
	Transform tf, temp1;
	Double w;
	Vector av = *axis;
	if ((w = VecLength2(av)) < 1e-15)
		return 1;
	w = 1.0 / sqrt(w);
	VecScaleSelf(av, w);
	memset(tf, 0, sizeof(tf));
	memset(temp1, 0, sizeof(tf));
	tf[0] = tf[4] = tf[8] = 1.0;
	tf[1] = tf[2] = tf[3] = tf[5] = tf[6] = tf[7] = 0.0;
	tf[9] = -center->x;
	tf[10] = -center->y;
	tf[11] = -center->z;
	MatrixRotation((Double *)temp1, &av, angle);
	temp1[9] = center->x;
	temp1[10] = center->y;
	temp1[11] = center->z;
	TransformMul(dst, temp1, tf);
	return 0;
}

Transform gIdentityTransform = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

#pragma mark ==== LAMatrix ====

typedef struct LAMatrixTempRecord {
	int size;
	LAMatrix *mat;
} LAMatrixTempRecord;

static Int sNTempRecords;
static LAMatrixTempRecord *sTempRecords = NULL;

LAMatrix *
LAMatrixAllocTempMatrix(int row, int column)
{
	int i, n;
	LAMatrixTempRecord *tp;
	/*  Look for already allocated records  */
	n = -1;
	for (i = 0, tp = sTempRecords; i < sNTempRecords; i++, tp++) {
		if (tp->size >= column * row) {
			tp->size = -tp->size;
			tp->mat->column = column;
			tp->mat->row = row;
			memset(tp->mat->data, 0, sizeof(__CLPK_doublereal) * column * row);
			return tp->mat;
		}
		if (tp->size > 0)
			n = i;  /*  Record the unused entry  */
	}
	if (n == -1) {
		tp = (LAMatrixTempRecord *)AssignArray(&sTempRecords, &sNTempRecords, sizeof(LAMatrixTempRecord), sNTempRecords, NULL);
		tp->mat = NULL;
	} else tp = &sTempRecords[n];
	tp->mat = (LAMatrix *)realloc(tp->mat, sizeof(LAMatrix) + sizeof(__CLPK_doublereal) * (column * row - 1));
	tp->size = -column * row;
	tp->mat->column = column;
	tp->mat->row = row;
	memset(tp->mat->data, 0, sizeof(__CLPK_doublereal) * column * row);
	return tp->mat;
}

void
LAMatrixReleaseTempMatrix(LAMatrix *mat)
{
	int i;
	LAMatrixTempRecord *tp;
	/*  Is this record found in sTempRecords?  */
	for (i = 0, tp = sTempRecords; i < sNTempRecords; i++, tp++) {
		if (tp->mat == mat) {
			tp->size = -tp->size;
			return;
		}
	}
	/*  If not, then just free it  */
	free(mat);
}

LAMatrix *
LAMatrixNew(int row, int column)
{
	LAMatrix *m = (LAMatrix *)calloc(sizeof(LAMatrix) + sizeof(__CLPK_doublereal) * (column * row - 1), 1);
	m->column = column;
	m->row = row;
	return m;
}

LAMatrix *
LAMatrixResize(LAMatrix *mat, int row, int column)
{
	if (mat == NULL || mat->row * mat->column < row * column)
		mat = (LAMatrix *)realloc(mat, sizeof(LAMatrix) + sizeof(__CLPK_doublereal) * (column * row - 1));
	mat->column = column;
	mat->row = row;
	return mat;
}

void
LAMatrixRelease(LAMatrix *mat)
{
	if (mat != NULL)
		free(mat);
}

LAMatrix *
LAMatrixNewFromMatrix(const LAMatrix *mat)
{
	LAMatrix *m = LAMatrixNew(mat->row, mat->column);
	memmove(m->data, mat->data, sizeof(m->data[0]) * mat->row * mat->column);
	return m;
}

void
LAMatrixMul(int trans1, int trans2, double scale1, const LAMatrix *mat1, const LAMatrix *mat2, double scale2, LAMatrix *mat3)
{
#if defined(__WXMAC__) || defined(__CMDMAC__)
	int m, n, k;
	if (trans1) {
		trans1 = CblasTrans;
		m = mat1->column;
		k = mat1->row;
	} else {
		trans1 = CblasNoTrans;
		m = mat1->row;
		k = mat1->column;
	}
	if (trans2) {
		trans2 = CblasTrans;
		n = mat2->row;
	} else {
		trans2 = CblasNoTrans;
		n = mat2->column;
	}
	cblas_dgemm(CblasColMajor, trans1, trans2, m, n, k, scale1, mat1->data, mat1->row, mat2->data, mat2->row, scale2, mat3->data, mat3->row);
#else
	char ctrans1, ctrans2;
	int m, n, k;
	if (trans1) {
		ctrans1 = 'T';
		m = mat1->column;
		k = mat1->row;
	} else {
		ctrans1 = 'N';
		m = mat1->row;
		k = mat1->column;
	}
	if (trans2) {
		ctrans2 = 'T';
		n = mat2->row;
	} else {
		ctrans2 = 'N';
		n = mat2->column;
	}
	dgemm_(&ctrans1, &ctrans2, &m, &n, &k, &scale1, mat1->data, &mat1->row, mat2->data, &mat2->row, &scale2, mat3->data, &mat3->row);
#endif
}

int
LAMatrixInvert(LAMatrix *mat1, const LAMatrix *mat2)
{
	LAMatrix *tmat1, *tmat2;
	__CLPK_integer m, n, lda, info, lwork;
	if (mat2->column != mat2->row || mat1->column * mat1->row < mat2->column * mat2->row)
		return -1;  /*  Wrong dimension  */
	lwork = m = n = lda = mat1->column;
	tmat1 = LAMatrixAllocTempMatrix(n, n);  /*  For work  */
	tmat2 = LAMatrixAllocTempMatrix(n, 1);  /*  For piv   */
	if (mat1 != mat2)
		memmove(mat1->data, mat2->data, sizeof(__CLPK_doublereal) * n * n);
	dgetrf_(&m, &n, mat1->data, &lda, (__CLPK_integer *)tmat2->data, &info);
	if (info == 0)
		dgetri_(&n, mat1->data, &lda, (__CLPK_integer *)tmat2->data, tmat1->data, &lwork, &info);
	LAMatrixReleaseTempMatrix(tmat1);
	LAMatrixReleaseTempMatrix(tmat2);
	mat1->column = mat1->row = m;
	return info;
}

Double
LAMatrixDeterminant(const LAMatrix *mat)
{
	LAMatrix *tmat1, *tmat2;
	Double det = 1.0;
	__CLPK_integer m, n, mn, i, lda, info;
	lda = m = mat->row;
	n = mat->column;
	if (m < n)
		mn = m;
	else mn = n;
	tmat1 = LAMatrixAllocTempMatrix(m, n);
	tmat2 = LAMatrixAllocTempMatrix(mn, 1);  /* For piv */
	memmove(tmat1->data, mat->data, sizeof(__CLPK_doublereal) * m * n);
	dgetrf_(&m, &n, tmat1->data, &lda, (__CLPK_integer *)tmat2->data, &info);
	if (info == 0) {
		for (i = 0; i < mn - 1; i++) {
			if (((__CLPK_integer *)tmat2->data)[i] != i + 1)
				det = -det;
		}
		for (i = 0; i < mn; i++)
			det *= tmat1->data[(m + 1) * i];
	} else det = 0.0;
	LAMatrixReleaseTempMatrix(tmat1);
	LAMatrixReleaseTempMatrix(tmat2);
	return det;
}

void
LAMatrixTranspose(LAMatrix *mat1, const LAMatrix *mat2)
{
	LAMatrix *mat3;
	int i, j;
	if (mat1 == mat2) {
		mat3 = LAMatrixAllocTempMatrix(mat2->row, mat2->column);
		memmove(mat3->data, mat2->data, sizeof(__CLPK_doublereal) * mat2->column * mat2->row);
	} else mat3 = (LAMatrix *)mat2;
	if (mat1->column * mat1->row >= mat3->column * mat3->row) {
		mat1->column = mat3->row;
		mat1->row = mat3->column;
		for (i = 0; i < mat1->column; i++) {
			for (j = 0; j < mat1->row; j++) {
				mat1->data[i * mat1->row + j] = mat3->data[j * mat1->column + i];
			}
		}
	}
	if (mat1 == mat2)
		LAMatrixReleaseTempMatrix(mat3);
}

/*  Diagonalize the symmetric matrix  */
/*  eigenValues = (m, 1) Matrix, eigenVectors = (m, m) Matrix (column vectors are the eigenvectors),
    mat2 = (m, m) symmetric Matrix  */
int
LAMatrixSymDiagonalize(LAMatrix *eigenValues, LAMatrix *eigenVectors, const LAMatrix *mat)
{
	__CLPK_integer n, lda, lwork, info;
	__CLPK_doublereal dwork;
	LAMatrix *tmat1;
	if (mat->column != mat->row || eigenVectors->column * eigenVectors->row < mat->column * mat->row || eigenValues->column * eigenValues->row < mat->column)
		return -1;  /*  Illegal dimension  */
	n = lda = mat->column;
	memmove(eigenVectors->data, mat->data, sizeof(__CLPK_doublereal) * n * n);
	lwork = -1;  /*  workspace query  */
	dsyev_("V", "U", &n, eigenVectors->data, &lda, eigenValues->data, &dwork, &lwork, &info);
	if (info == 0) {
		lwork = dwork;
		tmat1 = LAMatrixAllocTempMatrix(lwork, 1);
		dsyev_("V", "U", &n, eigenVectors->data, &lda, eigenValues->data, tmat1->data, &lwork, &info);
		LAMatrixReleaseTempMatrix(tmat1);
	}
	eigenValues->row = n;
	eigenValues->column = 1;
	eigenVectors->row = eigenVectors->column = n;
	return info;
}

#pragma mark ==== Array ====

/*  Assign a value to an array. An array is represented by two fields; count and base,
 *  where base is a pointer to an array and count is the number of items.
 *  The memory block of the array is allocated by 8*item_size. If the index exceeds
 *  that limit, then a new memory block is allocated.  */
void *
AssignArray(void *base, Int *count, int item_size, int idx, const void *value)
{
	void **bp = (void **)base;
	if (*count == 0 || idx / 8 > (*count - 1) / 8) {
		int new_size = (idx / 8 + 1) * 8;
		if (*bp == NULL)
			*bp = calloc(item_size, new_size);
		else
			*bp = realloc(*bp, new_size * item_size);
		if (*bp == NULL)
			return NULL;
		memset((char *)*bp + *count * item_size, 0, (new_size - *count) * item_size);
	}
	if (idx >= *count)
		*count = idx + 1;
	if (value != NULL)
		memcpy((char *)*bp + idx * item_size, value, item_size);
	return (char *)*bp + idx * item_size;
}

/*  Allocate a new array. This works consistently with AssignArray().
 *  Don't mix calloc()/malloc() with AssignArray(); that causes disasters!
 *  (free() is OK though).  */
void *
NewArray(void *base, Int *count, int item_size, int nitems)
{
	void **bp = (void *)base;
	*bp = NULL;
	*count = 0;
	if (nitems > 0)
		return AssignArray(base, count, item_size, nitems - 1, NULL);
	else return NULL;
}

/*  Insert items to an array.  */
void *
InsertArray(void *base, Int *count, int item_size, int idx, int nitems, const void *value)
{
	void **bp = (void *)base;
	void *p;
	int ocount = *count;
	if (nitems <= 0)
		return NULL;
	/*  Allocate storage  */
	p = AssignArray(base, count, item_size, *count + nitems - 1, NULL);
	if (p == NULL)
		return NULL;
	/*  Move items if necessary  */
	if (idx < ocount)
		memmove((char *)*bp + (idx + nitems) * item_size, (char *)*bp + idx * item_size, (ocount - idx) * item_size);
	/*  Copy items  */
	if (value != NULL)
		memmove((char *)*bp + idx * item_size, value, nitems * item_size);
	else
		memset((char *)*bp + idx * item_size, 0, nitems * item_size);
	return (char *)*bp + idx * item_size;
}

void *
DeleteArray(void *base, Int *count, int item_size, int idx, int nitems, void *outValue)
{
	void **bp = (void *)base;
	if (nitems <= 0 || idx < 0 || idx >= *count)
		return NULL;
	if (nitems > *count - idx)
		nitems = *count - idx;
	/*  Copy items  */
	if (outValue != NULL)
		memmove(outValue, (char *)*bp + idx * item_size, nitems * item_size);
	/*  Move items  */
	if (idx + nitems < *count)
		memmove((char *)*bp + idx * item_size, (char *)*bp + (idx + nitems) * item_size, (*count - idx - nitems) * item_size);
	*count -= nitems;
	if (*count == 0) {
		free(*bp);
		*bp = NULL;
	}
	return NULL;
}

int
ReadLine(char *buf, int size, FILE *stream, int *lineNumber)
{
	int i, c;
	i = 0;
	c = 0;
	while (i < size - 1) {
		c = getc(stream);
		if (c == EOF)
			break;
		buf[i++] = c;
		if (c == '\n')
			break;
		else if (c == '\r') {
			c = getc(stream);
			if (c != '\n')
				ungetc(c, stream);
			c = '\n';
			break;
		}
	}
	buf[i] = 0;
	if (c != '\n' && c != EOF) {
		/*  Skip until the end of line  */
		while (c != '\n' && c != '\r' && c != EOF)
			c = getc(stream);
		if (c == '\r') {
			c = getc(stream);
			if (c != '\n')
				ungetc(c, stream);
		}
	}
	if (lineNumber != NULL)
		(*lineNumber)++;
	return i;
}

int
ReadFormat(const char *str, const char *fmt,...)
{
	va_list ap;
	char buf[64];
	int c, n, count, len;
	Int *ip;
	Double *fp;
	char *sp;
	va_start(ap, fmt);
	count = 0;
	len = strlen(str);
	while (*fmt != 0 && *str != 0) {
		if (isspace(*fmt)) {
			fmt++;
			continue;
		}
		c = tolower(*fmt++);
		if (isdigit(*fmt)) {
			n = strtol(fmt, (char **)&fmt, 0);
			if (n > 63)
				n = 63;
			else if (n < 1)
				n = 1;
		} else n = 1;
		if (len < n)
			n = len;
		strncpy(buf, str, n);
		buf[n] = 0;
		str += n;
		len -= n;
		switch (c) {
			case 'i':
				ip = va_arg(ap, Int *);
				*ip = atoi(buf);
				count++;
				break;
			case 'f':
				fp = va_arg(ap, Double *);
				*fp = atof(buf);
				count++;
				break;
			case 's':
				sp = va_arg(ap, char *);
				sscanf(buf, " %s", sp);
				count++;
				break;
			default:
				break;
		}
	}
	va_end(ap);
	return count;
}
