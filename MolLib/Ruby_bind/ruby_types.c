/*
 *  ruby_types.c
 *  Molby
 *
 *  Created by Toshi Nagata on 09/01/24.
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

#include "Molby.h"

#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#pragma mark ====== Global Values ======

VALUE rb_cVector3D, rb_cTransform, rb_cIntGroup;

#pragma mark ====== Utility functions (Vector/Matrix) ======

/*
static VALUE
s_VectorClass(void)
{
	if (rb_cVector == 0)
		rb_cVector = rb_const_get(rb_cObject, rb_intern("Vector"));
	return rb_cVector;
}
*/

/*
VALUE
ValueFromVector(const Vector *vp)
{
	static ID mname = 0;
	if (mname == 0)
		mname = rb_intern("[]");
	return rb_funcall(rb_cVector3D, mname, 3, rb_float_new(vp->x), rb_float_new(vp->y), rb_float_new(vp->z));
}
*/

#pragma mark ====== Vector3D Class ======

void
VectorFromValue(VALUE val, Vector *vp)
{
	Vector *vp1;
	if (rb_obj_is_kind_of(val, rb_cVector3D)) {
		Data_Get_Struct(val, Vector, vp1);
		*vp = *vp1;
	} else {
		static ID mname = 0;
		if (mname == 0)
			mname = rb_intern("[]");
		vp->x = NUM2DBL(rb_Float(rb_funcall(val, mname, 1, INT2FIX(0))));
		vp->y = NUM2DBL(rb_Float(rb_funcall(val, mname, 1, INT2FIX(1))));
		vp->z = NUM2DBL(rb_Float(rb_funcall(val, mname, 1, INT2FIX(2))));
	}
}

static VALUE
s_Vector3D_Alloc(VALUE klass)
{
	Vector *vp = ALLOC(Vector);
	vp->x = vp->y = vp->z = 0.0;
	return Data_Wrap_Struct(klass, 0, -1, vp);
}

VALUE
ValueFromVector(const Vector *vp)
{
	Vector *vp1;
	VALUE retval = s_Vector3D_Alloc(rb_cVector3D);
	Data_Get_Struct(retval, Vector, vp1);
	*vp1 = *vp;
	return retval;
}

/*
 *  call-seq:
 *     Vector3D.new
 *     Vector3D.new(vector3d)
 *     Vector3D.new(ary)
 *
 *  Returns a new Vector3D object. In the first form, a zero vector
 *  is returned. In the second form, the given vector3d is duplicated.
 *  In the third form, a vector [ary[0], ary[1], ary[2]] is returned.
 *  The argument ary can be anything that responds to '[]' method.
 */
static VALUE
s_Vector3D_Initialize(int argc, VALUE *argv, VALUE self)
{
	VALUE val;
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	rb_scan_args(argc, argv, "01", &val);
	if (!NIL_P(val))
		VectorFromValue(val, vp);
	return Qnil;
}

/*
 *  call-seq:
 *     vector3d.size -> int
 *  
 *  Returns 3. This method is present only to be consistent with classes like
 *  Array or Vector.
 */
static VALUE
s_Vector3D_Size(VALUE self)
{
	return INT2FIX(3);
}

/* 
 *  call-seq:
 *     vector3d[index]             -> float
 *
 *  Element Reference---Returns the element at _index_. If _index_ is
 *  less than 0 or more than 2, an exception is thrown.
 */
static VALUE
s_Vector3D_ElementAtIndex(VALUE self, VALUE val)
{
	Vector *vp;
	Double w;
	int n = NUM2INT(val);
	Data_Get_Struct(self, Vector, vp);
	if (n < 0 || n >= 3)
		rb_raise(rb_eMolbyError, "index to Vector3D out of range");
	w = (n == 0 ? vp->x : (n == 1 ? vp->y : vp->z));
	return rb_float_new(w);
}

/* 
 *  call-seq:
 *     vector3d[index] = float
 *
 *  Element Assignment---Set the element at _index_. If _index_ is
 *  less than 0 or more than 2, an exception is thrown.
 */
static VALUE
s_Vector3D_SetElementAtIndex(VALUE self, VALUE idx, VALUE val)
{
	Vector *vp;
	Double w = NUM2DBL(rb_Float(val));
	int n = NUM2INT(idx);
	Data_Get_Struct(self, Vector, vp);
	if (n < 0 || n >= 3)
		rb_raise(rb_eMolbyError, "index to Vector3D out of range");
	if (n == 0)
		vp->x = w;
	else if (n == 1)
		vp->y = w;
	else
		vp->z = w;
	return rb_float_new(w);
}

/* 
 *  call-seq:
 *     vector3d == other_vector3d   ->   bool
 *
 *  Equality---Two vector3ds are equal if their elements are all equal.
 *  Usual caution about comparison between floating point numbers should be
 *  paid. Also consider using something like <code>vector3d.length < 1e-10</code>.
 */
static VALUE
s_Vector3D_IsEqual(VALUE self, VALUE val)
{
	Vector *vp1, v2;
	Data_Get_Struct(self, Vector, vp1);
	VectorFromValue(val, &v2);
	if (vp1->x == v2.x && vp1->y == v2.y && vp1->z == v2.z)
		return Qtrue;
	else return Qfalse;
}

/* 
 *  call-seq:
 *     vector3d + other_vector3d   ->   (new) vector3d
 *
 *  Add two vectors element by element.
 */
static VALUE
s_Vector3D_Add(VALUE self, VALUE val)
{
	Vector *vp1, v2;
	VALUE retval;
	Data_Get_Struct(self, Vector, vp1);
	VectorFromValue(val, &v2);
	retval = s_Vector3D_Alloc(rb_cVector3D);
	v2.x += vp1->x;
	v2.y += vp1->y;
	v2.z += vp1->z;
	return ValueFromVector(&v2);
}

/* 
 *  call-seq:
 *     vector3d - other_vector3d   ->   (new) vector3d
 *
 *  Subtract two vectors element by element.
 */
static VALUE
s_Vector3D_Subtract(VALUE self, VALUE val)
{
	Vector *vp1, v2;
	VALUE retval;
	Data_Get_Struct(self, Vector, vp1);
	VectorFromValue(val, &v2);
	retval = s_Vector3D_Alloc(rb_cVector3D);
	v2.x = vp1->x - v2.x;
	v2.y = vp1->y - v2.y;
	v2.z = vp1->z - v2.z;
	return ValueFromVector(&v2);
}

/* 
 *  call-seq:
 *     vector3d.dot(other_vector3d) ->  float
 *
 *  Calculate the dot (inner) product of the two vectors. See also <code>vector3d.*</code>.
 */
static VALUE
s_Vector3D_Dot(VALUE self, VALUE val)
{
	Vector *vp1, v2;
	Data_Get_Struct(self, Vector, vp1);
	VectorFromValue(val, &v2);
	return rb_float_new(vp1->x * v2.x + vp1->y * v2.y + vp1->z * v2.z);
}

/* 
 *  call-seq:
 *     vector3d * numeric           ->  (new) vector3d
 *     vector3d * other_vector3d    ->  float (the dot product)
 *
 *  In the first form, the vector is scaled by the numeric. In the second
 *  form, the dot (inner) product of the two vectors are returned (equivalent to
 *  <code>vector3d.dot(other_vector3d)</code>).
 */
static VALUE
s_Vector3D_Multiply(VALUE self, VALUE val)
{
	Vector *vp1, v2;
	Data_Get_Struct(self, Vector, vp1);
	if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		double w = NUM2DBL(rb_Float(val));
		v2.x = vp1->x * w;
		v2.y = vp1->y * w;
		v2.z = vp1->z * w;
		return ValueFromVector(&v2);
	} else return s_Vector3D_Dot(self, val);
}

/* 
 *  call-seq:
 *     vector3d / numeric           ->  (new) vector3d
 *
 *  The vector is scaled by the inverse of the given numeric.
 */
static VALUE
s_Vector3D_Divide(VALUE self, VALUE val)
{
	Vector *vp1, v2;
	double w = NUM2DBL(rb_Float(val));
	Data_Get_Struct(self, Vector, vp1);
	v2.x = vp1->x / w;
	v2.y = vp1->y / w;
	v2.z = vp1->z / w;
	return ValueFromVector(&v2);
}

/* 
 *  call-seq:
 *     vector3d.cross(other_vector3d) ->  (new) vector3d
 *
 *  Calculate the cross (outer) product of the two vectors.
 */
static VALUE
s_Vector3D_Cross(VALUE self, VALUE val)
{
	Vector *vp1, v2, v3;
	Data_Get_Struct(self, Vector, vp1);
	VectorFromValue(val, &v2);
	v3.x = vp1->y * v2.z - vp1->z * v2.y;
	v3.y = vp1->z * v2.x - vp1->x * v2.z;
	v3.z = vp1->x * v2.y - vp1->y * v2.x;
	return ValueFromVector(&v3);
}

/* 
 *  call-seq:
 *     vector3d.-@                ->  vector3d
 *
 *  Calculate the opposite vector. 
 */
static VALUE
s_Vector3D_UnaryMinus(VALUE self)
{
	Vector *vp, v1;
	Data_Get_Struct(self, Vector, vp);
	v1.x = -vp->x;
	v1.y = -vp->y;
	v1.z = -vp->z;
	return ValueFromVector(&v1);
}

/* 
 *  call-seq:
 *     vector3d.length                ->  float
 *
 *  Calculate the Pythagorean length of the vector. 
 *  Note that this method is <em>not</em> an alias of <code>vector3d.size</code>.
 */
static VALUE
s_Vector3D_Length(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	return rb_float_new(sqrt(vp->x * vp->x + vp->y * vp->y + vp->z * vp->z));
}

/* 
 *  call-seq:
 *     vector3d.length2                ->  float
 *
 *  Calculate the square of the Pythagorean length of the vector.
 */
static VALUE
s_Vector3D_Length2(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	return rb_float_new(vp->x * vp->x + vp->y * vp->y + vp->z * vp->z);
}

/* 
 *  call-seq:
 *     vector3d.normalize              ->  (new) vector3d
 *
 *  Returns a unit vector with the same direction. Raises an exception when the
 *  vector is a zero vector.
 */
static VALUE
s_Vector3D_Normalize(VALUE self)
{
	Vector *vp, v;
	double w;
	Data_Get_Struct(self, Vector, vp);
	w = 1.0 / sqrt(vp->x * vp->x + vp->y * vp->y + vp->z * vp->z);
	if (!isfinite(w))
		rb_raise(rb_eMolbyError, "trying to normalize a (nearly) zero vector");
	v.x = vp->x * w;
	v.y = vp->y * w;
	v.z = vp->z * w;
	return ValueFromVector(&v);
}

/* 
 *  call-seq:
 *     vector3d.to_a                    ->  Array
 *
 *  Returns <code>[self.x, self.y, self.z]</code>.
 */
static VALUE
s_Vector3D_ToArray(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);	
	return rb_ary_new3(3, rb_float_new(vp->x), rb_float_new(vp->y), rb_float_new(vp->z));
}

/* 
 *  call-seq:
 *     vector3d.each {|item| block }  -> vector3d (self)
 *
 *  Calls <i>block</i> once for x, y, z elements, passing that element as a parameter.
 */
static VALUE
s_Vector3D_Each(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	rb_yield(rb_float_new(vp->x));
	rb_yield(rb_float_new(vp->y));
	rb_yield(rb_float_new(vp->z));
	return self;
}

/* 
 *  call-seq:
 *     vector3d.x       -> float
 *
 *  Get the x element of the vector.
 */
static VALUE
s_Vector3D_GetX(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	return rb_float_new(vp->x);
}

/* 
 *  call-seq:
 *     vector3d.y       -> float
 *
 *  Get the y element of the vector.
 */
static VALUE
s_Vector3D_GetY(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	return rb_float_new(vp->y);
}

/* 
 *  call-seq:
 *     vector3d.z       -> float
 *
 *  Get the z element of the vector.
 */
static VALUE
s_Vector3D_GetZ(VALUE self)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	return rb_float_new(vp->z);
}

/* 
 *  call-seq:
 *     vector3d.x = float       -> float
 *
 *  Set the x element of the vector.
 */
static VALUE
s_Vector3D_SetX(VALUE self, VALUE val)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	vp->x = NUM2DBL(rb_Float(val));
	return rb_float_new(vp->x);
}

/* 
 *  call-seq:
 *     vector3d.y = float       -> float
 *
 *  Set the y element of the vector.
 */
static VALUE
s_Vector3D_SetY(VALUE self, VALUE val)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	vp->y = NUM2DBL(rb_Float(val));
	return rb_float_new(vp->y);
}

/* 
 *  call-seq:
 *     vector3d.z = float       -> float
 *
 *  Set the z element of the vector.
 */
static VALUE
s_Vector3D_SetZ(VALUE self, VALUE val)
{
	Vector *vp;
	Data_Get_Struct(self, Vector, vp);
	vp->z = NUM2DBL(rb_Float(val));
	return rb_float_new(vp->z);
}

/* 
 *  call-seq:
 *     Vector3d[fx, fy, fz]   -> (new) vector3d
 *
 *  Create a new vector3d object. Equivalent to <code>Vector3d.new([fx, fy, fz])</code>.
 */
static VALUE
s_Vector3D_Create(VALUE klass, VALUE args)
{
	VALUE val = s_Vector3D_Alloc(klass);
	s_Vector3D_Initialize(1, &args, val);
	return val;
}

/* 
 *  call-seq:
 *     vector3d.inspect        -> string
 *
 *  Create a readable string like "Vector3d[fx, fy, fz]".
 */
static VALUE
s_Vector3D_Inspect(VALUE self)
{
	/*  self.class.name << self.to_a.inspect  */
	VALUE klass = CLASS_OF(self);
	VALUE val = rb_funcall(klass, rb_intern("name"), 0);
	self = s_Vector3D_ToArray(self);
	return rb_funcall(val, rb_intern("<<"), 1, rb_funcall(self, rb_intern("inspect"), 0));
}

#pragma mark ====== Transform Class ======

static int s_index_order[16] = {0, 3, 6, 1, 4, 7, 2, 5, 8, 9, 10, 11};

void
TransformFromValue(VALUE val, Transform *tp)
{
	int i, j;
	Transform *tp1;
	if (rb_obj_is_kind_of(val, rb_cTransform)) {
		Data_Get_Struct(val, Transform, tp1);
		memmove(tp, tp1, sizeof(Transform));
	} else if (TYPE(val) == T_ARRAY) {
		/*  Array must be 
		 (1) [[a11 a21 a31] [a12 a22 a32] [a13 a23 a33] [a14 a24 a34]] or
		 (2) [a11 a21 a31 a12 a22 a32 a13 a23 a33 a14 a24 a34] */
		int len = RARRAY_LEN(val);
		VALUE *valp = RARRAY_PTR(val);
		if (len == 12) {
			for (i = 0; i < 12; i++)
				(*tp)[s_index_order[i]] = NUM2DBL(rb_Float(valp[i]));
			return;
		} else if (len == 4) {
			VALUE val2, *valp2;
			for (i = 0; i < 4; i++) {
				val2 = valp[i];
				if (TYPE(val2) != T_ARRAY)
					val2 = rb_funcall(val2, rb_intern("to_a"), 0);
				if (TYPE(val2) != T_ARRAY || RARRAY_LEN(val2) != 3)
					goto array_format_error;
				valp2 = RARRAY_PTR(val2);
				for (j = 0; j < 3; j++)
					(*tp)[s_index_order[3 * i + j]] = NUM2DBL(rb_Float(valp2[j]));
			}
			return;
		}
	array_format_error:
		rb_raise(rb_eMolbyError, "wrong array format; must be an array of either (1) four 3D column vectors or (2) 12 numerics");
	} else {
		static ID index_mid = 0, row_size_mid, column_size_mid;
		if (index_mid == 0) {
			index_mid = rb_intern("[]");
			row_size_mid = rb_intern("row_size");
			column_size_mid = rb_intern("column_size");
		}
		if (rb_respond_to(val, row_size_mid) && rb_respond_to(val, column_size_mid)) {
			/*  Matrix-type object  */
			for (i = 0; i < 4; i++) {
				for (j = 0; j < 3; j++)
					(*tp)[s_index_order[i * 3 + j]] = NUM2DBL(rb_Float(rb_funcall(val, index_mid, 2, INT2FIX(i), INT2FIX(j))));
			}
		} else {
			/*  Other "array-like" object  */
			for (i = 0; i < 12; i++)
				(*tp)[s_index_order[i]] = NUM2DBL(rb_Float(rb_funcall(val, index_mid, 1, INT2FIX(i))));
		}
	}
}

static VALUE
s_Transform_Alloc(VALUE klass)
{
	Transform *tp = ALLOC(Transform);
	memset(tp, 0, sizeof(Transform));
	(*tp)[0] = (*tp)[4] = (*tp)[8] = 1.0;
	return Data_Wrap_Struct(klass, 0, -1, tp);
}

VALUE
ValueFromTransform(Transform *tp)
{
	Transform *tp1;
	VALUE val = s_Transform_Alloc(rb_cTransform);
	Data_Get_Struct(val, Transform, tp1);
	memmove(tp1, tp, sizeof(Transform));
	return val;
}

/*
 *  call-seq:
 *     Transform.new
 *     Transform.new(array)
 *     Transform.new(matrix)
 *
 *  Returns a new Transform object.
 *
 *  In the first form, an identity transform is returned.
 *
 *  In the second form, the array must be either of the following
 *  two forms:
 *  (1) [[a11 a21 a31] [a12 a22 a32] [a13 a23 a33] [a14 a24 a34]]
 *  (2) [a11 a21 a31 a12 a22 a32 a13 a23 a33 a14 a24 a34]
 *  where [a11..a33] denotes the rotation part and [a14 a24 a34] denotes
 *  the translation part. All vectors in (1) are column vectors.
 *  
 *  In the third form, a new transform is built from a 3x4 matrix. The argument
 *  <code>matrix</code> must respond to a method call <code>matrix[col, row]</code>
 *  where <code>row</code> is in <code>0..2</code> and <code>col</code> in <code>0..3</code>. 
 */
static VALUE
s_Transform_Initialize(int argc, VALUE *argv, VALUE self)
{
	VALUE val;
	Transform *tp;
	Data_Get_Struct(self, Transform, tp);
	rb_scan_args(argc, argv, "01", &val);
	if (!NIL_P(val))
		TransformFromValue(val, tp);
	return Qnil;
}

static VALUE
s_Transform_NewFromTransform(Transform *tp)
{
	Transform *tp1;
	VALUE retval = s_Transform_Alloc(rb_cTransform);
	Data_Get_Struct(retval, Transform, tp1);
	memmove(tp1, tp, sizeof(Transform));
	return retval;
}

/*
 *  call-seq:
 *     Transform.from_columns(c1, c2, c3, c4)
 *
 *  Returns a new Transform object built from four column vectors. The arguments
 *  <code>c1..c4</code> are vectors of (at least) three-dimension. This is equivalent
 *  to <code>Transform.new([c1, c2, c3, c4])</code>.
 */
static VALUE
s_Transform_NewFromColumns(VALUE klass, VALUE val)
{
	Transform tr;
	int i, j, n;
	VALUE *valp;
	static ID to_a_mid = 0;
	if (to_a_mid == 0)
		to_a_mid = rb_intern("to_a");
	if (TYPE(val) != T_ARRAY)
		val = rb_funcall(val, to_a_mid, 0);
	memset(tr, 0, sizeof(tr));
	n = RARRAY_LEN(val);
	valp = RARRAY_PTR(val);
	for (i = 0; i < 4; i++) {
		int nn;
		VALUE *valpp;
		double w[3];
		w[0] = w[1] = w[2] = 0.0;
		if (i < n) {
			val = valp[i];
			if (TYPE(val) != T_ARRAY)
				val = rb_funcall(val, to_a_mid, 0);
			nn = RARRAY_LEN(val);
			valpp = RARRAY_PTR(val);
			for (j = 0; j < 3 && j < nn; j++)
				w[j] = NUM2DBL(rb_Float(valpp[j]));
		}
		for (j = 0; j < 3; j++)
			tr[s_index_order[i * 3 + j]] = w[j];
	}
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     Transform.from_rows(r1, r2, r3)
 *
 *  Returns a new Transform object built from three row vectors. The arguments
 *  <code>r1, r2, r3</code> are vectors of (at least) four-dimension.
 */
static VALUE
s_Transform_NewFromRows(VALUE klass, VALUE val)
{
	Transform tr;
	int i, j, n;
	VALUE *valp;
	static ID to_a_mid = 0;
	if (to_a_mid == 0)
		to_a_mid = rb_intern("to_a");
	if (TYPE(val) != T_ARRAY)
		val = rb_funcall(val, to_a_mid, 0);
	memset(tr, 0, sizeof(tr));
	n = RARRAY_LEN(val);
	valp = RARRAY_PTR(val);
	for (i = 0; i < 3; i++) {
		int nn;
		VALUE *valpp;
		double w[4];
		w[0] = w[1] = w[2] = w[3] = 0.0;
		if (i < n) {
			val = valp[i];
			if (TYPE(val) != T_ARRAY)
				val = rb_funcall(val, to_a_mid, 0);
			nn = RARRAY_LEN(val);
			valpp = RARRAY_PTR(val);
			for (j = 0; j < 4 && j < nn; j++)
				w[j] = NUM2DBL(rb_Float(valpp[j]));
		}
		for (j = 0; j < 4; j++)
			tr[s_index_order[j * 3 + i]] = w[j];
	}
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform[i, j]  -> float
 *
 *  Get the element (+i+,+j+) of the transform matrix, i.e. column +i+, row +j+.
 *  Be careful about the order of the arguments. It follows convention of multi-dimensional arrays
 *  rather than mathematical notation.
 */
static VALUE
s_Transform_ElementAtIndex(VALUE self, VALUE val1, VALUE val2)
{
	Transform *tp;
	double w;
	int n1 = NUM2INT(val1);
	int n2 = NUM2INT(val2);
	Data_Get_Struct(self, Transform, tp);
	if (n1 < 0 || n1 >= 4 || n2 < 0 || n2 >= 3)
		rb_raise(rb_eMolbyError, "index to Transform out of range");
	w = (*tp)[s_index_order[n1 * 3 + n2]];
	return rb_float_new(w);
}

/*
 *  call-seq:
 *     transform[i, j] = float  -> float
 *
 *  Set the element (+i+,+j+) of the transform matrix, i.e. column +i+, row +j+.
 *  Be careful about the order of the arguments. It follows convention of multi-dimensional arrays
 *  rather than mathematical notation.
 */
static VALUE
s_Transform_SetElementAtIndex(VALUE self, VALUE idx1, VALUE idx2, VALUE val)
{
	Transform *tp;
	double w;
	int n1 = NUM2INT(idx1);
	int n2 = NUM2INT(idx2);
	Data_Get_Struct(self, Transform, tp);
	if (n1 < 0 || n1 >= 4 || n2 < 0 || n2 >= 3)
		rb_raise(rb_eMolbyError, "index to Transform out of range");
	w = NUM2DBL(rb_Float(val));
	(*tp)[s_index_order[n1 * 3 + n2]] = w;
	return rb_float_new(w);
}

/*
 *  call-seq:
 *     transform == other_transform  -> bool
 *
 *  Returns +true+ if and only if all the corresponding elements are equal.
 *  Usual caution about the comparison of floating-point numbers should be paid.
 */
static VALUE
s_Transform_IsEqual(VALUE self, VALUE val)
{
	Transform *tp1, tr;
	int i;
	Data_Get_Struct(self, Transform, tp1);
	TransformFromValue(val, &tr);
	for (i = 0; i < 12; i++) {
		if ((*tp1)[i] != tr[i])
			return Qfalse;
	}
	return Qtrue;
}

/*
 *  call-seq:
 *     transform + other_transform  -> (new) transform
 *
 *  Returns a new transform corresponding to the sum of the two transform matrix.
 */
static VALUE
s_Transform_Add(VALUE self, VALUE val)
{
	Transform *tp1, tr;
	int i;
	Data_Get_Struct(self, Transform, tp1);
	TransformFromValue(val, &tr);
	for (i = 0; i < 12; i++)
		tr[i] += (*tp1)[i];
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform - other_transform  -> (new) transform
 *
 *  Returns a new transform corresponding to the difference of the two transform matrix.
 */
static VALUE
s_Transform_Subtract(VALUE self, VALUE val)
{
	Transform *tp1, tr;
	int i;
	Data_Get_Struct(self, Transform, tp1);
	TransformFromValue(val, &tr);
	for (i = 0; i < 12; i++)
		tr[i] = (*tp1)[i] - tr[i];
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform * numeric          -> (new) transform
 *     transform * vector3d         -> (new) vector3d
 *     transform * other_transform  -> (new) transform
 *
 *  Perform the matrix multiplication. In the first form, a new matrix with scaled elements
 *  is returned. In the second, the transformed vector is returned. In the third form,
 *  the multiple of the two matrices is returned.
 */
static VALUE
s_Transform_Multiply(VALUE self, VALUE val)
{
	Transform *tp1, tr;
	int i;
	Data_Get_Struct(self, Transform, tp1);
	if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		double w = NUM2DBL(rb_Float(val));
		for (i = 0; i < 12; i++)
			tr[i] = (*tp1)[i] * w;
		return s_Transform_NewFromTransform(&tr);
	} else {
		static ID size_mid = 0;
		if (size_mid == 0)
			size_mid = rb_intern("size");
		if (rb_respond_to(val, size_mid) && NUM2INT(rb_funcall(val, size_mid, 0)) == 3) {
			/*  A 3D vector  */
			Vector v;
			VectorFromValue(val, &v);
			TransformVec(&v, *tp1, &v);
			return ValueFromVector(&v);
		} else {
			/*  Transform  */
			TransformFromValue(val, &tr);
			TransformMul(tr, *tp1, tr);
			return s_Transform_NewFromTransform(&tr);
		}
	}
}

/*
 *  call-seq:
 *     Transform.identity  -> transform
 *
 *  Returns an identity transform, <code>[[1,0,0], [0,1,0], [0,0,1], [0,0,0]]</code>.
 */
static VALUE
s_Transform_Identity(VALUE klass)
{
	Transform tr;
	memset(tr, 0, sizeof(tr));
	tr[0] = tr[4] = tr[8] = 1.0;
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     Transform.zero  -> transform
 *
 *  Returns a zero transform, <code>[[0,0,0], [0,0,0], [0,0,0], [0,0,0]]</code>.
 */
static VALUE
s_Transform_Zero(VALUE klass)
{
	Transform tr;
	memset(tr, 0, sizeof(tr));
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     Transform.diagonal(array)
 *     Transform.diagonal(f1, f2 = nil, f3 = nil)
 *
 *  Returns a diagonal transform (the translational componets are all zero).
 *  In the first form, <code>array[0], array[1], array[2]</code> are for the
 *  x, y, z components, respectively. In the second form, <code>f1, f2, f3</code>
 *  are the x, y, z components. If <code>f3</code> is not given, the <code>f2</code>
 *  is used for the z components. If <code>f2</code> is not given, the <code>f1</code>
 *  is used for the y and z components.
 */
static VALUE
s_Transform_Diagonal(int argc, VALUE *argv, VALUE klass)
{
	Transform tr;
	VALUE arg1, arg2, arg3;
	memset(tr, 0, sizeof(tr));
	rb_scan_args(argc, argv, "12", &arg1, &arg2, &arg3);
	if (TYPE(arg1) == T_ARRAY) {
		/*  [d1, d2, d3]  */
		VALUE *valp;
		int len;
		len = RARRAY_LEN(arg1);
		valp = RARRAY_PTR(arg1);
		if (len >= 1)
			tr[0] = tr[4] = tr[8] = NUM2DBL(rb_Float(valp[0]));
		if (len >= 2)
			tr[4] = tr[8] = NUM2DBL(rb_Float(valp[1]));
		if (len >= 3)
			tr[8] = NUM2DBL(rb_Float(valp[2]));
	} else {
		tr[0] = tr[4] = tr[8] = NUM2DBL(rb_Float(arg1));
		if (!NIL_P(arg2)) {
			tr[4] = tr[8] = NUM2DBL(rb_Float(arg2));
			if (!NIL_P(arg3))
				tr[8] = NUM2DBL(rb_Float(arg3));
		}
	}
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform.inverse  -> (new) transform
 *
 *  Returns the inverse transform. If the matrix is not regular, an exception is raised.
 */
static VALUE
s_Transform_Inverse(VALUE self)
{
	Transform *tp1, tr;
	Data_Get_Struct(self, Transform, tp1);
	if (TransformInvert(tr, *tp1))
		rb_raise(rb_eMolbyError, "the transform matrix is not regular");
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform / other_transform -> (new) transform
 *
 *  Returns transform * other_transform.invert. If other_transform is not regular, 
 *  an exception is raised.
 */
static VALUE
s_Transform_Divide(VALUE self, VALUE val)
{
	Transform *tp1, tr;
	Data_Get_Struct(self, Transform, tp1);
	TransformFromValue(val, &tr);
	if (TransformInvert(tr, tr))
		rb_raise(rb_eMolbyError, "the transform matrix is not regular");
	TransformMul(tr, *tp1, tr);
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform.transpose -> (new) transform
 *
 *  Returns a new transform in which the rotation component is transposed from the original.
 */
static VALUE
s_Transform_Transpose(VALUE self)
{
	Transform *tp1, tr;
	Data_Get_Struct(self, Transform, tp1);
	TransformTranspose(tr, *tp1);
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     transform.determinant -> float
 *
 *  Returns the determinant of the transform.
 */
static VALUE
s_Transform_Determinant(VALUE self)
{
	Transform *tp1;
	Data_Get_Struct(self, Transform, tp1);
	return rb_float_new(TransformDeterminant(*tp1));
}

/*
 *  call-seq:
 *     transform.trace -> float
 *
 *  Returns the trace (sum of the diagonal elements) of the transform.
 */
static VALUE
s_Transform_Trace(VALUE self)
{
	Transform *tp1;
	Data_Get_Struct(self, Transform, tp1);
	return rb_float_new((*tp1)[0] + (*tp1)[4] + (*tp1)[8]);
}

/*
 *  call-seq:
 *     transform.column(index) -> vector3d
 *
 *  Returns the index-th (0..3) column vector.
 */
static VALUE
s_Transform_Column(VALUE self, VALUE val)
{
	Transform *tp1;
	Vector v;
	int n = NUM2INT(val);
	Data_Get_Struct(self, Transform, tp1);
	if (n < 0 || n >= 4)
		rb_raise(rb_eMolbyError, "row index out of range");
	v.x = (*tp1)[n];
	v.y = (*tp1)[n + 3];
	v.z = (*tp1)[n + 6];
	return ValueFromVector(&v);
}

/*
 *  call-seq:
 *     transform.eigenvalues -> [[k1, k2, k3], v1, v2, v3]
 *
 *  Calculate the eigenvalues and eigenvectors. The matrix must be symmetric.
 */
static VALUE
s_Transform_Eigenvalues(VALUE self)
{
	Transform *tp1;
	Vector v[3];
	Double d[3];
	int info;
	Data_Get_Struct(self, Transform, tp1);
	if ((info = MatrixSymDiagonalize(*((Mat33 *)tp1), d, v)) != 0)
		rb_raise(rb_eMolbyError, "cannot diagonalize the given matrix: info = %d", info);
	return rb_ary_new3(4, rb_ary_new3(3, rb_float_new(d[0]), rb_float_new(d[1]), rb_float_new(d[2])), ValueFromVector(v), ValueFromVector(v + 1), ValueFromVector(v + 2));
}

/*
 *  call-seq:
 *     Transform[*args]
 *
 *  Create a new transform. Equivalent to Transform.new(args).
 */
static VALUE
s_Transform_Create(VALUE klass, VALUE args)
{
	VALUE val = s_Transform_Alloc(klass);
	s_Transform_Initialize(1, &args, val);
	return val;
}

/*
 *  call-seq:
 *     transform.to_a  -> array
 *
 *  Convert a transform to an array of 12 float numbers.
 */
static VALUE
s_Transform_ToArray(VALUE self)
{
	Transform *tp1;
	VALUE val[12];
	int i;
	Data_Get_Struct(self, Transform, tp1);
	for (i = 0; i < 12; i++)
		val[i] = rb_float_new((*tp1)[s_index_order[i]]);
	return rb_ary_new4(12, val);
}

/*
 *  call-seq:
 *     transform.inspect  -> string
 *
 *  Convert a transform to a string like 
 *  "Transform[[a11,a21,a31],[a12,a22,a32],[a13,a23,a33],[a14,a24,a34]]".
 */
static VALUE
s_Transform_Inspect(VALUE self)
{
	Transform *tp;
	int i, j;
	VALUE klass = CLASS_OF(self);
	VALUE val = rb_funcall(klass, rb_intern("name"), 0);
	ID mid = rb_intern("<<");
	ID mid2 = rb_intern("inspect");
	rb_str_cat(val, "[", 1);
	Data_Get_Struct(self, Transform, tp);
	for (i = 0; i < 4; i++) {
		rb_str_cat(val, "[", 1);
		for (j = 0; j < 3; j++) {
			double f;
			f = (*tp)[s_index_order[i * 3 + j]];
			rb_funcall(val, mid, 1, rb_funcall(rb_float_new(f), mid2, 0));
			if (j < 2)
				rb_str_cat(val, ",", 1);
		}
		rb_str_cat(val, "]", 1);
		if (i < 3)
			rb_str_cat(val, ",", 1);
	}
	rb_str_cat(val, "]", 1);
	return val;
}

/*
 *  call-seq:
 *     Transform.translation(vec)
 *
 *  Returns a transform corresponding to translation along the given vector. Equivalent
 *  to <code>Transform[[1,0,0],[0,1,0],[0,0,1],vec]</code>.
 */
static VALUE
s_Transform_Translation(VALUE klass, VALUE vec)
{
	Transform *tp;
	Vector v;
	VALUE val = s_Transform_Alloc(klass);
	Data_Get_Struct(val, Transform, tp);
	VectorFromValue(vec, &v);
	(*tp)[9] = v.x;
	(*tp)[10] = v.y;
	(*tp)[11] = v.z;
	return val;
}

/*
 *  call-seq:
 *     Transform.rotation(axis, angle, center = [0,0,0])
 *
 *  Returns a transform corresponding to the rotation along the given axis and angle. If 
 *  center is also given, that point will be the center of rotation.
 */
static VALUE
s_Transform_Rotation(int argc, VALUE *argv, VALUE klass)
{
	Transform tr;
	VALUE axis, angle, center;
	Vector av, cv;
	double ang;
	rb_scan_args(argc, argv, "21", &axis, &angle, &center);
	VectorFromValue(axis, &av);
	if (NIL_P(center))
		cv.x = cv.y = cv.z = 0.0;
	else
		VectorFromValue(center, &cv);
	ang = NUM2DBL(rb_Float(angle));
	if (TransformForRotation(tr, &av, ang, &cv))
		rb_raise(rb_eMolbyError, "rotation axis cannot be a zero vector");
	return ValueFromTransform(&tr);
}

/*
 *  call-seq:
 *     Transform.reflection(axis, center = [0,0,0])
 *
 *  Returns a transform corresponding to the reflection along the given axis. If 
 *  center is also given, that point will be fixed.
 */
static VALUE
s_Transform_Reflection(int argc, VALUE *argv, VALUE klass)
{
	VALUE axis, center;
	Vector av, cv;
	Transform tr;
	rb_scan_args(argc, argv, "11", &axis, &center);
	VectorFromValue(axis, &av);
	if (NIL_P(center))
		cv.x = cv.y = cv.z = 0.0;
	else
		VectorFromValue(center, &cv);
	if (TransformForReflection(tr, &av, &cv))
		rb_raise(rb_eMolbyError, "reflection axis cannot be a zero vector");
	return ValueFromTransform(&tr);
}

/*
 *  call-seq:
 *     Transform.inversion(center = [0,0,0])
 *
 *  Returns a transform corresponding to the inversion along the given point.
 */
static VALUE
s_Transform_Inversion(int argc, VALUE *argv, VALUE klass)
{
	VALUE center;
	Vector cv;
	Transform tr;
	rb_scan_args(argc, argv, "01", &center);
	if (NIL_P(center))
		cv.x = cv.y = cv.z = 0.0;
	else
		VectorFromValue(center, &cv);
	TransformForInversion(tr, &cv);
	return ValueFromTransform(&tr);
}

#pragma mark ====== IntGroup Class ======

IntGroup *
IntGroupFromValue(VALUE val)
{
	IntGroup *ig;
	if (!rb_obj_is_kind_of(val, rb_cIntGroup))
		val = rb_funcall(rb_cIntGroup, rb_intern("new"), 1, val);
	Data_Get_Struct(val, IntGroup, ig);
	IntGroupRetain(ig);
	return ig;
}

VALUE
ValueFromIntGroup(IntGroup *ig)
{
	if (ig == NULL)
		return Qnil;
	IntGroupRetain(ig);
    return Data_Wrap_Struct(rb_cIntGroup, 0, (void (*)(void *))IntGroupRelease, ig);
}

void
IntGroup_RaiseIfError(int err)
{
	if (err != 0) {
		const char *s;
		switch (err) {
			case kIntGroupStatusOutOfMemory: s = "out of memory"; break;
			case kIntGroupStatusOutOfRange: s = "out of range"; break;
			default: s = ""; break;
		}
		rb_raise(rb_eMolbyError, "%s error occurred during IntGroup operation", s);
	}
}

/*  Allocator  */
VALUE
IntGroup_Alloc(VALUE klass)
{
	IntGroup *ig = IntGroupNew();
    return Data_Wrap_Struct(klass, 0, (void (*)(void *))IntGroupRelease, ig);
}

/*  Iterator block for initializer  */
static VALUE
s_IntGroup_Initialize_i(VALUE val, VALUE ig1)
{
	IntGroup_RaiseIfError(IntGroupAdd((IntGroup *)ig1, NUM2INT(val), 1));
	return Qnil;
}

static VALUE
s_IntGroup_Initialize(int argc, VALUE *argv, VALUE self)
{
	IntGroup *ig1;
	Data_Get_Struct(self, IntGroup, ig1);
	while (argc-- > 0) {
		VALUE arg = *argv++;
		int type = TYPE(arg);
		if (rb_obj_is_kind_of(arg, rb_cIntGroup))
			rb_funcall(rb_cIntGroup, rb_intern("merge"), 1, arg);
		else if (rb_obj_is_kind_of(arg, rb_cRange)) {
			int sp, ep;
			sp = NUM2INT(rb_funcall(arg, rb_intern("begin"), 0));
			ep = NUM2INT(rb_funcall(arg, rb_intern("end"), 0));
			if (RTEST(rb_funcall(arg, rb_intern("exclude_end?"), 0)))
				ep--;
			if (ep >= sp)
				IntGroup_RaiseIfError(IntGroupAdd(ig1, sp, ep - sp + 1));
		} else if (rb_respond_to(arg, rb_intern("each")) && type != T_STRING)
			rb_iterate(rb_each, arg, s_IntGroup_Initialize_i, (VALUE)ig1);
		else
			IntGroup_RaiseIfError(IntGroupAdd(ig1, NUM2INT(arg), 1));
	}
	if (rb_block_given_p()) {
		IntGroup *ig2 = IntGroupNew();
		int i, n;
		for (i = 0; (n = IntGroupGetNthPoint(ig1, i)) >= 0; i++) {
			n = NUM2INT(rb_yield(INT2NUM(n)));
			if (n >= 0)
				IntGroup_RaiseIfError(IntGroupAdd(ig2, n, 1));
		}
		IntGroup_RaiseIfError(IntGroupCopy(ig1, ig2));
	}
	return Qnil;
}

static VALUE
s_IntGroup_Clear(VALUE self)
{
	IntGroup *ig;
	Data_Get_Struct(self, IntGroup, ig);
	IntGroupClear(ig);
	return self;
}

static VALUE
s_IntGroup_InitializeCopy(VALUE self, VALUE val)
{
	IntGroup *ig1, *ig2;
	Data_Get_Struct(self, IntGroup, ig1);
	if (!rb_obj_is_kind_of(val, rb_cIntGroup))
		rb_raise(rb_eMolbyError, "IntGroup instance is expected");
    Data_Get_Struct(val, IntGroup, ig2);
	IntGroupCopy(ig1, ig2);
	return self;
}

static VALUE
s_IntGroup_Length(VALUE self)
{
	IntGroup *ig;
	Data_Get_Struct(self, IntGroup, ig);
	return INT2NUM(IntGroupGetCount(ig));
}

static VALUE
s_IntGroup_MemberP(VALUE self, VALUE val)
{
	IntGroup *ig;
	int n = NUM2INT(val);
	Data_Get_Struct(self, IntGroup, ig);
	return (IntGroupLookup(ig, n, NULL) ? Qtrue : Qfalse);
}

static VALUE
s_IntGroup_ElementAtIndex(VALUE self, VALUE val)
{
	IntGroup *ig;
	int n;
	int index = NUM2INT(rb_Integer(val));
	Data_Get_Struct(self, IntGroup, ig);
	n = IntGroupGetNthPoint(ig, index);
	return (n >= 0 ? INT2NUM(n) : Qnil);
}

static VALUE
s_IntGroup_Each(VALUE self)
{
	IntGroup *ig;
	int i, j, sp, ep;
	Data_Get_Struct(self, IntGroup, ig);
	for (i = 0; (sp = IntGroupGetStartPoint(ig, i)) >= 0; i++) {
		ep = IntGroupGetEndPoint(ig, i);
		for (j = sp; j < ep; j++) {
			rb_yield(INT2NUM(j));
		}
	}
	return self;
}

static VALUE
s_IntGroup_Add(VALUE self, VALUE val)
{
	IntGroup *ig, *ig2;
    if (OBJ_FROZEN(self))
		rb_error_frozen("IntGroup");
	Data_Get_Struct(self, IntGroup, ig);
	if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		int n = NUM2INT(rb_Integer(val));
		if (n < 0)
			rb_raise(rb_eMolbyError, "the integer group can contain only non-negative values");
		IntGroupAdd(ig, n, 1);
	} else {
		ig2 = IntGroupFromValue(val);
		IntGroupAddIntGroup(ig, ig2);
		IntGroupRelease(ig2);
	}
	return self;
}

static VALUE
s_IntGroup_Delete(VALUE self, VALUE val)
{
	IntGroup *ig, *ig2;
    if (OBJ_FROZEN(self))
		rb_error_frozen("IntGroup");
	Data_Get_Struct(self, IntGroup, ig);
	if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		int n = NUM2INT(rb_Integer(val));
		if (n >= 0 && IntGroupLookup(ig, n, NULL))
			IntGroupRemove(ig, n, 1);
	} else {
		ig2 = IntGroupFromValue(val);
		IntGroupRemoveIntGroup(ig, ig2);
		IntGroupRelease(ig2);
	}
	return self;
}

static VALUE
s_IntGroup_Binary(VALUE self, VALUE val, int (*func)(const IntGroup *, const IntGroup *, IntGroup *))
{
	IntGroup *ig1, *ig2, *ig3;
	VALUE retval;
	Data_Get_Struct(self, IntGroup, ig1);
	ig2 = IntGroupFromValue(val);
	retval = IntGroup_Alloc(rb_cIntGroup);
	Data_Get_Struct(retval, IntGroup, ig3);
	IntGroup_RaiseIfError(func(ig1, ig2, ig3));
	IntGroupRelease(ig2);
	return retval;
}

static VALUE
s_IntGroup_Union(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupUnion);
}

static VALUE
s_IntGroup_Intersection(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupIntersect);
}

static VALUE
s_IntGroup_Difference(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupDifference);
}

static VALUE
s_IntGroup_SymDifference(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupXor);
}

static VALUE
s_IntGroup_Convolute(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupConvolute);
}

static VALUE
s_IntGroup_Deconvolute(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupDeconvolute);
}

static VALUE
s_IntGroup_RangeAt(VALUE self, VALUE val)
{
	IntGroup *ig;
	int n = NUM2INT(val);
	int sp, ep;
	Data_Get_Struct(self, IntGroup, ig);
	sp = IntGroupGetStartPoint(ig, n);
	if (sp < 0)
		return Qnil;
	ep = IntGroupGetEndPoint(ig, n) - 1;
	return rb_funcall(rb_cRange, rb_intern("new"), 2, INT2NUM(sp), INT2NUM(ep));
}

static VALUE
s_IntGroup_Merge(VALUE self, VALUE val)
{
	IntGroup *ig1, *ig2;
	int i, sp, interval;
    if (OBJ_FROZEN(self))
		rb_error_frozen("IntGroup");
	Data_Get_Struct(self, IntGroup, ig1);
	ig2 = IntGroupFromValue(val);
	for (i = 0; (sp = IntGroupGetStartPoint(ig2, i)) >= 0; i++) {
		interval = IntGroupGetInterval(ig2, i);
		IntGroup_RaiseIfError(IntGroupAdd(ig1, sp, interval));
	}
	IntGroupRelease(ig2);
	return self;
}

static VALUE
s_IntGroup_Subtract(VALUE self, VALUE val)
{
	IntGroup *ig1, *ig2;
	int i, sp, interval;
    if (OBJ_FROZEN(self))
		rb_error_frozen("IntGroup");
	Data_Get_Struct(self, IntGroup, ig1);
	ig2 = IntGroupFromValue(val);
	for (i = 0; (sp = IntGroupGetStartPoint(ig2, i)) >= 0; i++) {
		interval = IntGroupGetInterval(ig2, i);
		IntGroup_RaiseIfError(IntGroupRemove(ig1, sp, interval));
	}
	IntGroupRelease(ig2);
	return self;
}

static VALUE
s_IntGroup_Offset(VALUE self, VALUE ofs)
{
	IntGroup *ig1, *ig2;
	int iofs;
	VALUE val;
	Data_Get_Struct(self, IntGroup, ig1);
	ig2 = IntGroupNewFromIntGroup(ig1);
	if (ig2 == NULL)
		rb_raise(rb_eMolbyError, "Cannot duplicate IntGroup");
	iofs = NUM2INT(ofs);
	if (IntGroupOffset(ig2, iofs) != 0)
		rb_raise(rb_eMolbyError, "Bad offset %d", iofs);
	val = ValueFromIntGroup(ig2);
	IntGroupRelease(ig2);
	return val;
}

static VALUE
s_IntGroup_Create(int argc, VALUE *argv, VALUE klass)
{
	VALUE val = IntGroup_Alloc(klass);
	s_IntGroup_Initialize(argc, argv, val);
	return val;
}

static VALUE
s_IntGroup_Inspect(VALUE self)
{
	int i, sp, ep;
	IntGroup *ig;
	char buf[64];
	VALUE klass = CLASS_OF(self);
	VALUE val = rb_funcall(klass, rb_intern("name"), 0);
	Data_Get_Struct(self, IntGroup, ig);
	rb_str_cat(val, "[", 1);
	for (i = 0; (sp = IntGroupGetStartPoint(ig, i)) >= 0; i++) {
		if (i > 0)
			rb_str_cat(val, ", ", 2);
		ep = IntGroupGetEndPoint(ig, i);
		if (ep > sp + 1)
			snprintf(buf, sizeof buf, "%d..%d", sp, ep - 1);
		else
			snprintf(buf, sizeof buf, "%d", sp);
		rb_str_cat(val, buf, strlen(buf));
	}
	rb_str_cat(val, "]", 1);
	return val;
}

void
Init_MolbyTypes(void)
{
	/*  class Vector3D  */
	rb_cVector3D = rb_define_class("Vector3D", rb_cObject);
	rb_define_alloc_func(rb_cVector3D, s_Vector3D_Alloc);
	rb_define_method(rb_cVector3D, "initialize", s_Vector3D_Initialize, -1);
	rb_define_method(rb_cVector3D, "size", s_Vector3D_Size, 0);
	rb_define_method(rb_cVector3D, "[]", s_Vector3D_ElementAtIndex, 1);
	rb_define_method(rb_cVector3D, "[]=", s_Vector3D_SetElementAtIndex, 2);
	rb_define_method(rb_cVector3D, "==", s_Vector3D_IsEqual, 1);
	rb_define_method(rb_cVector3D, "+", s_Vector3D_Add, 1);
	rb_define_method(rb_cVector3D, "-", s_Vector3D_Subtract, 1);
	rb_define_method(rb_cVector3D, "*", s_Vector3D_Multiply, 1);
	rb_define_method(rb_cVector3D, "/", s_Vector3D_Divide, 1);
	rb_define_method(rb_cVector3D, "dot", s_Vector3D_Dot, 1);
	rb_define_method(rb_cVector3D, "cross", s_Vector3D_Cross, 1);
	rb_define_method(rb_cVector3D, "-@", s_Vector3D_UnaryMinus, 0);
	rb_define_method(rb_cVector3D, "length", s_Vector3D_Length, 0);
	rb_define_alias(rb_cVector3D, "r", "length"); /* size and length are not synonym!  */
	rb_define_method(rb_cVector3D, "length2", s_Vector3D_Length2, 0);
	rb_define_alias(rb_cVector3D, "r2", "length2");
	rb_define_method(rb_cVector3D, "normalize", s_Vector3D_Normalize, 0);
	rb_define_method(rb_cVector3D, "to_a", s_Vector3D_ToArray, 0);
	rb_define_method(rb_cVector3D, "each", s_Vector3D_Each, 0);
	rb_define_method(rb_cVector3D, "x", s_Vector3D_GetX, 0);
	rb_define_method(rb_cVector3D, "y", s_Vector3D_GetY, 0);
	rb_define_method(rb_cVector3D, "z", s_Vector3D_GetZ, 0);
	rb_define_method(rb_cVector3D, "x=", s_Vector3D_SetX, 1);
	rb_define_method(rb_cVector3D, "y=", s_Vector3D_SetY, 1);
	rb_define_method(rb_cVector3D, "z=", s_Vector3D_SetZ, 1);
	rb_define_method(rb_cVector3D, "inspect", s_Vector3D_Inspect, 0);
	rb_define_alias(rb_cVector3D, "to_s", "inspect");
	rb_define_singleton_method(rb_cVector3D, "[]", s_Vector3D_Create, -2);

	/*  class Transform  */
	rb_cTransform = rb_define_class("Transform", rb_cObject);
	rb_define_alloc_func(rb_cTransform, s_Transform_Alloc);
	rb_define_method(rb_cTransform, "initialize", s_Transform_Initialize, -1);
	rb_define_method(rb_cTransform, "[]", s_Transform_ElementAtIndex, 2);
	rb_define_method(rb_cTransform, "[]=", s_Transform_SetElementAtIndex, 3);
	rb_define_method(rb_cTransform, "==", s_Transform_IsEqual, 1);
	rb_define_method(rb_cTransform, "+", s_Transform_Add, 1);
	rb_define_method(rb_cTransform, "-", s_Transform_Subtract, 1);
	rb_define_method(rb_cTransform, "*", s_Transform_Multiply, 1);
	rb_define_method(rb_cTransform, "/", s_Transform_Divide, 1);
	rb_define_method(rb_cTransform, "inverse", s_Transform_Inverse, 0);
	rb_define_method(rb_cTransform, "transpose", s_Transform_Transpose, 0);
	rb_define_method(rb_cTransform, "determinant", s_Transform_Determinant, 0);
	rb_define_method(rb_cTransform, "trace", s_Transform_Trace, 0);
	rb_define_method(rb_cTransform, "column", s_Transform_Column, 1);
	rb_define_method(rb_cTransform, "eigenvalues", s_Transform_Eigenvalues, 0);
	rb_define_method(rb_cTransform, "to_a", s_Transform_ToArray, 0);
	rb_define_method(rb_cTransform, "inspect", s_Transform_Inspect, 0);
	rb_define_alias(rb_cTransform, "to_s", "inspect");
	rb_define_singleton_method(rb_cTransform, "diagonal", s_Transform_Diagonal, -1);
	rb_define_singleton_method(rb_cTransform, "[]", s_Transform_Create, -2);
	rb_define_singleton_method(rb_cTransform, "from_columns", s_Transform_NewFromColumns, -2);
	rb_define_singleton_method(rb_cTransform, "from_rows", s_Transform_NewFromRows, -2);
	rb_define_singleton_method(rb_cTransform, "identity", s_Transform_Identity, 0);
	rb_define_singleton_method(rb_cTransform, "zero", s_Transform_Zero, 0);
	rb_define_singleton_method(rb_cTransform, "translation", s_Transform_Translation, 1);
	rb_define_singleton_method(rb_cTransform, "rotation", s_Transform_Rotation, -1);
	rb_define_singleton_method(rb_cTransform, "reflection", s_Transform_Reflection, -1);
	rb_define_singleton_method(rb_cTransform, "inversion", s_Transform_Inversion, -1);

	/*  class IntGroup  */
	rb_cIntGroup = rb_define_class("IntGroup", rb_cObject);
	rb_include_module(rb_cIntGroup, rb_mEnumerable);
	rb_define_alloc_func(rb_cIntGroup, IntGroup_Alloc);
	rb_define_method(rb_cIntGroup, "clear", s_IntGroup_Clear, 0);
	rb_define_method(rb_cIntGroup, "initialize", s_IntGroup_Initialize, -1);
	rb_define_method(rb_cIntGroup, "initialize_copy", s_IntGroup_InitializeCopy, 1);
	rb_define_method(rb_cIntGroup, "length", s_IntGroup_Length, 0);
	rb_define_alias(rb_cIntGroup, "size", "length");
	rb_define_method(rb_cIntGroup, "member?", s_IntGroup_MemberP, 1);
	rb_define_alias(rb_cIntGroup, "include?", "member?");
	rb_define_method(rb_cIntGroup, "each", s_IntGroup_Each, 0);
	rb_define_method(rb_cIntGroup, "[]", s_IntGroup_ElementAtIndex, 1);
	rb_define_method(rb_cIntGroup, "add", s_IntGroup_Add, 1);
	rb_define_alias(rb_cIntGroup, "<<", "add");
	rb_define_method(rb_cIntGroup, "delete", s_IntGroup_Delete, 1);
	rb_define_method(rb_cIntGroup, "union", s_IntGroup_Union, 1);
	rb_define_method(rb_cIntGroup, "difference", s_IntGroup_Difference, 1);
	rb_define_method(rb_cIntGroup, "intersection", s_IntGroup_Intersection, 1);
	rb_define_method(rb_cIntGroup, "sym_difference", s_IntGroup_SymDifference, 1);
	rb_define_method(rb_cIntGroup, "convolute", s_IntGroup_Convolute, 1);
	rb_define_method(rb_cIntGroup, "deconvolute", s_IntGroup_Deconvolute, 1);
	rb_define_method(rb_cIntGroup, "offset", s_IntGroup_Offset, 1);
	rb_define_alias(rb_cIntGroup, "+", "union");
	rb_define_alias(rb_cIntGroup, "|", "union");
	rb_define_alias(rb_cIntGroup, "-", "difference");
	rb_define_alias(rb_cIntGroup, "&", "intersection");
	rb_define_alias(rb_cIntGroup, "^", "sym_difference");
	rb_define_method(rb_cIntGroup, "range_at", s_IntGroup_RangeAt, 1);
	rb_define_method(rb_cIntGroup, "merge", s_IntGroup_Merge, -1);
	rb_define_method(rb_cIntGroup, "subtract", s_IntGroup_Subtract, -1);
	rb_define_method(rb_cIntGroup, "inspect", s_IntGroup_Inspect, 0);
	rb_define_alias(rb_cIntGroup, "to_s", "inspect");
	rb_define_singleton_method(rb_cIntGroup, "[]", s_IntGroup_Create, -1);
	{
		VALUE igval = IntGroup_Alloc(rb_cIntGroup);
		IntGroup *ig;
		Data_Get_Struct(igval, IntGroup, ig);
		IntGroupAdd(ig, 0, ATOMS_MAX_NUMBER);
		rb_define_global_const("All", igval);
	}
}
