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

VALUE rb_cVector3D, rb_cTransform, rb_cLAMatrix, rb_cIntGroup;

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
	} else if (rb_obj_is_kind_of(val, rb_cLAMatrix)) {
		LAMatrix *mp;
		Data_Get_Struct(val, LAMatrix, mp);
		if ((mp->row == 1 && mp->column >= 3) || (mp->row >= 3 && mp->column == 1)) {
			vp->x = mp->data[0];
			vp->y = mp->data[1];
			vp->z = mp->data[2];
		}
	} else if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		/*  Vector3D[1] is rejected; this is desirable because Integer implements 
		    the index method ([]), which could cause confusion.  */
		rb_raise(rb_eMolbyError, "single number cannot be converted to a Vector3D");
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
 *     new
 *     new(Vector3D)
 *     new(Array)
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
 *     size -> Integer
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
 *     self[index]        -> Float
 *
 *  Element Reference---Returns the element at the given index. If the index is
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
 *     self[index] = val
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
 *     self == val   ->   bool
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
 *     self + val   ->   (new) Vector3D
 *
 *  Add two vectors element by element. Val is converted to Vector3D.
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
 *     self - val   ->   (new) Vector3D
 *
 *  Subtract two vectors element by element. Val is converted to Vector3D.
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
 *     self.dot(val) ->  Float
 *
 *  Calculate the dot (inner) product of the two vectors. Val is converted to Vector3D.
 *
 *  <b>See Also:</b> Vector3D.*
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
 *     self * numeric           ->  (new) Vector3D
 *     self * val    ->  Float
 *
 *  In the first form, the vector is scaled by the numeric. In the second
 *  form, the dot (inner) product of the two vectors are returned, which is equivalent to
 *  self.dot(val).
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
 *     self / numeric           ->  (new) Vector3D
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
 *     self.cross(val) ->  (new) Vector3D
 *
 *  Calculate the cross (outer) product of the two vectors. Val is converted to Vector3D.
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
 *     -self                ->  (new) Vector3D
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
 *     length                ->  Float
 *
 *  Calculate the Pythagorean length of the vector. 
 *  Note that this method is <em>not</em> an alias of Vector3D#size, which returns 3.
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
 *     length2                ->  Float
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
 *     normalize              ->  (new) Vector3D
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
 *     to_a                    ->  Array
 *
 *  Returns [self.x, self.y, self.z].
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
 *     each {|item| ...}
 *
 *  Calls block for x, y, z elements, passing that element as a parameter.
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
 *     x       -> Float
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
 *     y       -> Float
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
 *     z       -> Float
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
 *     x = val
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
 *     y = val
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
 *     z = val
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
 *     Vector3d[vx] -> (new) Vector3D.
 *     Vector3d[fx, fy, fz]   -> (new) Vector3D
 *
 *  Create a new vector3d object. The first form is equivalent to Vector3D#new(vx), and
 *  the second form is equivalent to Vector3D#new([fx, fy, fz]).
 */
static VALUE
s_Vector3D_Create(VALUE klass, VALUE args)
{
	VALUE val = s_Vector3D_Alloc(klass);
	if (RARRAY_LEN(args) == 1)
		s_Vector3D_Initialize(RARRAY_LEN(args), RARRAY_PTR(args), val);
	else s_Vector3D_Initialize(1, &args, val);
	return val;
}

/* 
 *  call-seq:
 *     inspect        -> String
 *
 *  Create a readable string like "Vector3D[fx, fy, fz]".
 */
static VALUE
s_Vector3D_Inspect(VALUE self)
{
	/*  self.class.name << self.to_a.inspect  */
/*	VALUE klass = CLASS_OF(self); */
/*	VALUE val = rb_funcall(klass, rb_intern("name"), 0); */
	VALUE val = rb_str_new2("Vector3D");
	self = s_Vector3D_ToArray(self);
	return rb_funcall(val, rb_intern("<<"), 1, rb_funcall(self, rb_intern("inspect"), 0));
}

#pragma mark ====== Transform Class ======

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
				(*tp)[i] = NUM2DBL(rb_Float(valp[i]));
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
					(*tp)[3 * i + j] = NUM2DBL(rb_Float(valp2[j]));
			}
			return;
		}
	array_format_error:
		rb_raise(rb_eMolbyError, "wrong array format; must be an array of either (1) four 3D column vectors or (2) 12 numerics");
	} else if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		/*  Transform[1] is rejected; this is desirable because Integer implements 
		 the index method ([]), which could cause confusion.  */
		rb_raise(rb_eMolbyError, "single number cannot be converted to a Transform");
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
					(*tp)[i * 3 + j] = NUM2DBL(rb_Float(rb_funcall(val, index_mid, 2, INT2FIX(i), INT2FIX(j))));
			}
		} else {
			/*  Other "array-like" object  */
			for (i = 0; i < 12; i++)
				(*tp)[i] = NUM2DBL(rb_Float(rb_funcall(val, index_mid, 1, INT2FIX(i))));
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
 *     new
 *     new(array)
 *     new(matrix)
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
 *  +matrix+ must respond to a method call <tt>matrix[col, row]</tt>
 *  where <tt>row</tt> is in <tt>0..2</tt> and <tt>col</tt> in <tt>0..3</tt>. 
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
 *     from_columns(c1, c2, c3, c4)
 *
 *  Returns a new Transform object built from four column vectors. The arguments
 *  <tt>c1..c4</tt> are vectors of (at least) three-dimension. This is equivalent
 *  to <tt>Transform.new([c1, c2, c3, c4])</tt>.
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
			tr[i * 3 + j] = w[j];
	}
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     from_rows(r1, r2, r3)
 *
 *  Returns a new Transform object built from three row vectors. The arguments
 *  <tt>r1, r2, r3</tt> are vectors of (at least) four-dimension.
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
			tr[j * 3 + i] = w[j];
	}
	return s_Transform_NewFromTransform(&tr);
}

/*
 *  call-seq:
 *     self[i, j]  -> Float
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
	w = (*tp)[n1 * 3 + n2];
	return rb_float_new(w);
}

/*
 *  call-seq:
 *     self[i, j] = val
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
	(*tp)[n1 * 3 + n2] = w;
	return rb_float_new(w);
}

/*
 *  call-seq:
 *     self == val  -> bool
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
 *     self + val  -> (new) Transform
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
 *     self - val  -> (new) Transform
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
 *     self * numeric          -> (new) Transform
 *     self * Vector3D         -> (new) Vector3D
 *     self * other_transform  -> (new) Transform
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
 *     identity  -> Transform
 *
 *  Returns an identity transform, <tt>[[1,0,0], [0,1,0], [0,0,1], [0,0,0]]</tt>.
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
 *     zero  -> Transform
 *
 *  Returns a zero transform, <tt>[[0,0,0], [0,0,0], [0,0,0], [0,0,0]]</tt>.
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
 *     diagonal(Array)
 *     diagonal(f1, f2 = nil, f3 = nil)
 *
 *  Returns a diagonal transform (the translational componets are all zero).
 *  In the first form, <tt>array[0], array[1], array[2]</tt> are for the
 *  x, y, z components, respectively. In the second form, <tt>f1, f2, f3</tt>
 *  are the x, y, z components. If <tt>f3</tt> is not given, the <tt>f2</tt>
 *  is used for the z components. If <tt>f2</tt> is not given, the <tt>f1</tt>
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
 *     inverse  -> (new) Transform
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
 *     self / val -> (new) Transform
 *
 *  Returns self * val.invert. If val is not a regular transform,
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
 *     transpose -> (new) Transform
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
 *     determinant -> Float
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
 *     trace -> Float
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
 *     column(index) -> Vector3D
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
	v.x = (*tp1)[n * 3];
	v.y = (*tp1)[n * 3 + 1];
	v.z = (*tp1)[n * 3 + 2];
	return ValueFromVector(&v);
}

/*
 *  call-seq:
 *     eigenvalues -> [[k1, k2, k3], v1, v2, v3]
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
 *     Transform[*args] -> (new) Transform
 *
 *  Create a new transform. Equivalent to Transform.new(args).
 *  If only one argument is given that is not a number, then Transform.new(args[0])
 *  is called instead of Transform.new(args).
 */
static VALUE
s_Transform_Create(VALUE klass, VALUE args)
{
	VALUE val = s_Transform_Alloc(klass);
	if (RARRAY_LEN(args) == 1 && !rb_obj_is_kind_of(RARRAY_PTR(args)[0], rb_cNumeric))
		s_Transform_Initialize(RARRAY_LEN(args), RARRAY_PTR(args), val);
	else
		s_Transform_Initialize(1, &args, val);
	return val;
}

/*
 *  call-seq:
 *     to_a  -> Array
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
		val[i] = rb_float_new((*tp1)[i]);
	return rb_ary_new4(12, val);
}

/*
 *  call-seq:
 *     inspect  -> String
 *
 *  Convert a transform to a string like 
 *  "Transform[[a11,a21,a31],[a12,a22,a32],[a13,a23,a33],[a14,a24,a34]]".
 */
static VALUE
s_Transform_Inspect(VALUE self)
{
	Transform *tp;
	int i, j;
/*	VALUE klass = CLASS_OF(self); */
/*	VALUE val = rb_funcall(klass, rb_intern("name"), 0); */
	VALUE val = rb_str_new2("Transform");
	ID mid = rb_intern("<<");
	ID mid2 = rb_intern("inspect");
	rb_str_cat(val, "[", 1);
	Data_Get_Struct(self, Transform, tp);
	for (i = 0; i < 4; i++) {
		rb_str_cat(val, "[", 1);
		for (j = 0; j < 3; j++) {
			double f;
			f = (*tp)[i * 3 + j];
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
 *     translation(vec) -> (new) Transform
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
 *     rotation(axis, angle, center = [0,0,0]) -> (new) Transform
 *
 *  Returns a transform corresponding to the rotation along the given axis and angle.
 *  Angle is given in degree. If center is also given, that point will be the center of rotation.
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
	ang = NUM2DBL(rb_Float(angle)) * kDeg2Rad;
	if (TransformForRotation(tr, &av, ang, &cv))
		rb_raise(rb_eMolbyError, "rotation axis cannot be a zero vector");
	return ValueFromTransform(&tr);
}

/*
 *  call-seq:
 *     reflection(axis, center = [0,0,0]) -> (new)Transform
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
 *     inversion(center = [0,0,0]) -> (new) Transform
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

#pragma mark ====== LAMatrix Class ======

/*  Coerce the value to LAMatrix if necessary */
LAMatrix *
LAMatrixFromValue(VALUE val, int *needsRelease, int rowHint, int columnHint)
{
	int i, j, row, column;
	LAMatrix *mp1;
	VALUE *valp1, *valp2, val2;
	static ID index_mid = 0, row_size_mid, column_size_mid, to_a_mid;
	if (rb_obj_is_kind_of(val, rb_cLAMatrix)) {
		if (needsRelease != NULL)
			*needsRelease = 0;
		return (LAMatrix *)DATA_PTR(val);
	} else if (rb_obj_is_kind_of(val, rb_cTransform)) {
		/*  Transform is described as a 4x4 matrix
		    a11 a12 a13 a14
		    a21 a22 a23 a24
		    a31 a32 a33 a34
			0   0   0   1  
		    or a 3x3 matrix, depending on the value of rowHint and columnHint.  */
		Transform *tp;
		Data_Get_Struct(val, Transform, tp);
		if (rowHint == 3 || columnHint == 3)
			row = 3;
		else row = 4;
		mp1 = LAMatrixNew(row, row);
		for (i = 0; i < row; i++) {
			for (j = 0; j < 3; j++) {
				mp1->data[i * row + j] = (*tp)[i * 3 + j];
			}
		}
		if (row == 4) {
			mp1->data[3] = mp1->data[7] = mp1->data[11] = 0.0;
			mp1->data[15] = 1.0;
		}
		if (needsRelease != NULL)
			*needsRelease = 1;
		return mp1;
	} else if (rb_obj_is_kind_of(val, rb_cVector3D)) {
		/*  Vector3D is described as a column vector (vx vy vz 1).transpose or (vx vy vz).transpose,
		    depending on the rowHint. */
		Vector *vp;
		Data_Get_Struct(val, Vector, vp);
		if (rowHint == 3)
			row = 3;
		else row = 4;
		mp1 = LAMatrixNew(row, 1);
		mp1->data[0] = vp->x;
		mp1->data[1] = vp->y;
		mp1->data[2] = vp->z;
		if (row == 4)
			mp1->data[3] = 1.0;
		if (needsRelease != NULL)
			*needsRelease = 1;
		return mp1;
	}
	if (index_mid == 0) {
		index_mid = rb_intern("[]");
		row_size_mid = rb_intern("row_size");
		column_size_mid = rb_intern("column_size");
		to_a_mid = rb_intern("to_a");
	}
	if (rb_respond_to(val, row_size_mid) && rb_respond_to(val, column_size_mid)) {
		/*  Matrix-type object  */
		column = NUM2INT(rb_Integer(rb_funcall(val, column_size_mid, 0)));
		row = NUM2INT(rb_Integer(rb_funcall(val, row_size_mid, 0)));
		if (column <= 0)
			rb_raise(rb_eMolbyError, "Bad column dimension (%d) for creating LAMatrix", column);
		if (row <= 0)
			rb_raise(rb_eMolbyError, "Bad row dimension (%d) for creating LAMatrix", row);
		mp1 = LAMatrixNew(row, column);
		for (i = 0; i < column; i++) {
			for (j = 0; j < row; j++)
				mp1->data[i * row + j] = NUM2DBL(rb_Float(rb_funcall(val, index_mid, 2, INT2FIX(i), INT2FIX(j))));
		}		
	} else {
		/*  Array  */
		if (TYPE(val) != T_ARRAY)
			val = rb_funcall(val, to_a_mid, 0);
		/*  Set of column vectors  */
		column = RARRAY_LEN(val);
		if (column == 0)
			rb_raise(rb_eMolbyError, "Cannot convert empty array to LAMatrix");
		valp1 = RARRAY_PTR(val);
		if (rb_obj_is_kind_of(valp1[0], rb_cNumeric)) {
			/*  A single column vector  */
			mp1 = LAMatrixNew(column, 1);  /* The first argument is actually the number of rows */
			for (i = 0; i < column; i++)
				mp1->data[i] = NUM2DBL(rb_Float(valp1[i]));
		} else {
			/*  Array of column vectors  */
			row = 0;
			valp2 = ALLOC_N(VALUE, column);
			for (i = 0; i < column; i++) {
				val2 = valp1[i];
				if (TYPE(val2) != T_ARRAY)
					val2 = rb_funcall(val2, to_a_mid, 0);
				if (RARRAY_LEN(val2) > row)
					row = RARRAY_LEN(val2);
				valp2[i] = val2;
			}
			if (row == 0)
				rb_raise(rb_eMolbyError, "Cannot convert array containing empty array to LAMatrix");
			mp1 = LAMatrixNew(row, column);
			for (i = 0; i < column; i++) {
				val2 = valp2[i];
				for (j = 0; j < RARRAY_LEN(val2); j++)
					mp1->data[i * row + j] = NUM2DBL(rb_Float((RARRAY_PTR(val2))[j]));
			}
			xfree(valp2);
		}
	}
	if (needsRelease != NULL)
		*needsRelease = 1;
	return mp1;
}

static VALUE
s_LAMatrix_Alloc(VALUE klass)
{
	return Data_Wrap_Struct(klass, 0, -1, NULL);
}

/*
 *  call-seq:
 *     new(column, row)
 *     new(array)
 *     new(matrix)
 *
 *  Returns a new LAMatrix object.
 *  In the first form, a zero LAMatrix of given size is returned. Note the order of the
 *  arguments; they are opposite to the mathematical convention.
 *  In the second form, the array must be either of an array (the column vector),
 *  or an array of arrays (a set of column vectors).
 *  In the third form, a new transform is built from a matrix. The argument
 *  +matrix+ must respond to a method call <tt>matrix[col, row]</tt>,
 *  <tt>row_size</tt> and <tt>column_size</tt>.
 */
static VALUE
s_LAMatrix_Initialize(int argc, VALUE *argv, VALUE self)
{
	LAMatrix *mp;
	if (argc == 2) {
		int row, column;
		row = NUM2INT(rb_Integer(argv[1]));
		column = NUM2INT(rb_Integer(argv[0]));
		if (column <= 0)
			rb_raise(rb_eMolbyError, "Bad column dimension (%d) for creating LAMatrix", column);
		if (row <= 0)
			rb_raise(rb_eMolbyError, "Bad row dimension (%d) for creating LAMatrix", row);
		mp = DATA_PTR(self);
		if (mp != NULL && mp->column * mp->row >= column * row) {
			mp->column = column;
			mp->row = row;
			memset(mp->data, 0, sizeof(mp->data[0]) * column * row);
		} else {
			if (mp != NULL)
				free(mp);
			mp = LAMatrixNew(row, column);
			DATA_PTR(self) = mp;
		}
	} else if (argc == 0 || argc >= 3) {
		rb_raise(rb_eArgError, "Wrong number of arguments: expecting two integers (row, column) or a single Array, LAMatrix, Vector3D, or Transform");
	} else {
		int needsRelease;
		mp = LAMatrixFromValue(argv[0], &needsRelease, 0, 0);
		if (!needsRelease) {
			/*  Needs duplicate  */
			mp = LAMatrixNewFromMatrix(mp);
		}
		if (DATA_PTR(self) != NULL)
			free(DATA_PTR(self));
		DATA_PTR(self) = mp;
	}
	return Qnil;
}

/*
 *  call-seq:
 *     initialize_copy(arg) -> self
 *
 *  Duplicate the matrix. Also used for implementation of dup.
 */
static VALUE
s_LAMatrix_InitializeCopy(VALUE self, VALUE val)
{
	LAMatrix *mp1;
	int needsRelease;
	mp1 = LAMatrixFromValue(val, &needsRelease, 0, 0);
	if (!needsRelease)
		mp1 = LAMatrixNewFromMatrix(mp1);
	if (DATA_PTR(self) != NULL)
		free(DATA_PTR(self));
	DATA_PTR(self) = mp1;
	return self;
}

/*
 *  call-seq:
 *     from_columns(r1,...)
 *
 *  Returns a new LAMatrix object built from column vectors.
 */
static VALUE
s_LAMatrix_NewFromColumns(VALUE klass, VALUE val)
{
	LAMatrix *mp;
	int needsRelease;
	mp = LAMatrixFromValue(val, &needsRelease, 0, 0);
	if (!needsRelease)
		mp = LAMatrixNewFromMatrix(mp);
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp);
}

/*
 *  call-seq:
 *     from_rows(r1,...)
 *
 *  Returns a new LAMatrix object built from row vectors.
 */
static VALUE
s_LAMatrix_NewFromRows(VALUE klass, VALUE val)
{
	LAMatrix *mp;
	int needsRelease;
	mp = LAMatrixFromValue(val, &needsRelease, 0, 0);
	if (!needsRelease)
		mp = LAMatrixNewFromMatrix(mp);
	LAMatrixTranspose(mp, mp);
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp);
}

/*
 *  call-seq:
 *     self[i, j]  -> Float
 *     self[i] -> Array
 *
 *  Get the element (+i+,+j+) of the matrix, i.e. column +i+, row +j+.
 *  In the second form, the i-th column vector is returned as an array of floating numbers.
 *  (If you want to get the column vector as a LAMatrix, use self.column(i) instead.)
 *  Be careful about the order of the arguments. It follows convention of multi-dimensional arrays
 *  rather than mathematical notation.
 */
static VALUE
s_LAMatrix_ElementAtIndex(int argc, VALUE *argv, VALUE self)
{
	LAMatrix *mp;
	double w;
	VALUE val1, val2;
	int n1, n2;
	rb_scan_args(argc, argv, "11", &val1, &val2);
	n1 = NUM2INT(val1);
	if (val2 != Qnil)
		n2 = NUM2INT(val2);
	else n2 = -1;
	Data_Get_Struct(self, LAMatrix, mp);
	if (n1 < 0 || n1 >= mp->column || (val2 != Qnil && (n2 < 0 || n2 >= mp->row)))
		rb_raise(rb_eMolbyError, "index to LAMatrix out of range");
	if (n2 >= 0) {
		w = mp->data[n1 * mp->row + n2];
		return rb_float_new(w);
	} else {
		val2 = rb_ary_new2(mp->row);
		for (n2 = 0; n2 < mp->row; n2++)
			rb_ary_store(val2, n2, rb_float_new(mp->data[n1 * mp->row + n2]));
		return val2;
	}
}

/*
 *  call-seq:
 *     self[i, j] = val
 *     self[i] = array
 *
 *  Set the element (+i+,+j+) of the matrix, i.e. column +i+, row +j+.
 *  In the second form, the i-th column vector is replaced with the given array. For convenience,
 *  if the given array is an array of array, then it is flattened before use.
 *  This allows use of a column vector or a row vector as the argument.
 *  Be careful about the order of the arguments. It follows convention of multi-dimensional arrays
 *  rather than mathematical notation.
 */
static VALUE
s_LAMatrix_SetElementAtIndex(int argc, VALUE *argv, VALUE self)
{
	LAMatrix *mp;
	double w;
	int n1, n2;
	VALUE idx1, idx2, val;
	rb_scan_args(argc, argv, "21", &idx1, &idx2, &val);
	if (argc == 2) {
		val = idx2;
		idx2 = Qnil;
	}
	n1 = NUM2INT(rb_Integer(idx1));
	if (idx2 != Qnil)
		n2 = NUM2INT(rb_Integer(idx2));
	else n2 = -1;
	Data_Get_Struct(self, LAMatrix, mp);
	if (n1 < 0 || n1 >= mp->column || (idx2 != Qnil && (n2 < 0 || n2 >= mp->row)))
		rb_raise(rb_eRangeError, "index to LAMatrix out of range");
	if (idx2 == Qnil) {
		val = rb_ary_to_ary(val);
		if (RARRAY_LEN(val) > 0 && rb_obj_is_kind_of(RARRAY_PTR(val)[0], rb_cArray)) {
			/*  Flatten  */
			if (RARRAY_LEN(val) == 1)
				val = RARRAY_PTR(val)[0];
			else {
				static ID sFlatten = 0;
				if (sFlatten == 0)
					sFlatten = rb_intern("flatten");
				val = rb_funcall(val, sFlatten, 0);
			}
		}
		if (RARRAY_LEN(val) != mp->row)
			rb_raise(rb_eArgError, "the row dimension does not match the given array");
		for (n2 = 0; n2 < mp->row; n2++) {
			mp->data[n1 * mp->row + n2] = NUM2DBL(rb_Float(RARRAY_PTR(val)[n2]));
		}
	} else {
		val = rb_Float(val);
		w = NUM2DBL(val);
		mp->data[n1 * mp->row + n2] = w;
	}
	return val;
}

/*
 *  call-seq:
 *     self == val  -> bool
 *
 *  Returns +true+ if and only if both dimensions and all the corresponding elements are equal.
 *  Usual caution about the comparison of floating-point numbers should be paid.
 */
static VALUE
s_LAMatrix_IsEqual(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2;
	int i, needsRelease, n;
	VALUE retval = Qtrue;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixFromValue(val, &needsRelease, mp1->row, mp1->column);
	if (mp1->column != mp2->column || mp1->row != mp2->row)
		retval = Qfalse;
	else {
		n = mp1->column * mp1->row;
		for (i = 0; i < n; i++) {
			if (mp1->data[i] != mp2->data[i]) {
				retval = Qfalse;
				break;
			}
		}
	}
	if (needsRelease)
		free(mp2);
	return retval;
}

/*
 *  call-seq:
 *     self + val  -> (new) LAMatrix
 */
static VALUE
s_LAMatrix_Add(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2, *mp3;
	int i, n, needsRelease;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixFromValue(val, &needsRelease, mp1->row, mp1->column);
	if (mp1->column != mp2->column && mp1->row != mp2->row)
		rb_raise(rb_eArgError, "mismatch dimensions of LAMatrix");
	mp3 = LAMatrixNew(mp1->row, mp1->column);
	n = mp1->row * mp1->column;
	for (i = 0; i < n; i++)
		mp3->data[i] = mp1->data[i] + mp2->data[i];
	if (needsRelease)
		free(mp2);
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp3);
}

/*
 *  call-seq:
 *     self - val  -> (new) LAMatrix
 */
static VALUE
s_LAMatrix_Subtract(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2, *mp3;
	int i, n, needsRelease;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixFromValue(val, &needsRelease, mp1->row, mp1->column);
	if (mp1->column != mp2->column && mp1->row != mp2->row)
		rb_raise(rb_eArgError, "mismatch dimensions of LAMatrix");
	mp3 = LAMatrixNew(mp1->row, mp1->column);
	n = mp1->row * mp1->column;
	for (i = 0; i < n; i++)
		mp3->data[i] = mp1->data[i] - mp2->data[i];
	if (needsRelease)
		free(mp2);
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp3);
}

/*
 *  call-seq:
 *     self * numeric          -> (new) LAMatrix
 *     self * other_transform  -> (new) LAMatrix
 *
 *  Perform the matrix multiplication. In the first form, a new matrix with scaled elements
 *  is returned. In the second form, the multiple of the two matrices is returned.
 */
static VALUE
s_LAMatrix_Multiply(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2, *mp3;
	int needsRelease, i, n;
	Data_Get_Struct(self, LAMatrix, mp1);
	if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		double w = NUM2DBL(rb_Float(val));
		mp2 = LAMatrixNewFromMatrix(mp1);
		n = mp2->column * mp2->row;
		for (i = 0; i < n; i++)
			mp2->data[i] *= w;
		return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2);
	} else {
		mp2 = LAMatrixFromValue(val, &needsRelease, mp1->column, 0);
		if (mp1->column != mp2->row)
			rb_raise(rb_eArgError, "mismatch dimensions in LAMatrix multiplication");
		mp3 = LAMatrixNew(mp1->row, mp2->column);
		LAMatrixMul(0, 0, 1.0, mp1, mp2, 0.0, mp3);
		if (needsRelease)
			free(mp2);
		return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp3);
	}
}

/*
 *  call-seq:
 *     self.add!(val) -> self
 *
 *  Add a matrix to self.
 */
static VALUE
s_LAMatrix_Add_Bang(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2;
	int i, n, needsRelease;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixFromValue(val, &needsRelease, mp1->row, mp1->column);
	if (mp1->column != mp2->column && mp1->row != mp2->row)
		rb_raise(rb_eArgError, "mismatch dimensions of LAMatrix");
	n = mp1->row * mp1->column;
	for (i = 0; i < n; i++)
		mp1->data[i] = mp1->data[i] + mp2->data[i];
	if (needsRelease)
		free(mp2);
	return self;
}

/*
 *  call-seq:
 *     self.sub!(val)  -> self
 *
 *  Subtract a matrix from self.
 */
static VALUE
s_LAMatrix_Subtract_Bang(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2;
	int i, n, needsRelease;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixFromValue(val, &needsRelease, mp1->row, mp1->column);
	if (mp1->column != mp2->column && mp1->row != mp2->row)
		rb_raise(rb_eArgError, "mismatch dimensions of LAMatrix");
	n = mp1->row * mp1->column;
	for (i = 0; i < n; i++)
		mp1->data[i] = mp1->data[i] - mp2->data[i];
	if (needsRelease)
		free(mp2);
	return self;
}

/*
 *  call-seq:
 *     self.multiply!(arg1, arg2, ...) -> self
 *
 *  Multiply the arguments to self.
 *  If argN is a string "t", then the subsequent matrix is
 *  transposed (the object is not modified; it is just transposed during calculation).
 *  If argN is a string "i", then the subsequent matrix is inverted (the object is not modified).
 *  The "t" and "i" can be combined as "ti" or "it".
 *  If argN is a number, then scalar multiplication is performed.
 *  Otherwise, argN must be a matrix, which has the same row-size as the column-size of the
 *  last result.
 */
static VALUE
s_LAMatrix_Multiply_Bang(int argc, VALUE *argv, VALUE self)
{
	LAMatrix *mp1, *mp2;
	__CLPK_doublereal mul;
	int needsRelease, n, inverseFlag, transposeFlag, temp1Flag, temp2Flag;
	if (self == Qnil) {
		/*  This is for handling the singleton method LAMatrix.multiply(arg1, arg2, ...)  */
		mp1 = NULL;
	} else {
		Data_Get_Struct(self, LAMatrix, mp1);
	}
	inverseFlag = transposeFlag = temp1Flag = 0;
	mul = 1.0;
	while (argc > 0) {
		VALUE val = argv[0];
		if (rb_obj_is_kind_of(val, rb_cNumeric)) {
			double w = NUM2DBL(rb_Float(val));
			mul *= w;
			if (argc == 1) {
				/*  This is the last argument: perform scalar multiplication  */
				if (!temp1Flag) {
					mp2 = LAMatrixAllocTempMatrix(mp1->row, mp1->column);
					memmove(mp2->data, mp1->data, sizeof(mp1->data[0]) * mp1->row * mp1->column);
					mp1 = mp2;
					temp1Flag = 1;
				}
				for (n = mp1->row * mp1->column - 1; n >= 0; n--) {
					mp1->data[n] *= mul;
				}
				mul = 1.0;  /*  Not necessary, but just to be consistent  */
			}
		} else if (rb_obj_is_kind_of(val, rb_cString)) {
			char *p = StringValuePtr(val);
			while (*p != 0) {
				if (*p == 'i' || *p == 'I')
					inverseFlag = !inverseFlag;
				else if (*p == 't' || *p == 'T')
					transposeFlag = !transposeFlag;
				else {
					if (temp1Flag)
						LAMatrixReleaseTempMatrix(mp1);
					rb_raise(rb_eArgError, "wrong character: only i or t is recognized");
				}
				p++;
			}
		} else {
			LAMatrix *mptemp;
			int sizeHint = (mp1 == NULL ? 0 : mp1->column);
			temp2Flag = 0;
			mp2 = LAMatrixFromValue(val, &needsRelease, (!transposeFlag ? sizeHint : 0), (transposeFlag ? sizeHint : 0));
			if (inverseFlag) {
				if (mp2->column != mp2->row)
					rb_raise(rb_eStandardError, "cannot invert non-square matrix");
				if (needsRelease)
					n = LAMatrixInvert(mp2, mp2);  /*  Invert in place  */
				else {
					mptemp = LAMatrixAllocTempMatrix(mp2->row, mp2->column);
					n = LAMatrixInvert(mptemp, mp2);
					mp2 = mptemp;   /*  Original mp2 is no longer needed  */
					temp2Flag = 1;  /*  mp2 should be later deallocated by LAMatrixReleaseTempMatrix  */
				}
				if (n != 0) {
					if (temp2Flag)
						LAMatrixReleaseTempMatrix(mp2);
					else if (needsRelease)
						LAMatrixRelease(mp2);
					if (temp1Flag)
						LAMatrixReleaseTempMatrix(mp1);
					rb_raise(rb_eStandardError, "cannot invert matrix");
				}
			}
			if (mp1 != NULL) {
				/*  Perform multiply  */
				if ((transposeFlag && (mp1->column != mp2->column)) || (!transposeFlag && (mp1->column != mp2->row)))
					rb_raise(rb_eArgError, "mismatch dimensions in LAMatrix multiplication");
				mptemp = LAMatrixAllocTempMatrix(mp1->row, (transposeFlag ? mp2->row : mp2->column));
				LAMatrixMul(0, transposeFlag, mul, mp1, mp2, 0.0, mptemp);
				if (temp1Flag)
					LAMatrixReleaseTempMatrix(mp1);
			} else {
				/*  mp1 = mp2 * mul  */
				int i;
				mptemp = LAMatrixAllocTempMatrix(mp2->row, mp2->column);
				if (transposeFlag) 
					LAMatrixTranspose(mptemp, mp2);
				else
					memmove(mptemp->data, mp2->data, sizeof(mp2->data[0]) * mp2->row * mp2->column);
				for (i = mp2->row * mp2->column - 1; i >= 0; i--)
					mptemp->data[i] *= mul;
			}
			mp1 = mptemp;
			temp1Flag = 1;
			if (temp2Flag)
				LAMatrixReleaseTempMatrix(mp2);
			else if (needsRelease)
				LAMatrixRelease(mp2);
			inverseFlag = transposeFlag = 0;
			mul = 1.0;
		}
		argc--;
		argv++;
	}
	if (self != Qnil) {
		/*  Copy the result back to self  */
		Data_Get_Struct(self, LAMatrix, mp2);
		if (mp1 != mp2)
			mp2 = LAMatrixResize(mp2, mp1->row, mp1->column);
		memmove(mp2->data, mp1->data, sizeof(mp1->data[0]) * mp1->row * mp1->column);
		DATA_PTR(self) = mp2;
		if (temp1Flag)
			LAMatrixReleaseTempMatrix(mp1);
		return self;
	} else {
		/*  Create a new LAMatrix from mp1  */
		mp2 = LAMatrixNewFromMatrix(mp1);
		if (temp1Flag)
			LAMatrixReleaseTempMatrix(mp1);
		return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2);
	}
}

/*
 *  call-seq:
 *     LAMatrix.multiply(arg1, arg2, ...) -> (new) LAMatrix
 *
 *  Returns a new LAMatrix that is the product of the arguments.
 *  The arguments are treated in the same way as LAMatrix.multiply!.
 */
static VALUE
s_LAMatrix_Multiply_Singleton(int argc, VALUE *argv, VALUE self)
{
	return s_LAMatrix_Multiply_Bang(argc, argv, Qnil);
}

/*
 *  call-seq:
 *     identity(size)  -> LAMatrix
 *
 *  Returns an identity matrix of size x size.
 */
static VALUE
s_LAMatrix_Identity(VALUE klass, VALUE val)
{
	LAMatrix *mp;
	int i, n;
	n = NUM2INT(rb_Integer(val));
	if (n <= 0)
		rb_raise(rb_eArgError, "invalid matrix dimension");
	mp = LAMatrixNew(n, n);
	for (i = 0; i < n; i++)
		mp->data[i * n + i] = 1.0;
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp);
}

/*
 *  call-seq:
 *     zero(column [, row])  -> LAMatrix
 *
 *  Returns a zero matrix of the specified size.
 */
static VALUE
s_LAMatrix_Zero(int argc, VALUE *argv, VALUE klass)
{
	LAMatrix *mp;
	VALUE val1, val2;
	int n1, n2;
	rb_scan_args(argc, argv, "11", &val1, &val2);
	n1 = NUM2INT(rb_Integer(val1));
	n2 = (val2 == Qnil ? n1 : NUM2INT(rb_Integer(val2)));
	if (n1 <= 0 || n2 <= 0)
		rb_raise(rb_eArgError, "invalid matrix dimension");
	mp = LAMatrixNew(n2, n1);
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp);
}

/*
 *  call-seq:
 *     diagonal(Array)
 *     diagonal(size, num)
 *
 *  Returns a diagonal matrix.
 *  In the first form, the dimension is defined by the size of the array and
 *  the numbers in the array become the diagonal elements.
 *  In the second form, a square matrix of size x size is created with 
 *  all diagonal elements equal to num.
 */
static VALUE
s_LAMatrix_Diagonal(int argc, VALUE *argv, VALUE klass)
{
	LAMatrix *mp;
	VALUE val1, val2, *valp;
	int i, n;
	rb_scan_args(argc, argv, "11", &val1, &val2);
	if (argc == 1) {
		val1 = rb_ary_to_ary(val1);
		n = RARRAY_LEN(val1);
		valp = RARRAY_PTR(val1);
		if (n == 0)
			rb_raise(rb_eArgError, "bad argument: empty array");
		mp = LAMatrixNew(n, n);
		for (i = 0; i < n; i++) {
			mp->data[i * n + i] = NUM2DBL(rb_Float(valp[i]));
		}
	} else {
		double w = NUM2DBL(rb_Float(val2));
		n = NUM2INT(rb_Integer(val1));
		mp = LAMatrixNew(n, n);
		for (i = 0; i < n; i++)
			mp->data[i * n + i] = w;
	}
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp);
}

/*
 *  call-seq:
 *     inverse  -> (new) LAMatrix
 *
 *  Returns the inverse matrix. If the matrix is not regular, an exception is raised.
 */
static VALUE
s_LAMatrix_Inverse(VALUE self)
{
	LAMatrix *mp1, *mp2;
	Data_Get_Struct(self, LAMatrix, mp1);
	if (mp1->column != mp1->row)
		rb_raise(rb_eArgError, "the matrix is not square");
	mp2 = LAMatrixNew(mp1->row, mp1->column);
	if (LAMatrixInvert(mp2, mp1))
		rb_raise(rb_eMolbyError, "singular matrix");
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2);
}

/*
 *  call-seq:
 *     inverse!  -> self
 *
 *  Replace self with the inverse matrix and returns self. 
 *  If the matrix is not regular, an exception is raised.
 */
static VALUE
s_LAMatrix_Inverse_Bang(VALUE self)
{
	LAMatrix *mp1;
	Data_Get_Struct(self, LAMatrix, mp1);
	if (mp1->column != mp1->row)
		rb_raise(rb_eArgError, "the matrix is not square");
	if (LAMatrixInvert(mp1, mp1))
		rb_raise(rb_eMolbyError, "singular matrix");
	return self;
}

/*
 *  call-seq:
 *     self / val -> (new) LAMatrix
 *
 *  Returns self * val.invert. If val is not a regular matrix,
 *  an exception is raised.
 */
static VALUE
s_LAMatrix_Divide(VALUE self, VALUE val)
{
	LAMatrix *mp1, *mp2, *mp3;
	int needsRelease, i, n;
	Data_Get_Struct(self, LAMatrix, mp1);
	if (rb_obj_is_kind_of(val, rb_cNumeric)) {
		double w = NUM2DBL(rb_Float(val));
		mp2 = LAMatrixNewFromMatrix(mp1);
		n = mp2->column * mp2->row;
		for (i = 0; i < n; i++)
			mp2->data[i] /= w;
		return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2);
	} else {
		mp2 = LAMatrixFromValue(val, &needsRelease, mp1->column, mp1->column);
		if (mp2->column != mp2->row)
			rb_raise(rb_eArgError, "cannot invert non-square matrix");
		if (mp1->column != mp2->row)
			rb_raise(rb_eArgError, "mismatch dimensions in LAMatrix multiplication");
		if (!needsRelease)
			mp2 = LAMatrixNewFromMatrix(mp2);
		if (LAMatrixInvert(mp2, mp2))
			rb_raise(rb_eArgError, "singular matrix");
		mp3 = LAMatrixNew(mp1->row, mp2->column);
		LAMatrixMul(0, 0, 1.0, mp1, mp2, 0.0, mp3);
		free(mp2);
		return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp3);
	}
}

/*
 *  call-seq:
 *     transpose -> (new) LAMatrix
 *
 *  Returns a new transposed matrix.
 */
static VALUE
s_LAMatrix_Transpose(VALUE self)
{
	LAMatrix *mp1, *mp2;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixNew(mp1->column, mp1->row);
	LAMatrixTranspose(mp2, mp1);
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2);
}

/*
 *  call-seq:
 *     transpose! -> self
 *
 *  Replace self with the transposed matrix and return self.
 */
static VALUE
s_LAMatrix_Transpose_Bang(VALUE self)
{
	LAMatrix *mp1;
	Data_Get_Struct(self, LAMatrix, mp1);
	LAMatrixTranspose(mp1, mp1);
	return self;
}

/*
 *  call-seq:
 *     determinant -> Float
 *
 *  Returns the determinant of the matrix.
 */
static VALUE
s_LAMatrix_Determinant(VALUE self)
{
	LAMatrix *mp1;
	Data_Get_Struct(self, LAMatrix, mp1);
	return rb_float_new(LAMatrixDeterminant(mp1));
}

/*
 *  call-seq:
 *     trace -> Float
 *
 *  Returns the trace (sum of the diagonal elements) of the matrix.
 */
static VALUE
s_LAMatrix_Trace(VALUE self)
{
	LAMatrix *mp1;
	Double tr;
	int i, n;
	Data_Get_Struct(self, LAMatrix, mp1);
	n = mp1->column;
	if (n > mp1->row)
		n = mp1->row;
	tr = 0.0;
	for (i = 0; i < n; i++)
		tr += mp1->data[i * mp1->row + i];
	return rb_float_new(tr);
}

/*
 *  call-seq:
 *     fnorm -> Float
 *
 *  Returns the Frobenius norm of the matrix.
 */
static VALUE
s_LAMatrix_Fnorm(VALUE self)
{
	LAMatrix *mp1;
	Double norm;
	int i, n;
	Data_Get_Struct(self, LAMatrix, mp1);
	norm = 0.0;
	n = mp1->row * mp1->column;
	for (i = 0; i < n; i++)
		norm += mp1->data[i] * mp1->data[i];
	norm = sqrt(norm);
	return rb_float_new(norm);
}

/*
 *  call-seq:
 *     fnorm2 -> Float
 *
 *  Returns the square of the Frobenius norm of the matrix.
 */
static VALUE
s_LAMatrix_Fnorm2(VALUE self)
{
	LAMatrix *mp1;
	Double norm;
	int i, n;
	Data_Get_Struct(self, LAMatrix, mp1);
	norm = 0.0;
	n = mp1->row * mp1->column;
	for (i = 0; i < n; i++)
		norm += mp1->data[i] * mp1->data[i];
	return rb_float_new(norm);
}

static VALUE
s_LAMatrix_Submatrix_sub(VALUE self, int rowpos, int columnpos, int row, int column)
{
	LAMatrix *mp1, *mp2;
	int i, j;
	Data_Get_Struct(self, LAMatrix, mp1);
	if (rowpos < 0 || rowpos >= mp1->row)
		rb_raise(rb_eArgError, "row number out of range");
	if (columnpos < 0 || columnpos >= mp1->column)
		rb_raise(rb_eArgError, "column number out of range");
	if (row == -1)
		row = mp1->row - rowpos;
	else if (row <= 0 || rowpos + row > mp1->row)
		rb_raise(rb_eArgError, "number of rows out of range");
	if (column == -1)
		column = mp1->column - columnpos;
	else if (column <= 0 || columnpos + column > mp1->column)
		rb_raise(rb_eArgError, "number of columns out of range");
	mp2 = LAMatrixNew(row, column);
	for (i = 0; i < row; i++) {
		for (j = 0; j < column; j++) {
			mp2->data[j * row + i] = mp1->data[(j + columnpos) * mp1->row + i + rowpos];
		}
	}
	return Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2);
}

/*
 *  call-seq:
 *     submatrix(columnpos, rowpos, column, row) -> LAMatrix
 *
 *  Returns the submatrix beginning from (columnpos, rowpos) and size (column, row).
 *  Note the order of the arguments; they are opposite to the mathematical convention.
 *  If -1 is specified for row or column, all elements in that direction are used.
 */
static VALUE
s_LAMatrix_Submatrix(VALUE self, VALUE columnposval, VALUE rowposval, VALUE columnval, VALUE rowval)
{
	return s_LAMatrix_Submatrix_sub(self, NUM2INT(rb_Integer(rowposval)), NUM2INT(rb_Integer(columnposval)), NUM2INT(rb_Integer(rowval)), NUM2INT(rb_Integer(columnval)));
}

/*
 *  call-seq:
 *     column(index) -> LAMatrix
 *
 *  Returns the index-th column vector as a (N, 1) matrix.
 */
static VALUE
s_LAMatrix_Column(VALUE self, VALUE val)
{
	return s_LAMatrix_Submatrix_sub(self, 0, NUM2INT(rb_Integer(val)), -1, 1);
}

/*
 *  call-seq:
 *     row(index) -> LAMatrix
 *
 *  Returns the index-th row vector as a (1, N) matrix.
 */
static VALUE
s_LAMatrix_Row(VALUE self, VALUE val)
{
	return s_LAMatrix_Submatrix_sub(self, NUM2INT(rb_Integer(val)), 0, 1, -1);
}

/*
 *  call-seq:
 *     column_size -> LAMatrix
 *
 *  Returns the column size.
 */
static VALUE
s_LAMatrix_ColumnSize(VALUE self)
{
	LAMatrix *mp;
	Data_Get_Struct(self, LAMatrix, mp);
	return INT2NUM(mp->column);
}

/*
 *  call-seq:
 *     row_size -> LAMatrix
 *
 *  Returns the row size.
 */
static VALUE
s_LAMatrix_RowSize(VALUE self)
{
	LAMatrix *mp;
	Data_Get_Struct(self, LAMatrix, mp);
	return INT2NUM(mp->row);
}

/*
 *  call-seq:
 *     eigenvalues -> [eigenvalues, eigenvectors]
 *
 *  Calculate the eigenvalues and eigenvectors. The matrix must be symmetric.
 */
static VALUE
s_LAMatrix_Eigenvalues(VALUE self)
{
	LAMatrix *mp1, *mp2, *mp3;
	int info;
	Data_Get_Struct(self, LAMatrix, mp1);
	if (mp1->column != mp1->row)
		rb_raise(rb_eArgError, "cannot get eigenvectors for non-square matrix");
	mp2 = LAMatrixNew(mp1->row, 1);
	mp3 = LAMatrixNew(mp1->row, mp1->column);
	if ((info = LAMatrixSymDiagonalize(mp2, mp3, mp1)) != 0)
		rb_raise(rb_eArgError, "cannot diagonalize");
	return rb_ary_new3(2, Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2), Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp3));
}

/*
 *  call-seq:
 *     svd -> [left_matrix, singular_values, right_matrix]
 *
 *  Decompose the given (m,n) matrix to a product of three matrices, U, S, V.
 *  U is a (m,m) orthogonal matrix, S is a (m,n) matrix which is zero except for
 *  min(m,n) diagonal elements, and V is a (n,n) orthogonal matrix.
 *  (Usually SVD is defined as M = U*S*transpose(V), but this methods returns
 *  transpose(V) rather than V.) The singular_values is a min(m,n) dimension
 *  column vector.
 */
static VALUE
s_LAMatrix_SVD(VALUE self)
{
	LAMatrix *mp1, *mp2, *mp3, *mp4;
	int info, n;
	Data_Get_Struct(self, LAMatrix, mp1);
	mp2 = LAMatrixNew(mp1->row, mp1->row);
	n = (mp1->column > mp1->row ? mp1->row : mp1->column);
	mp3 = LAMatrixNew(n, 1);
	mp4 = LAMatrixNew(mp1->column, mp1->column);
	if ((info = LAMatrixSingularValueDecomposition(mp2, mp3, mp4, mp1)) != 0)
		rb_raise(rb_eArgError, "cannot perform singular value decomposition");
	return rb_ary_new3(3, Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp2), Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp3), Data_Wrap_Struct(rb_cLAMatrix, 0, -1, mp4));
}

/*
 *  call-seq:
 *     LAMatrix[f1, f2, ..., fn] -> (new) LAMatrix
 *     LAMatrix[a1, a2, ..., an] -> (new) LAMatrix
 *     LAMatrix[obj] -> (new) LAMatrix
 *
 *  Create a new matrix.
 *  In the first form, f1...fn must be numbers, and a (n, 1) matrix (a column vector; LAMatrix[1, n] in our convention) is created.
 *  In the second form, a1...an must be array-like objects, and a (m, n) matrix (m is the maximum 
 *  dimension of a1...an; LAMatrix[n, m] in our convention) is created.
 *  In the third form, obj must be either an array, a Vector3D, a Transform, or an LAMatrix.
 */
static VALUE
s_LAMatrix_Create(VALUE klass, VALUE args)
{
	VALUE val = s_LAMatrix_Alloc(klass);
	int argc = RARRAY_LEN(args);
	VALUE *argv = RARRAY_PTR(args);
	if (argc == 0)
		rb_raise(rb_eArgError, "the argument should not be empty");
	if (rb_obj_is_kind_of(argv[0], rb_cNumeric) || rb_obj_is_kind_of(argv[0], rb_cArray)) {
		/*  First form and second form  */
		s_LAMatrix_Initialize(1, &args, val);
	} else {
		/*  Third form: type check is done within initialize  */
		s_LAMatrix_Initialize(argc, argv, val);
	}
	return val;
}

/*
 *  call-seq:
 *     to_a  -> Array
 *
 *  Convert a matrix to a nested array (array of column vectors).
 */
static VALUE
s_LAMatrix_ToArray(VALUE self)
{
	LAMatrix *mp1;
	VALUE val;
	int i, j;
	Data_Get_Struct(self, LAMatrix, mp1);
	val = rb_ary_new2(mp1->column);
	for (i = 0; i < mp1->column; i++) {
		VALUE val2 = rb_ary_new2(mp1->row);
		for (j = 0; j < mp1->row; j++) {
			rb_ary_store(val2, j, rb_float_new(mp1->data[i * mp1->row + j]));
		}
		rb_ary_store(val, i, val2);
	}
	return val;
}

/*
 *  call-seq:
 *     inspect  -> String
 *
 *  Convert a matrix to a string like 
 *  "LAMatrix[[a11,a21,a31],[a12,a22,a32],[a13,a23,a33],[a14,a24,a34]]".
 */
static VALUE
s_LAMatrix_Inspect(VALUE self)
{
	LAMatrix *mp;
	int i, j;
	/*	VALUE klass = CLASS_OF(self); */
	/*	VALUE val = rb_funcall(klass, rb_intern("name"), 0); */
	VALUE val = rb_str_new2("LAMatrix");
	ID mid = rb_intern("<<");
	ID mid2 = rb_intern("inspect");
	rb_str_cat(val, "[", 1);
	Data_Get_Struct(self, LAMatrix, mp);
	for (i = 0; i < mp->column; i++) {
		rb_str_cat(val, "[", 1);
		for (j = 0; j < mp->row; j++) {
			double f;
			f = mp->data[i * mp->row + j];
			rb_funcall(val, mid, 1, rb_funcall(rb_float_new(f), mid2, 0));
			if (j < mp->row - 1)
				rb_str_cat(val, ",", 1);
		}
		rb_str_cat(val, "]", 1);
		if (i < mp->column - 1)
			rb_str_cat(val, ",", 1);
	}
	rb_str_cat(val, "]", 1);
	return val;
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

/*
 *  call-seq:
 *     new(arg1, arg2,...)
 *     new(arg1, arg2,...) {|i| ...}
 *
 *  Create a new integer group. If no arguments are given, an empty group is returned.
 *  The arguments are either IntGroup, Range, Enumerable, or Numeric. In either case,
 *  the non-negative integers included in the arguments are added to the result.
 *  If a block is given, the block is called with the each integer in the given arguments,
 *  and the integer group consisting with the returned integers is returned.
 */
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

/*
 *  call-seq:
 *     clear -> self
 *
 *  Discard all integers included in self.
 */
static VALUE
s_IntGroup_Clear(VALUE self)
{
	IntGroup *ig;
	Data_Get_Struct(self, IntGroup, ig);
	IntGroupClear(ig);
	return self;
}

/*
 *  call-seq:
 *     dup(IntGroup) -> (new) IntGroup
 *
 *  (Deep) copy the given IntGroup.
 */
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

/*
 *  call-seq:
 *     length -> Integer
 *     size -> Integer
 *
 *  Returns the number of integers included in self.
 */
static VALUE
s_IntGroup_Length(VALUE self)
{
	IntGroup *ig;
	Data_Get_Struct(self, IntGroup, ig);
	return INT2NUM(IntGroupGetCount(ig));
}

/*
 *  call-seq:
 *     member?(val) -> bool
 *     include?(val) -> bool
 *
 *  Check whether the val is included in self.
 */
static VALUE
s_IntGroup_MemberP(VALUE self, VALUE val)
{
	IntGroup *ig;
	int n = NUM2INT(val);
	Data_Get_Struct(self, IntGroup, ig);
	return (IntGroupLookup(ig, n, NULL) ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     self[index] -> Integer or nil
 *
 *  Get the index-th point in self. If the index is out of range, nil is returned.
 */
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

/*
 *  call-seq:
 *     each {|i| ...}
 *
 *  Call the block with each integer in self.
 */
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

/*
 *  call-seq:
 *     add(IntGroup) -> self
 *
 *  Add the points in the given group.
 */
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

/*
 *  call-seq:
 *     delete(IntGroup) -> self
 *
 *  Remove the points in the given group.
 */
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

/*
 *  call-seq:
 *     union(val) -> (new)IntGroup
 *     self + val -> (new)IntGroup
 *     self | val -> (new)IntGroup
 *
 *  Returns a union group.
 */
static VALUE
s_IntGroup_Union(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupUnion);
}

/*
 *  call-seq:
 *     intersection(val) -> (new)IntGroup
 *     self & val -> (new)IntGroup
 *
 *  Returns an intersection group.
 */
static VALUE
s_IntGroup_Intersection(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupIntersect);
}

/*
 *  call-seq:
 *     difference(val) -> (new)IntGroup
 *     self - val -> (new)IntGroup
 *
 *  Returns a difference group.
 */
static VALUE
s_IntGroup_Difference(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupDifference);
}

/*
 *  call-seq:
 *     sym_difference(val) -> (new)IntGroup
 *     self ^ val -> (new)IntGroup
 *
 *  Returns a symmetric-difference group (i.e. a group containing elements that are included
 *  in either self or val but not both).
 */
static VALUE
s_IntGroup_SymDifference(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupXor);
}

/*
 *  call-seq:
 *     convolute(val) -> (new)IntGroup
 *
 *  For each element n in self, get the n-th point in val, and return the result as a new group.
 *  If n is out of range, then that point is ignored.
 *
 *  <b>See Also:</b> IntGroup#deconvolute. If all points in self are within the range of val, then 
 *  an equation <tt>self.convolute(val).deconvolute(val) == self</tt> holds.
 */
static VALUE
s_IntGroup_Convolute(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupConvolute);
}

/*
 *  call-seq:
 *     deconvolute(val) -> (new)IntGroup
 *
 *  For each element n in self, find the point n in val, and return the found indices as a new group.
 *  If n is not found in val, then that point is ignored.
 *
 *  <b>See Also:</b> IntGroup#convolute. If all points in self are found in val, then 
 *  an equation <tt>self.deconvolute(val).convolute(val) == self</tt> holds.
 */
static VALUE
s_IntGroup_Deconvolute(VALUE self, VALUE val)
{
	return s_IntGroup_Binary(self, val, IntGroupDeconvolute);
}

/*
 *  call-seq:
 *     range_at(val) -> Range
 *
 *  Split self into consecutive chunks of integers, and return the val-th chunk as a Range.
 *  This method is relatively efficient, because it directly uses the internal representation
 *  of IntGroup.
 */
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

/*
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
*/

/*
 *  call-seq:
 *     offset(val) -> self
 *
 *  Move all points by an integer value. A negative val is allowed, but it
 *  must be no smaller than -(self[0]), otherwise an exception is thrown.
 */
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

/*
 *  call-seq:
 *     inspect -> String
 *
 *  Create a String in the form "IntGroup[...]".
 */
static VALUE
s_IntGroup_Inspect(VALUE self)
{
	int i, sp, ep;
	IntGroup *ig;
	char buf[64];
/*	VALUE klass = CLASS_OF(self);
	VALUE val = rb_funcall(klass, rb_intern("name"), 0); */
/*	rb_str_cat(val, "[", 1); */
	VALUE val = rb_str_new2("IntGroup[");
	Data_Get_Struct(self, IntGroup, ig);
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
	rb_cVector3D = rb_define_class_under(rb_mMolby, "Vector3D", rb_cObject);
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
	rb_cTransform = rb_define_class_under(rb_mMolby, "Transform", rb_cObject);
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

	/*  class LAMatrix  */
	rb_cLAMatrix = rb_define_class_under(rb_mMolby, "LAMatrix", rb_cObject);
	rb_define_alloc_func(rb_cLAMatrix, s_LAMatrix_Alloc);
	rb_define_method(rb_cLAMatrix, "initialize", s_LAMatrix_Initialize, -1);
	rb_define_private_method(rb_cLAMatrix, "initialize_copy", s_LAMatrix_InitializeCopy, 1);
	rb_define_method(rb_cLAMatrix, "[]", s_LAMatrix_ElementAtIndex, -1);
	rb_define_method(rb_cLAMatrix, "[]=", s_LAMatrix_SetElementAtIndex, -1);
	rb_define_method(rb_cLAMatrix, "==", s_LAMatrix_IsEqual, 1);
	rb_define_method(rb_cLAMatrix, "+", s_LAMatrix_Add, 1);
	rb_define_method(rb_cLAMatrix, "-", s_LAMatrix_Subtract, 1);
	rb_define_method(rb_cLAMatrix, "*", s_LAMatrix_Multiply, 1);
	rb_define_method(rb_cLAMatrix, "/", s_LAMatrix_Divide, 1);
	rb_define_method(rb_cLAMatrix, "add!", s_LAMatrix_Add_Bang, 1);
	rb_define_method(rb_cLAMatrix, "sub!", s_LAMatrix_Subtract_Bang, 1);
	rb_define_method(rb_cLAMatrix, "multiply!", s_LAMatrix_Multiply_Bang, -1);
	rb_define_method(rb_cLAMatrix, "inverse", s_LAMatrix_Inverse, 0);
	rb_define_method(rb_cLAMatrix, "inverse!", s_LAMatrix_Inverse_Bang, 0);
	rb_define_method(rb_cLAMatrix, "transpose", s_LAMatrix_Transpose, 0);
	rb_define_method(rb_cLAMatrix, "transpose!", s_LAMatrix_Transpose_Bang, 0);
	rb_define_method(rb_cLAMatrix, "determinant", s_LAMatrix_Determinant, 0);
	rb_define_method(rb_cLAMatrix, "trace", s_LAMatrix_Trace, 0);
	rb_define_method(rb_cLAMatrix, "fnorm", s_LAMatrix_Fnorm, 0);
	rb_define_method(rb_cLAMatrix, "fnorm2", s_LAMatrix_Fnorm2, 0);
	rb_define_method(rb_cLAMatrix, "submatrix", s_LAMatrix_Submatrix, 4);
	rb_define_method(rb_cLAMatrix, "column", s_LAMatrix_Column, 1);
	rb_define_method(rb_cLAMatrix, "row", s_LAMatrix_Row, 1);
	rb_define_method(rb_cLAMatrix, "column_size", s_LAMatrix_ColumnSize, 0);
	rb_define_method(rb_cLAMatrix, "row_size", s_LAMatrix_RowSize, 0);
	rb_define_method(rb_cLAMatrix, "eigenvalues", s_LAMatrix_Eigenvalues, 0);
	rb_define_method(rb_cLAMatrix, "svd", s_LAMatrix_SVD, 0);
	rb_define_method(rb_cLAMatrix, "to_a", s_LAMatrix_ToArray, 0);
	rb_define_method(rb_cLAMatrix, "inspect", s_LAMatrix_Inspect, 0);
	rb_define_alias(rb_cLAMatrix, "to_s", "inspect");
	rb_define_singleton_method(rb_cLAMatrix, "diagonal", s_LAMatrix_Diagonal, -1);
	rb_define_singleton_method(rb_cLAMatrix, "[]", s_LAMatrix_Create, -2);
	rb_define_singleton_method(rb_cLAMatrix, "from_columns", s_LAMatrix_NewFromColumns, -2);
	rb_define_singleton_method(rb_cLAMatrix, "from_rows", s_LAMatrix_NewFromRows, -2);
	rb_define_singleton_method(rb_cLAMatrix, "identity", s_LAMatrix_Identity, 1);
	rb_define_singleton_method(rb_cLAMatrix, "zero", s_LAMatrix_Zero, -1);
	rb_define_singleton_method(rb_cLAMatrix, "multiply", s_LAMatrix_Multiply_Singleton, -1);

	/*  class IntGroup  */
	rb_cIntGroup = rb_define_class_under(rb_mMolby, "IntGroup", rb_cObject);
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
/*	rb_define_method(rb_cIntGroup, "merge", s_IntGroup_Merge, -1);
	rb_define_method(rb_cIntGroup, "subtract", s_IntGroup_Subtract, -1); */
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
