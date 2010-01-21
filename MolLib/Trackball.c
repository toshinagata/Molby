/*
 *  Trackball.c
 *
 *  Created by Toshi Nagata on 2005/09/17.
 *  Copyright 2005-2008 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "Trackball.h"
#include "Types.h"
#include <string.h>

#pragma mark ====== Internal functions ======

static void
multiplyQuat(const float *a, const float *b, float *c)
{
    c[0] = a[0] * b[0] - (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
    c[1] = b[0] * a[1] + a[0] * b[1] + (a[2] * b[3] - a[3] * b[2]);
    c[2] = b[0] * a[2] + a[0] * b[2] + (a[3] * b[1] - a[1] * b[3]);
    c[3] = b[0] * a[3] + a[0] * b[3] + (a[1] * b[2] - a[2] * b[1]);
}

static void
normalizeQuat(float *a)
{
    float f;
    f = 1.0 / sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
    a[0] *= f;
    a[1] *= f;
    a[2] *= f;
    a[3] *= f;
}

static void
getRotationToCenter(float x, float y, float radius, float *quat)
{
    const float sqrt2 = 1.41421356 / 2;
    float r, r2;
    r2 = radius * radius;
    r = x * x + y * y;
    if (r >= r2) {
        /*  Rotation of 90 degree along (y, -x, 0) */
        r = 1.0 / sqrt(r);
        x *= r;
        y *= r;
        quat[0] = sqrt2;
        quat[1] = sqrt2 * y;
        quat[2] = -sqrt2 * x;
        quat[3] = 0.0;
    } else {
        float z, cosAng, sinAng;
        z = sqrt(r2 - r);	/* z coordinate of the point on the sphere */
        /*  Rotation of arccos(z/radius) along (y, -x, 0)  */
        if (r < 1e-6) {
            quat[0] = 1.0;
            quat[1] = quat[2] = quat[3] = 0.0;
        } else {
            cosAng = sqrt(0.5 * (1.0 + z / radius));
            sinAng = sqrt(0.5 * (1.0 - z / radius));
            quat[0] = cosAng;
            r = 1.0 / sqrt(r);
            quat[1] = sinAng * y * r;
            quat[2] = -sinAng * x * r;
            quat[3] = 0.0;
        }
    }
}

static void
rotation2Quat(const float *A, float *q)
{
    float ang2;  /* The half-angle */
    float sinAng2; /* sin(half-angle) */
    
    /* Convert a GL-style rotation to a quaternion.  The GL rotation looks like this:
     * {angle, x, y, z}, the corresponding quaternion looks like this:
     * {{v}, cos(angle/2)}, where {v} is {x, y, z} / sin(angle/2).  */
    
    ang2 = A[0] * kDeg2Rad * 0.5;  /* Convert from degrees to radians, get the half-angle. */
    sinAng2 = sin(ang2);
    q[0] = A[1] * sinAng2; q[1] = A[2] * sinAng2; q[2] = A[3] * sinAng2;
    q[3] = cos(ang2);
}

#pragma mark ====== New/Retain/Release ======

Trackball *
TrackballNew(void)
{
	Trackball *track = (Trackball *)calloc(sizeof(Trackball), 1);
	MALLOC_CHECK(track, "allocating a trackball record");
    track->quat[0] = track->tempQuat[0] = 1.0;
	return track;
}

void
TrackballRetain(Trackball *track)
{
	if (track != NULL)
		track->refCount++;
}

void
TrackballRelease(Trackball *track)
{
	if (track != NULL && --(track->refCount) == 0) {
		free(track);
	}
}


/*  Get the current rotation (in GL format)  */
//- (void)getRotate:(float *)a
void
TrackballGetRotate(const Trackball *track, float *a)
{
    float w[4], t, f;
	NULL_CHECK(track, "TrackballGetRotate");	
    /*  Rotate: multiply the two quaternions and convert it to a GL rotater  */
    multiplyQuat(track->tempQuat, track->quat, w);
    normalizeQuat(w);
    if (fabs(fabs(w[0]) - 1.0) < 1e-7) {
        /*  Identity rotation  */
        a[0] = 0.0;
        a[1] = 1.0;
        a[2] = a[3] = 0.0;
    } else {
        /*  Quaternion: (cos(theta/2), sin(theta/2)*(x,y,z))  */
        /*  GL rotater: (theta, x, y, z)  */
        t = acos(w[0]);
        a[0] = t * 2.0 * kRad2Deg;
        f = 1.0 / sin(t);
        a[1] = w[1] * f;
        a[2] = w[2] * f;
        a[3] = w[3] * f;
    }
}

/*  Get the current translation.  The result values should be multiplied with
 * the dimension of the model.  */
void
TrackballGetTranslate(const Trackball *track, float *a)
//- (void)getTranslate:(float *)a
{
	NULL_CHECK(track, "TrackballGetTranslate");
    a[0] = track->tempTrans[0] + track->trans[0];
    a[1] = track->tempTrans[1] + track->trans[1];
    a[2] = track->tempTrans[2] + track->trans[2];
/*    f = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (f > 1.0) {
        f = 1.0 / sqrt(f);
        a[0] *= f;
        a[1] *= f;
        a[2] *= f;
    } */
}

/*  Get the perspective parameter (fovy and distance)
 *  Fovy (a[0]) is in degree.
 *  Distance (a[1]) should be multiplied with the dimension of the model. */
void
TrackballGetPerspective(const Trackball *track, float *a)
//- (void)getPerspective:(float *)a
{
    float f;
	NULL_CHECK(track, "TrackballGetPerspective");
    /*  Get the current scale factor (-1.0 to 1.0 for x1/10 to x10)  */
    f = track->tempScale + track->scale;
    if (f < 0) {
        /*  Shrink: fovy is fixed at 30 degree, and distance is enlarged  */
        if (f < -5.0)
            f = -5.0;
        a[1] = kCot15Deg * exp(-kLog10 * f);
        a[0] = 30.0;
    } else {
        /*  Expand: fovy is reduced while distance fixed at kCot15Deg  */
        if (f > 5.0)
            f = 5.0;
        a[0] = 30.0 * exp(-kLog10 * f);
        a[1] = kCot15Deg;
    }
}

void
TrackballReset(Trackball *track)
{
	NULL_CHECK(track, "TrackballReset");
	memset(track->trans, 0, sizeof(track->trans));
	memset(track->tempTrans, 0, sizeof(track->trans));
	memset(track->quat, 0, sizeof(track->trans));
	memset(track->tempQuat, 0, sizeof(track->trans));
    track->quat[0] = track->tempQuat[0] = 1.0;
	track->tempScale = 0;
	track->scale = 0;
}

void
TrackballSetScale(Trackball *track, float scale)
{
	NULL_CHECK(track, "TrackballSetScale");
	track->tempScale = 0;
	track->scale = scale;
}

void
TrackballSetTranslate(Trackball *track, const float *a)
{
	NULL_CHECK(track, "TrackballSetTranslate");
	track->tempTrans[0] = track->tempTrans[1] = track->tempTrans[2] = 0;
	track->trans[0] = a[0];
	track->trans[1] = a[1];
	track->trans[2] = a[2];
}

#pragma mark ====== Mouse operations ======

void
TrackballStartDragging(Trackball *track, const float *mousePos, TrackballMode mode)
//- (void)startAt:(NSPoint)pt sender:(NSView *)sender mode:(int)inMode
{
	NULL_CHECK(track, "TrackballStartDragging");
    
    /*  1: rotate, 2: translate, 3: scale  */
    track->mode = mode;

    /*  The start point vector  */
	if (mousePos != NULL) {
		track->start[0] = mousePos[0];
		track->start[1] = mousePos[1];
	} else {
		track->start[0] = track->start[1] = 0.0;
	}
    
    /*  The rotation from the start point to the center, for trackball operation  */
    getRotationToCenter(track->start[0], track->start[1], 1.0, track->startQuat);
}

void
TrackballSetTemporaryRotation(Trackball *track, const float *q)
{
	memmove(track->tempQuat, q, sizeof(float) * 4);
}

void
TrackballDrag(Trackball *track, const float *mousePos)
//- (void)dragTo:(NSPoint)pt sender:(NSView *)sender
{
    float rot[4];
    float w1[4], w2[4];
    float w;
	NULL_CHECK(track, "TrackballDrag");

    if (track->mode == kTrackballRotateMode) {

        /*  Rotation from the center to the end point  */
        getRotationToCenter(mousePos[0], mousePos[1], 1.0, w1);
        w1[1] = -w1[1];
        w1[2] = -w1[2];
        w1[3] = -w1[3];
        
        //  Accumulate 'start to center' and 'center to end' in this order
        multiplyQuat(track->startQuat, w1, track->tempQuat);
    
    } else if (track->mode == kTrackballTranslateMode) {

        /*  Make an (x, y, 0) vector in viewport coordinate  */
        TrackballGetPerspective(track, rot);
        w = (rot[1] * tan(rot[0] * 0.5 * kDeg2Rad) * 2.0);
        rot[0] = 0.0;
        rot[1] = (mousePos[0] - track->start[0]) * w;
        rot[2] = (mousePos[1] - track->start[1]) * w;
        rot[3] = 0.0;
        
        /*  Rotate with the viewport transform to give world coordinate
         *  y = q* . x . q, where x = original vector, y = transformed vector,
         *  q = quaternion, q* = conjugate quaternion  */
        multiplyQuat(rot, track->quat, w1);
        w2[0] = track->quat[0];
        w2[1] = -track->quat[1];
        w2[2] = -track->quat[2];
        w2[3] = -track->quat[3];
        multiplyQuat(w2, w1, rot);
        track->tempTrans[0] = rot[1];
        track->tempTrans[1] = rot[2];
        track->tempTrans[2] = rot[3];
    } else if (track->mode == kTrackballScaleMode) {
        w = (mousePos[0] - track->start[0]);
        if (w > 5.0)
            w = 5.0;
        else if (w < -5.0)
            w = -5.0;
        track->tempScale = w;
    }
}

void
TrackballEndDragging(Trackball *track, const float *mousePos)
//- (void)endAt:(NSPoint)pt sender:(NSView *)sender
{
    float w[4];
	NULL_CHECK(track, "TrackballEndDragging");
	if (mousePos != NULL)
		TrackballDrag(track, mousePos);
    multiplyQuat(track->tempQuat, track->quat, w);
    normalizeQuat(w);
    track->quat[0] = w[0];
    track->quat[1] = w[1];
    track->quat[2] = w[2];
    track->quat[3] = w[3];
    track->tempQuat[0] = 1.0;
    track->tempQuat[1] = track->tempQuat[2] = track->tempQuat[3] = 0.0;
    track->trans[0] += track->tempTrans[0];
    track->trans[1] += track->tempTrans[1];
    track->trans[2] += track->tempTrans[2];
/*    f = track->trans[0] * track->trans[0] + track->trans[1] * track->trans[1] + track->trans[2] * track->trans[2];
    if (f > 1.0) {
        f = 1.0 / sqrt(f);
        track->trans[0] *= f;
        track->trans[1] *= f;
        track->trans[2] *= f;
    } */
    track->tempTrans[0] = track->tempTrans[1] = track->tempTrans[2] = 0.0;
    track->scale += track->tempScale;
    if (track->scale > 5.0)
        track->scale = 5.0;
    else if (track->scale < -5.0)
        track->scale = -5.0;
    track->tempScale = 0.0;
}

