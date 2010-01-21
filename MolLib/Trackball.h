/*
 *  Trackball.h
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

#ifndef __Trackball_H__
#define __Trackball_H__

#include <math.h>

#ifndef kRad2Deg
#define kRad2Deg  (180./3.14159265358979)
#endif
#ifndef kDeg2Rad
#define kDeg2Rad  (3.14159265358979 / 180.)
#endif
#ifndef kLog10
#define kLog10 2.3025851
#endif
#ifndef kCot15Deg
#define kCot15Deg 3.73205075
#endif

#ifdef __cplusplus
extern "C" {
#endif
	
typedef enum TrackballMode {
	kTrackballRotateMode = 1,
	kTrackballTranslateMode,
	kTrackballScaleMode,
	/*  The following "modes" are not related to the trackball operation,
	    but maybe useful for use in molecular modeling applications  */
	kTrackballSelectionMode,
	kTrackballCreateMode,
	kTrackballEraseMode
} TrackballMode;

/*  Simulate the trackball operation. The screen coordinates should be converted
 *  so that the unit circle represents the trackball (i.e. the center of the
 *  view is (0, 0), and the half of the screen height (or width) is 1.0)  */
typedef struct Trackball {
	int refCount;
    float start[2];     /*  Drag started here  */
    float startQuat[4];	/*  Rotation from the start point to the center (of the trackball) */
    TrackballMode   mode;

    float quat[4];		/*  Rotation (in quaternion)  */
    float tempQuat[4];	/*  Temporary rotation during dragging of trackball  */
    float trans[3];		/*  Translation  */
    float tempTrans[3];	/*  Temporary translation  */
    float scale;		/*  Scale  */
    float tempScale;	/*  Temporary scale */
} Trackball;

Trackball *TrackballNew(void);
void TrackballRetain(Trackball *track);
void TrackballRelease(Trackball *track);

void TrackballGetRotate(const Trackball *track, float *a);
void TrackballGetTranslate(const Trackball *track, float *a);
void TrackballGetPerspective(const Trackball *track, float *a);

void TrackballReset(Trackball *track);
void TrackballSetScale(Trackball *track, float scale);
void TrackballSetTranslate(Trackball *track, const float *a);

void TrackballStartDragging(Trackball *track, const float *mousePos, TrackballMode mode);
void TrackballSetTemporaryRotation(Trackball *track, const float *q);
void TrackballDrag(Trackball *track, const float *mousePos);
void TrackballEndDragging(Trackball *track, const float *mousePos);

#ifdef __cplusplus
}
#endif
		
#endif /* __Trackball_H__  */
