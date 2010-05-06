/*
 *  MainView.c
 *
 *  Created by Toshi Nagata on 06/07/30.
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

#include "MolLib.h"

#include "MD/MDCore.h"
#include "MD/MDGraphite.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define biso2radius(r) ((r) > 0.5 ? sqrt((r) / 78.9568352087147) : 0.08)

/*  Invalid bond/angle/torsion value, used in internal cache  */
const Double kInvalidValue = -10000000.0;

#pragma mark ==== MainView public methods ====

MainView *
MainView_newMainView(void *ref)
{
	static GLdouble sIdentity[16] = {
		1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
	};
	MainView *mview = (MainView *)malloc(sizeof(MainView));
	if (mview != NULL) {
		memset(mview, 0, sizeof(MainView));
		mview->ref = ref;
		mview->track = TrackballNew();
		mview->mode = kTrackballRotateMode;
		mview->tempAtoms[0] = mview->tempAtoms[1] = -1;
		memmove(mview->modelview_matrix, sIdentity, sizeof(sIdentity));
		memmove(mview->projection_matrix, sIdentity, sizeof(sIdentity));
		mview->atomRadius = 0.4;
		mview->bondRadius = 0.1;
		mview->probabilityScale = 1.5382;
		mview->dimension = 10.0;
	/*	MainView_resizeToFit(mview); */
		mview->showHydrogens = mview->showDummyAtoms = mview->showExpandedAtoms = 1;
		mview->showPeriodicBox = 1;
		mview->showGraphite = 5;
		mview->tableCache = IntGroupNew();
		mview->tableSelection = IntGroupNew();
	}
	return mview;
}

void
MainView_release(MainView *mview)
{
	int i;
	if (mview != NULL) {
		MainView_setMolecule(mview, NULL);
		TrackballRelease(mview->track);
		IntGroupRelease(mview->tableCache);
		IntGroupRelease(mview->tableSelection);
		if (mview->nlabels > 0) {
			for (i = 0; i < mview->nlabels; i++) {
				MainViewCallback_releaseLabel(mview->labels[i].label);
			}
			free(mview->labels);
			free(mview->sortedLabels);
		}
		if (mview->rotateFragment != NULL)
			IntGroupRelease(mview->rotateFragment);
		if (mview->rotateFragmentOldPos != NULL)
			free(mview->rotateFragmentOldPos);
		if (mview->visibleFlags != NULL)
			free(mview->visibleFlags);
	}
	free(mview);
}

void
MainView_setMolecule(MainView *mview, struct Molecule *mol)
{
	if (mview == NULL || mview->mol == mol)
		return;
	if (mview->mol != NULL) {
		mview->mol->mview = NULL;  /*  No need to release  */
		MoleculeRelease(mview->mol);
	}
	mview->mol = mol;
	if (mol != NULL) {
		MoleculeRetain(mol);
		mol->mview = mview;  /*  No retain  */
		MainViewCallback_moleculeReplaced(mview, mol);
		MainView_resizeToFit(mview);
		MoleculeCallback_notifyModification(mol, 0);
	/*	MainViewCallback_setNeedsDisplay(mview, 1); */
	}
}

int
MainView_isAtomHidden(MainView *mview, int index)
{
	if (index < 0 || index >= mview->mol->natoms)
		return 1;
	else if (mview->visibleFlags == NULL)
		return 0;
	else return (mview->visibleFlags[index] == 0);
/*	Atom *ap = ATOM_AT_INDEX(mview->mol->atoms, index);
	if (ap->exflags & kAtomHiddenFlag)
		return 1;
	if (!mview->showHydrogens && ap->atomicNumber == 1)
		return 1;
	if (!mview->showDummyAtoms && ap->atomicNumber == 0)
		return 1;
	return 0; */
}

void
MainView_refreshCachedInfo(MainView *mview)
{
	Molecule *mol;
	int i, j, n1, n2, n3, n4;
	Byte f1, f2, f3, f4;
	Atom *ap;

	if (mview == NULL || (mol = mview->mol) == NULL)
		return;

	/*  Rebuild internal caches  */
	
	/*  Visible flags  */
	if (mview->visibleFlags == NULL)
		mview->visibleFlags = (Byte *)malloc(mol->natoms);
	else
		mview->visibleFlags = (Byte *)realloc(mview->visibleFlags, mol->natoms);
	memset(mview->visibleFlags, 0, mol->natoms);
	mview->countHidden = 0;
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->exflags & kAtomHiddenFlag) {
			mview->countHidden++;
			continue;
		}
		if (!mview->showHydrogens && ap->atomicNumber == 1)
			continue;
		if (!mview->showDummyAtoms && ap->atomicNumber == 0)
			continue;
		mview->visibleFlags[i] = 1;
	}
	
	/*  Selection (only temporary used within this function)  */
	{
		IntGroup *sel = MoleculeGetSelection(mol);
		if (sel != NULL && IntGroupGetCount(sel) > 0) {
			for (i = 0; (n1 = IntGroupGetStartPoint(sel, i)) >= 0; i++) {
				n2 = IntGroupGetEndPoint(sel, i);
				if (n2 > mol->natoms)
					n2 = mol->natoms;
				while (n1 < n2)
					mview->visibleFlags[n1++] |= 2;
			}
		}
	}
	
	/*  Table caches  */
	if (mview->tableCache == NULL)
		mview->tableCache = IntGroupNew();
	if (mview->tableSelection == NULL)
		mview->tableSelection = IntGroupNew();
	IntGroupClear(mview->tableCache);
	IntGroupClear(mview->tableSelection);
		
	if (mview->tableIndex == kMainViewAtomTableIndex) {  /* Atoms */
		for (i = j = 0; i < mol->natoms; i++) {
			f1 = mview->visibleFlags[i];
			if ((f1 & 1) != 0) {
				IntGroupAdd(mview->tableCache, i, 1);
				if ((f1 & 2) != 0)
					IntGroupAdd(mview->tableSelection, j, 1);
				j++;
			}
		}
	} else if (mview->tableIndex == kMainViewBondTableIndex) {  /* Bonds */
		for (i = j = 0; i < mol->nbonds; i++) {
			n1 = mol->bonds[i * 2];
			n2 = mol->bonds[i * 2 + 1];
			f1 = mview->visibleFlags[n1];
			f2 = mview->visibleFlags[n2];
			if ((f1 & 1) != 0 && (f2 & 1) != 0) {
				IntGroupAdd(mview->tableCache, i, 1);
				if ((f1 & 2) != 0 && (f2 & 2) != 0)
					IntGroupAdd(mview->tableSelection, j, 1);
				j++;
			}
		}
	} else if (mview->tableIndex == kMainViewAngleTableIndex) {  /* Angles */
		for (i = j = 0; i < mol->nangles; i++) {
			n1 = mol->angles[i * 3];
			n2 = mol->angles[i * 3 + 1];
			n3 = mol->angles[i * 3 + 2];
			f1 = mview->visibleFlags[n1];
			f2 = mview->visibleFlags[n2];
			f3 = mview->visibleFlags[n3];
			if ((f1 & 1) != 0 && (f2 & 1) != 0 && (f3 & 1) != 0) {
				IntGroupAdd(mview->tableCache, i, 1);
				if ((f1 & 2) != 0 && (f2 & 2) != 0 && (f3 & 2) != 0)
					IntGroupAdd(mview->tableSelection, j, 1);
				j++;
			}
		}
	} else if (mview->tableIndex == kMainViewDihedralTableIndex) {  /* Dihedrals */
		for (i = j = 0; i < mol->ndihedrals; i++) {
			n1 = mol->dihedrals[i * 4];
			n2 = mol->dihedrals[i * 4 + 1];
			n3 = mol->dihedrals[i * 4 + 2];
			n4 = mol->dihedrals[i * 4 + 3];
			f1 = mview->visibleFlags[n1];
			f2 = mview->visibleFlags[n2];
			f3 = mview->visibleFlags[n3];
			f4 = mview->visibleFlags[n4];
			if ((f1 & 1) != 0 && (f2 & 1) != 0 && (f3 & 1) != 0 && (f4 & 1) != 0) {
				IntGroupAdd(mview->tableCache, i, 1);
				if ((f1 & 2) != 0 && (f2 & 2) != 0 && (f3 & 2) != 0 && (f4 & 2) != 0)
					IntGroupAdd(mview->tableSelection, j, 1);
				j++;
			}
		}
	} else if (mview->tableIndex == kMainViewImproperTableIndex) {  /* Impropers */
		for (i = j = 0; i < mol->nimpropers; i++) {
			n1 = mol->impropers[i * 4];
			n2 = mol->impropers[i * 4 + 1];
			n3 = mol->impropers[i * 4 + 2];
			n4 = mol->impropers[i * 4 + 3];
			f1 = mview->visibleFlags[n1];
			f2 = mview->visibleFlags[n2];
			f3 = mview->visibleFlags[n3];
			f4 = mview->visibleFlags[n4];
			if ((f1 & 1) != 0 && (f2 & 1) != 0 && (f3 & 1) != 0 && (f4 & 1) != 0) {
				IntGroupAdd(mview->tableCache, i, 1);
				if ((f1 & 2) != 0 && (f2 & 2) != 0 && (f3 & 2) != 0 && (f4 & 2) != 0)
					IntGroupAdd(mview->tableSelection, j, 1);
				j++;
			}
		}
	} else if (mview->tableIndex == kMainViewParameterTableIndex) {  /* Parameter infos */
		/*  Do nothing (tableCache will not be used)  */
	} else if (mview->tableIndex == kMainViewMOTableIndex) {  /* MO infos  */
		/*  Really no need to cache info, but create it anyway to simplify code  */
		if (mol->bset != NULL && mol->bset->nmos > 0)
			IntGroupAdd(mview->tableCache, 0, mol->bset->nmos);
	}
	
	/*  Clear internal selection flag  */
	for (i = 0; i < mol->natoms; i++) {
		mview->visibleFlags[i] &= ~2;
	}
}

#pragma mark ====== 2D/3D transform operations ======

void
MainView_resizeToFit(MainView *mview)
{
	Vector p;
	float f[4];
	if (mview == NULL || mview->mol == NULL)
		return;
	if (mview->mol->natoms == 0) {
		TrackballReset(mview->track);
		return;
	}
	MoleculeCenterOfMass(mview->mol, &p, NULL);
/*	if (mview->mol->is_xtal_coord)
		TransformVec(&p, mview->mol->cell->tr, &p); */
	f[0] = -p.x / mview->dimension;
	f[1] = -p.y / mview->dimension;
	f[2] = -p.z / mview->dimension;
	TrackballSetTranslate(mview->track, f);

	/*  Set scale
		r0: the longest distance from the center of mass
		r0 > dimension: scale = -log(r0/dimension)/log(10)  (negative)
		r0 < dimension: scale = -log(atan2(r0, dimension*cot(15deg))*180deg/pi*2/30deg)/log(10) (positive)
	*/
	{
		int i;
		Vector q;
		Atom *ap;
		double r0 = 0.0, r1, scale;
		for (i = 0, ap = mview->mol->atoms; i < mview->mol->natoms; i++, ap = ATOM_NEXT(ap)) {
			q = ap->r;
		/*	if (mview->mol->is_xtal_coord)
				TransformVec(&q, mview->mol->cell->tr, &q); */
			VecDec(q, p);
			r1 = VecLength(q);
			if (r1 > r0)
				r0 = r1;
		}
		r0 /= mview->dimension;
		if (r0 < 1e-6)
			scale = 0.0;
		else if (r0 < 1.0)
			scale = -log(atan2(r0, kCot15Deg) * kRad2Deg * 2 / 30.0) / kLog10;
		else
			scale = -log(r0) / kLog10;
		TrackballSetScale(mview->track, scale);
	}

	MainViewCallback_setNeedsDisplay(mview, 1);

#if 0
	GLdouble mm[4][4];
	GLdouble p0[4];
	int i, j, natoms;
	const Atom *ap;
	const ElementPar *app;
	GLdouble *pp, *p1;
	GLdouble max[3], min[3], center[3];
	float frame[4], width, height, cot15, rsin15, cot_th, rsin_th, d, trans[3];
	Molecule *mol;
	
	if (mview == NULL || mview->mol == NULL)
		return;
	mol = mview->mol;
	
	/*  Transform the coordinates of all atoms with the rotation part of
	 *  the model-view matrix  */
	memmove((GLdouble *)mm, mview->modelview_matrix, sizeof(float) * 16);
	natoms = mol->natoms;
	pp = (GLdouble *)calloc(sizeof(GLdouble), 4 * natoms);
	MALLOC_CHECK(pp, "resizing the model to fit");
	for (i = 0; i < natoms; i++) {
		ap = ATOM_AT_INDEX(mol->atoms, i);
		if (ap == NULL)
			continue;
		p0[0] = ap->r.x; p0[1] = ap->r.y; p0[2] = ap->r.z;
		p0[3] = 0.0;  /*  This disables the translation part  */
		p1 = &pp[i * 4];
		MAT_DOT_VEC_4X4(p1, mm, p0);  /* #define'd in GLUT/vvector.h */
	}
	
	/*  Determine the center for each axis  */
	max[0] = max[1] = max[2] = -1e20;
	min[0] = min[1] = min[2] = 1e20;
	for (i = 0; i < natoms; i++) {
		p1 = &pp[i * 4];
		for (j = 0; j < 3; j++) {
			if (p1[j] > max[j])
				max[j] = p1[j];
			if (p1[j] < min[j])
				min[j] = p1[j];
		}
	}
	for (j = 0; j < 3; j++)
		center[j] = (min[j] + max[j]) * 0.5;
	
	/*  Get the frame size  */
	MainViewCallback_frame(mview, frame);
	width = frame[2] - frame[0];
	height = frame[3] - frame[1];

	/*  Calculate the minimum distance from which the view angle of
	 *  the atom becomes no more than 15 degree (for the y direction)
	 *  and theta degree (for the x direction, theta = arctan(h/w*tan(15 deg))) */
	cot15 = 3.73205080756888;
	rsin15 = 3.86370330515628;
	cot_th = width / height * cot15;
	rsin_th = sqrt(1.0 + cot_th * cot_th);
	d = 0;
	for (i = 0; i < natoms; i++) {
		float dx, dy, z, r;
		ap = ATOM_AT_INDEX(mol->atoms, i);
		if (ap == NULL)
			continue;
		app = &(gBuiltinParameters->atomPars[ap->atomicNumber]);
		if (app == NULL)
			continue;
		r = app->radius * mview->atomRadius;
		z = pp[i * 4 + 2] - center[2];
		dx = pp[i * 4] * cot_th + r * rsin_th + z;
		dy = pp[i * 4 + 1] * cot15 + r * rsin15 + z;
		if (d < dx)
			d = dx;
		if (d < dy)
			d = dy;
	}
	mview->dimension = d / cot15;
	TrackballSetScale(mview->track, 1.0);
	trans[0] = center[0] / mview->dimension;
	trans[1] = center[1] / mview->dimension;
	trans[2] = center[2] / mview->dimension;
	TrackballSetTranslate(mview->track, trans);

	free(pp);
#endif
}

int
MainView_convertScreenPositionToObjectPosition(MainView *mview, const GLfloat *screenPos, GLfloat *objectPos)
{
	float rect[4];
    GLint viewport[4], n;
    GLfloat winZ;
    GLdouble posX, posY, posZ;
	if (mview == NULL)
		return 0;
	MainViewCallback_frame(mview, rect);
    viewport[0] = viewport[1] = 0;
    viewport[2] = (GLint)(rect[2] - rect[0]);
    viewport[3] = (GLint)(rect[3] - rect[1]);
	MainViewCallback_lockFocus(mview);
    if (screenPos[2] >= 0.0)
        winZ = screenPos[2];
    else
        glReadPixels(screenPos[0], screenPos[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
    gluUnProject(screenPos[0], screenPos[1], winZ, mview->modelview_matrix, mview->projection_matrix, viewport, &posX, &posY, &posZ);
    n = glGetError();
	MainViewCallback_unlockFocus(mview);
	objectPos[0] = posX;
	objectPos[1] = posY;
	objectPos[2] = posZ;
    if (n != GL_NO_ERROR || winZ == 1.0)
        return 0;
    else
        return 1;
}

int
MainView_convertObjectPositionToScreenPosition(MainView *mview, const GLfloat *objectPos, GLfloat *screenPos)
{
	float rect[4];
    GLint viewport[4];
    GLdouble objX, objY, objZ;
	if (mview == NULL)
		return 0;
	MainViewCallback_frame(mview, rect);
    viewport[0] = viewport[1] = 0;
    viewport[2] = (GLint)(rect[2] - rect[0]);
    viewport[3] = (GLint)(rect[3] - rect[1]);
    gluProject(objectPos[0], objectPos[1], objectPos[2], mview->modelview_matrix, mview->projection_matrix, viewport, &objX, &objY, &objZ);
    if (glGetError() == GL_NO_ERROR) {
		screenPos[0] = objX;
		screenPos[1] = objY;
		screenPos[2] = objZ;
	/*	fprintf(stderr, "object(%.3f,%.3f,%.3f) screen(%.3f,%.3f,%.3f)\n", objectPos[0], objectPos[1], objectPos[2], screenPos[0], screenPos[1], screenPos[2]); */
		return 1;
	} else return 0;
}

int
MainView_findObjectAtPoint(MainView *mview, const float *mousePos, int *outIndex1, int *outIndex2, int mouseCheck, int ignoreExpandedAtoms)
{
	float screenPos[3], op[3], oq[3], pqlen, pqlen2;
	Vector pq, pa, v1, r1, r2;
    int i, natoms, nbonds;
	float r, d2, z;
	const Atom *ap, *bp;
	const ElementPar *ep;
	const int *ip;
	Molecule *mol;
	float minDepth;
	int n1, n2;

	if (mview == NULL || mview->mol == NULL) {
		*outIndex1 = *outIndex2 = -1;
		return 0;
	}
	mol = mview->mol;

#if 1
	screenPos[0] = mousePos[0];
	screenPos[1] = mousePos[1];
	screenPos[2] = -1.0;
	if (MainView_convertScreenPositionToObjectPosition(mview, screenPos, op) == 0)
		return 0;  /*  Nothing is here  */

	/*  PQ is the part of the eyesight line in the visible area  */
	screenPos[2] = 0.0;
	MainView_convertScreenPositionToObjectPosition(mview, screenPos, op);
	screenPos[2] = 1.0;
	MainView_convertScreenPositionToObjectPosition(mview, screenPos, oq);
	pq.x = oq[0] - op[0];
	pq.y = oq[1] - op[1];
	pq.z = oq[2] - op[2];
	pqlen2 = VecLength2(pq);
	pqlen = sqrt(pqlen2);
	natoms = mol->natoms;
	n1 = n2 = -1;
	minDepth = 100.0;
	for (i = 0; i < natoms; i++) {
		Vector pq1, pa1;
		float pq1len2, pq1len;
		if (mouseCheck && i % 50 == 0 && MainViewCallback_mouseCheck(mview))
			return 0;  /*  If mouse event is detected return immediately  */
		/*  Examine if an atom is visible or not  */
		/*  The distance of the atom center (A) from line PQ: */
		/*    d = |VecCross(PA,PQ)|/|PQ|  */
		/*  z = VecDot(PA,PQ)/|PQ|^2 - sqrt(r^2 - d^2)/|PQ|  */
        ap = ATOM_AT_INDEX(mol->atoms, i);
		if (ap == NULL)
			continue;
		if (MainView_isAtomHidden(mview, i))
			continue;
		if (ignoreExpandedAtoms && SYMOP_ALIVE(ap->symop))
			continue;
		r1 = ap->r;
	/*	if (mol->is_xtal_coord)
			TransformVec(&r1, mol->cell->tr, &r1); */
		pa.x = r1.x - op[0];
		pa.y = r1.y - op[1];
		pa.z = r1.z - op[2];
		if (ap->aniso != NULL) {
			/*  Convert to ellipsoid principal axes  */
			Mat33 m1;
			Aniso *anp = ap->aniso;
			MatrixInvert(m1, anp->pmat);
			MatrixVec(&pq1, m1, &pq);
			MatrixVec(&pa1, m1, &pa);
			r = mview->probabilityScale;
			pq1len2 = VecLength2(pq1);
			pq1len = sqrt(pq1len2);
		} else {
			if (mview->showEllipsoids) {
				r = biso2radius(ap->tempFactor) * mview->probabilityScale;
			} else {
				ep = &(gElementParameters[ap->atomicNumber]);
				if (ep == NULL)
					continue;
				r = ep->radius * mview->atomRadius;
			}
			pa1 = pa;
			pq1 = pq;
			pq1len2 = pqlen2;
			pq1len = pqlen;
		}
		VecCross(v1, pa1, pq1);
		d2 = VecLength2(v1) / pq1len2;
		if (d2 > r * r)
			continue;  /*  Not visible  */
		z = VecDot(pa1, pq1) / pq1len2 - sqrt(r * r - d2) / pq1len;
		if (z < 0.0 || z > 1.0)
			continue;  /*  Out of viewing volume  */
		if (z < minDepth) {
			minDepth = z;
			n1 = i;
		}
	}
	nbonds = mol->nbonds;
	for (i = 0; i < nbonds; i++) {
		Vector vx, vy, vz, vv, vp;
		Double wb, wa, t, wx, blen;
		if (mouseCheck && i % 50 == 0 && MainViewCallback_mouseCheck(mview))
			return 0;  /*  If mouse event is detected return immediately  */
		/*  Examine if a bond is visible or not  */
		ip = &(mol->bonds[i * 2]);
		ap = ATOM_AT_INDEX(mol->atoms, ip[0]);
		bp = ATOM_AT_INDEX(mol->atoms, ip[1]);
		if (MainView_isAtomHidden(mview, ip[0]) || MainView_isAtomHidden(mview, ip[1]))
			continue;
		if (ignoreExpandedAtoms && SYMOP_ALIVE(ap->symop) && SYMOP_ALIVE(bp->symop))
			continue;
		/*  vx/vy/vz is a orthonormal base in which AB parallels the x-axis
		    and AP in in the xy plane  */
		/*  vp and vv is AP and PQ in that coordinate system  */
		r1 = ap->r;
		r2 = bp->r;
	/*	if (mol->is_xtal_coord) {
			TransformVec(&r1, mol->cell->tr, &r1);
			TransformVec(&r2, mol->cell->tr, &r2);
		} */
		v1.x = op[0] - r1.x;
		v1.y = op[1] - r1.y;
		v1.z = op[2] - r1.z;
		VecSub(vx, r2, r1);
		blen = sqrt(VecLength2(vx));
		if (blen < 1e-10)
			continue;
		vx.x /= blen;
		vx.y /= blen;
		vx.z /= blen;
		VecCross(vz, vx, v1);
		if (NormalizeVec(&vz, &vz))
			continue;
		VecCross(vy, vz, vx);
		vp.x = VecDot(v1, vx);
		vp.y = VecDot(v1, vy);
		vp.z = VecDot(v1, vz);
		vv.x = VecDot(pq, vx);
		vv.y = VecDot(pq, vy);
		vv.z = VecDot(pq, vz);
		/*  The bond surface is y^2 + z^2 = r^2, 0 <= x <= 1  */
		/*  The eyesight line is (x,y,z) = vp + t * vv, 0 <= t <= 1  */
		/*  The crossing point: t = (-vv.y*vp.y - sqrt((vv.y*vp.y)^2 - (vv.y^2+vv.z^2)(vp.y^2-r^2)))/(vv.y^2+vv.z^2)  */
		/*  (Note that vp.z = 0 by definition)  */
		r = mview->bondRadius;
		wb = vv.y * vp.y;
		wa = vv.y * vv.y + vv.z * vv.z;
		d2 = wb * wb - wa * (vp.y * vp.y - r * r);
		if (d2 < 0)
			continue;  /*  Not visible  */
		t = (-wb - sqrt(d2)) / wa;
		if (t < 0 || t > 1)
			continue;  /*  Out of visible volume  */
		wx = vp.x + t * vv.x;
		if (wx < 0 || wx > blen)
			continue;  /*  Outside of bond  */
		if (t < minDepth) {
			minDepth = t;
			n1 = ip[0];
			n2 = ip[1];
		}
	}
	*outIndex1 = n1;
	*outIndex2 = n2;
	return (n1 >= 0 || n2 >= 0);
		
#else
	/*  Examine if anything is drawn at the point  */
	screenPos[0] = mousePos[0];
	screenPos[1] = mousePos[1];
	screenPos[2] = -1.0;
	if (MainView_convertScreenPositionToObjectPosition(mview, screenPos, objectPos) == 0)
		return 0;  /*  Nothing is here  */
	op[0] = objectPos[0];
	op[1] = objectPos[1];
	op[2] = objectPos[2];

	/*  Examine the distance from the atom center with op, and select the one
	    with the smallest difference with the radius on the screen  */
	atomRadius = mview->atomRadius;
	bondRadius = mview->bondRadius;
	natoms = mol->natoms;
	mindr = 100.0;
	n1 = n2 = -1;
    for (i = 0; i < natoms; i++) {
        ap = ATOM_AT_INDEX(mol->atoms, i);
		if (ap == NULL)
			continue;
		dp = &(gBuiltinParameters->atomPars[ap->atomicNumber]);
		if (dp == NULL)
			continue;
		p[0] = objectPos[0] - ap->r.x;
		p[1] = objectPos[1] - ap->r.y;
		p[2] = objectPos[2] - ap->r.z;
		rr = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) / (dp->radius * atomRadius);
		dr = fabs(1.0 - rr);
		if (dr < mindr) {
			mindr = dr;
			n1 = i;
		}
	}
	if (mindr > 0.1)
		n1 = -1;
	nbonds = mol->nbonds;
    for (i = 0; i < nbonds; i++) {
		float q;
		ip = &(mol->bonds[i * 2]);
		ap = ATOM_AT_INDEX(mol->atoms, ip[0]);
		bp = ATOM_AT_INDEX(mol->atoms, ip[1]);
		p[0] = objectPos[0] - ap->r.x;
		p[1] = objectPos[1] - ap->r.y;
		p[2] = objectPos[2] - ap->r.z;
		b[0] = bp->r.x - ap->r.x;
		b[1] = bp->r.y - ap->r.y;
		b[2] = bp->r.z - ap->r.z;
		/*  ap.r + q*(b[0],b[1],b[2]) is the closest point on line AB to objectPos */
		/*  If q is outside [0, 1] then the point is outside the line segment AB
		    so it is ignored  */
		q = (p[0]*b[0]+p[1]*b[1]+p[2]*b[2])/(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
		if (q < 0 || q > 1)
			continue;
		p[0] -= q * b[0];
		p[1] -= q * b[1];
		p[2] -= q * b[2];
		dr = fabs(bondRadius - sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]));
		if (dr < mindr) {
			mindr = dr;
			n1 = ip[0];
			n2 = ip[1];
		}
    }
	*outIndex1 = n1;
	*outIndex2 = n2;
	return (n1 >= 0 || n2 >= 0);
#endif
}

int
MainView_screenCenterPointOfAtom(MainView *mview, int index, float *outScreenPos)
{
	const Atom *ap;
	const ElementPar *dp;
	Vector cv, pv, v;
	Double rad, w;
	float p[3];

	if (mview == NULL || mview->mol == NULL || index < 0 || index >= mview->mol->natoms)
		return 0;
	ap = ATOM_AT_INDEX(mview->mol->atoms, index);
	
	/*  Camera position in object coordinates  */
	cv = mview->camera;
	
	/*  The atom position (in Cartesian)  */
	v = ap->r;
/*	if (mview->mol->is_xtal_coord)
		TransformVec(&v, mview->mol->cell->tr, &v); */

	/*  The vector from atom center to camera  */
	VecSub(pv, cv, v);
	
	/*  Get the surface point of the ellipsoid/sphere along the camera vector  */
	if (mview->showEllipsoids) {
		Mat33 m1;
		Aniso *anp = ap->aniso;
		if (anp != NULL) {
			Vector vx, vy, vz;
			MatrixInvert(m1, anp->pmat);
			/*  Convert the 'screen plane' vectors to ellipsoid principal axes  */
			vy = mview->up;
			VecCross(vx, mview->lookto, vy);
			MatrixVec(&vx, m1, &vx);
			MatrixVec(&vy, m1, &vy);
			/*  Touching point of the 'screen plane' to the ellipsoid  */
			VecCross(vz, vx, vy);
			w = mview->probabilityScale / VecLength(vz) * 1.1;
			VecScaleSelf(vz, w);
			MatrixVec(&vz, anp->pmat, &vz);
			/*  The crossing point of the camera vector with the touching plane */
			w = fabs(VecDot(pv, vz) / VecLength2(pv));
			VecScaleSelf(pv, w);
		} else {
			w = mview->probabilityScale * biso2radius(ap->tempFactor) / VecLength(pv) * 1.1;
			VecScaleSelf(pv, w);
		}
	} else {
		dp = &(gElementParameters[ap->atomicNumber]);
		rad = dp->radius * mview->atomRadius;
		w = rad / VecLength(pv) * 1.1;
		VecScaleSelf(pv, w);
	}
	VecInc(v, pv);
	
	/*  Transform to screen coordinate  */
	p[0] = v.x;
	p[1] = v.y;
	p[2] = v.z;
	return MainView_convertObjectPositionToScreenPosition(mview, p, outScreenPos);
}

void
MainView_getCamera(MainView *mview, Vector *outCamera, Vector *outLookAt, Vector *outUp)
{
	if (mview != NULL) {
		*outCamera = mview->camera;
		*outLookAt = mview->lookat;
		*outUp = mview->up;
	}
}

#pragma mark ====== Draw model ======

void
MainView_initializeOpenGLView(MainView *mview)
{
	static GLfloat ambient[] = {0.6, 0.6, 0.6, 1.0};  // Some white ambient light.
	static GLfloat diffuse[] = {1.0, 1.0, 1.0, 1.0};  // A white light.
	
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	
	//  Set the ambient light
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
	
	//  Set the light and switch it on
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glEnable(GL_LIGHT0);
	
	//  Enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if (mview != NULL)
		mview->isInitialized = 1;
}

/*  Get orthogonal unit vectors  */
static int
getOrthogonalVectors(const GLfloat *ap, GLfloat *bp, GLfloat *cp)
{
    Double ra, rb;
    ra = sqrt(ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2]);
    if (ra < 1e-20)
        return 0;
    ra = 1 / ra;
    if (fabs(ap[0]) < fabs(ap[1])) {
        if (fabs(ap[0]) < fabs(ap[2])) {
            bp[0] = 0;
            bp[1] = -ap[2];
            bp[2] = ap[1];
        } else {
            bp[0] = ap[1];
            bp[1] = -ap[0];
            bp[2] = 0;
        }
    } else {
        if (fabs(ap[1]) < fabs(ap[2])) {
            bp[0] = -ap[2];
            bp[1] = 0;
            bp[2] = ap[0];
        } else {
            bp[0] = ap[1];
            bp[1] = -ap[0];
            bp[2] = 0;
        }
    }
    rb = 1 / sqrt(bp[0] * bp[0] + bp[1] * bp[1] + bp[2] * bp[2]);
    bp[0] *= rb;
    bp[1] *= rb;
    bp[2] *= rb;
    cp[0] = ra * (ap[1] * bp[2] - ap[2] * bp[1]);
    cp[1] = ra * (ap[2] * bp[0] - ap[0] * bp[2]);
    cp[2] = ra * (ap[0] * bp[1] - ap[1] * bp[0]);
/*    printf("a = (%f, %f, %f) b = (%f, %f, %f) c = (%f, %f, %f)\n",
        ap[0], ap[1], ap[2], bp[0], bp[1], bp[2], cp[0], cp[1], cp[2]); */
    return 1;
}

static GLfloat sSinCache[81];
static int sSinCacheSect = 0;

static int
setSinCache(int sect)
{
    int n, m, i;
    m = sect / 4;
    n = m * 4;
    if (n >= 64)
        n = 64;
    if (n != sSinCacheSect) {
        sSinCacheSect = n;
        for (i = 0; i <= m * 5; i++)
            sSinCache[i] = sin(3.14159265358979 * 2 / n * i);
    }
    return n;
}

static void
drawCylinder(const GLfloat *a, const GLfloat *b, GLfloat r, int sect)
{
    GLfloat *c, *s;
    int n, i;
    GLfloat d[3], v[3], w[3];
    n = setSinCache(sect);
    if (n <= 0)
        return;
    s = sSinCache;
    c = &sSinCache[n/4];
    d[0] = b[0] - a[0];
    d[1] = b[1] - a[1];
    d[2] = b[2] - a[2];
    if (getOrthogonalVectors(d, v, w) == 0)
        return;
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i <= n; i++) {
        float nx, ny, nz;
        nx = v[0] * c[i] + w[0] * s[i];
        ny = v[1] * c[i] + w[1] * s[i];
        nz = v[2] * c[i] + w[2] * s[i];
        glNormal3f(nx, ny, nz);
        glVertex3f(a[0] + r * nx, a[1] + r * ny, a[2] + r * nz);
        glVertex3f(b[0] + r * nx, b[1] + r * ny, b[2] + r * nz);
    }
    glEnd();    
}

static void
drawSphere(const GLfloat *p, GLfloat r, int sect)
{
    GLfloat *c, *s;
    int n, i, j;
    n = setSinCache(sect);
    if (n <= 0)
        return;
    s = sSinCache;
    c = &sSinCache[n/4];
    for (i = 0; i <= n; i++) {
        glBegin(GL_QUAD_STRIP);
        for (j = 1; j <= n / 2 - 1; j++) {
            glNormal3f(s[j] * c[i], s[j] * s[i], c[j]);
            glVertex3f(r * s[j] * c[i] + p[0], r * s[j] * s[i] + p[1], r * c[j] + p[2]);
            glNormal3f(s[j] * c[i+1], s[j] * s[i+1], c[j]);
            glVertex3f(r * s[j] * c[i+1] + p[0], r * s[j] * s[i+1] + p[1], r * c[j] + p[2]);
        }
        glEnd();
    }
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex3f(p[0], p[1], r + p[2]);
    for (i = n; i >= 0; i--) {
        glNormal3f(s[1] * c[i], s[1] * s[i], c[1]);
        glVertex3f(r * s[1] * c[i] + p[0], r * s[1] * s[i] + p[1], r * c[1] + p[2]);
    }
    glEnd();
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, -1);
    glVertex3f(p[0], p[1], -r + p[2]);
    for (i = 0; i <= n; i++) {
        glNormal3f(s[1] * c[i], s[1] * s[i], -c[1]);
        glVertex3f(r * s[1] * c[i] + p[0], r * s[1] * s[i] + p[1], -r * c[1] + p[2]);
    }
    glEnd();
}

static void
drawEllipsoid(const GLfloat *p, const GLfloat *v1, const GLfloat *v2, const GLfloat *v3, int sect)
{
	GLfloat mat[16];
	static const GLfloat origin[3] = {0, 0, 0};
	mat[0] = v1[0]; mat[1] = v1[1]; mat[2] = v1[2]; mat[3] = 0;
	mat[4] = v2[0]; mat[5] = v2[1]; mat[6] = v2[2]; mat[7] = 0;
	mat[8] = v3[0]; mat[9] = v3[1]; mat[10] = v3[2]; mat[11] = 0;
	mat[12] = p[0]; mat[13] = p[1]; mat[14] = p[2]; mat[15] = 1;
    glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glMultMatrixf(mat);
	glEnable(GL_NORMALIZE);
	//	glutSolidSphere(1, sect, sect); /* Is this faster than my code? */
	drawSphere(origin, 1, sect);
	glDisable(GL_NORMALIZE);
	glPopMatrix();
}

static char *
temporarySelection(MainView *mview, int flags, int clickCount, int ignoreExpandedAtoms)
{
	char *selectFlags;
	int i, natoms;
	const Atom *ap;
	Double rect[4];
	natoms = mview->mol->natoms;
	selectFlags = (char *)calloc(sizeof(char), natoms);
	if (selectFlags == NULL)
		return NULL;
	if (clickCount > 0) {
		int n1, n2;
		if (MainView_findObjectAtPoint(mview, mview->dragStartPos, &n1, &n2, 0, 0)) {
			if (n1 >= 0 && n1 < natoms)
				selectFlags[n1] = 1;
			if (n2 >= 0 && n2 < natoms)
				selectFlags[n2] = 1;
		}
	} else {
		if (mview->dragStartPos[0] < mview->dragEndPos[0]) {
			rect[0] = mview->dragStartPos[0];
			rect[2] = mview->dragEndPos[0];
		} else {
			rect[0] = mview->dragEndPos[0];
			rect[2] = mview->dragStartPos[0];
		}
		if (mview->dragStartPos[1] < mview->dragEndPos[1]) {
			rect[1] = mview->dragStartPos[1];
			rect[3] = mview->dragEndPos[1];
		} else {
			rect[1] = mview->dragEndPos[1];
			rect[3] = mview->dragStartPos[1];
		}
		for (i = 0; i < natoms; i++) {
			ap = ATOM_AT_INDEX(mview->mol->atoms, i);
			if (ap == NULL)
				continue;
			if (MainView_isAtomHidden(mview, i))
				continue;
			if (ignoreExpandedAtoms && SYMOP_ALIVE(ap->symop))
				continue;
			if (mview->draggingMode == kMainViewSelectingRegion) {
				/*  Check if this atom is within the selection rectangle  */
				GLfloat objectPos[3];
				GLfloat screenPos[3];
				Vector r1;
				r1 = ap->r;
			/*	if (mview->mol->is_xtal_coord)
					TransformVec(&r1, mview->mol->cell->tr, &r1); */
				objectPos[0] = r1.x;
				objectPos[1] = r1.y;
				objectPos[2] = r1.z;
				if (MainView_convertObjectPositionToScreenPosition(mview, objectPos, screenPos) && screenPos[0] >= rect[0] && screenPos[0] <= rect[2] && screenPos[1] >= rect[1] && screenPos[1] <= rect[3])
					selectFlags[i] = 1;
				else selectFlags[i] = 0;
			}
		}
	}
	if (flags & kShiftKeyMask) {
		for (i = 0; i < natoms; i++)
			selectFlags[i] ^= (MoleculeIsAtomSelected(mview->mol, i) != 0);
	}
	return selectFlags;
}

static void
drawGraphite(MainView *mview)
{
	static GLfloat sDarkCyanColor[] = {0, 0.75, 0.75, 1};
	MDArena *arena;
	MDGraphiteArena *graphite;
	if ((arena = mview->mol->arena) != NULL && (graphite = arena->graphite) != NULL) {
		Vector xaxis, yaxis, zaxis, origin;
		Double R;
		int i, j, i0, i1, j0, j1, ir;
		Double x, dx, y, dy, xx, yy;
		GLfloat p[12];
		graphite_get_axes(graphite, &origin, &xaxis, &yaxis, &zaxis, &R);
		i0 = -(mview->showGraphite / 2) - 1;
		i1 = i0 + mview->showGraphite + 1;
		j0 = -(mview->showGraphite / 2);
		j1 = j0 + mview->showGraphite;
		dx = 0.5 * R;
		dy = 0.866025403784439 * R;
		glDisable(GL_LIGHTING);	
		glColor3fv(sDarkCyanColor);
		for (i = i0; i <= i1; i++) {
			for (j = j0; j <= j1; j++) {
				Byte f1, f2, f3;
				ir = (i % 2 == 0 ? 0 : 1);
				x = 3 * i * dx;
				y = (2 * j + ir) * dy;
				yy = y - dy;
				xx = x - 2 * dx;
				p[0] = xaxis.x * xx + yaxis.x * y + origin.x;
				p[1] = xaxis.y * xx + yaxis.y * y + origin.y;
				p[2] = xaxis.z * xx + yaxis.z * y + origin.z;
				xx += dx;
				p[3] = xaxis.x * xx + yaxis.x * yy + origin.x;
				p[4] = xaxis.y * xx + yaxis.y * yy + origin.y;
				p[5] = xaxis.z * xx + yaxis.z * yy + origin.z;
				xx += 2 * dx;
				p[6] = xaxis.x * xx + yaxis.x * yy + origin.x;
				p[7] = xaxis.y * xx + yaxis.y * yy + origin.y;
				p[8] = xaxis.z * xx + yaxis.z * yy + origin.z;
				xx += dx;
				p[9] = xaxis.x * xx + yaxis.x * y + origin.x;
				p[10] = xaxis.y * xx + yaxis.y * y + origin.y;
				p[11] = xaxis.z * xx + yaxis.z * y + origin.z;
				f1 = f2 = f3 = 1;
				if (i == i0) {
					f1 = f2 = 0;
					if ((ir == 0 && j == j0) || (ir == 1 && j == j1))
						continue;
				} else if (i == i1) {
					f2 = f3 = 0;
					if ((ir == 0 && j == j0) || (ir == 1 && j == j1))
						continue;
				} else if (j == j1) {
					if (ir == 1) {
						f1 = f3 = 0;
					} else if (i == i0 + 1) {
						f1 = 0;
					} else if (i == i1 - 1) {
						f3 = 0;
					}
				}
				glBegin(GL_LINES);
				if (f1) { 
					glVertex3fv(p);
					glVertex3fv(p + 3);
				}
				if (f2) {
					glVertex3fv(p + 3);
					glVertex3fv(p + 6);
				}
				if (f3) {
					glVertex3fv(p + 6);
					glVertex3fv(p + 9);
				}
				glEnd();
			}
		}
		glEnable(GL_LIGHTING);	
	}
}

static GLfloat sRedColor[] = {1, 0, 0, 1};

static void
drawAtom(MainView *mview, int i1, int selected, const Vector *dragOffset, const Vector *periodicOffset)
{
	const Atom *ap;
	const ElementPar *dp;
	int an1;
	int expanded = 0;
	Vector r1;
	GLfloat p[6];
	char label[16];
	GLfloat rgba[4];
	Transform *trp = NULL;
	int natoms = mview->mol->natoms;
	if (i1 >= natoms) {
		/*  Extra 2 atoms for the bond being newly created  */
		if (mview->draggingMode != kMainViewCreatingBond)
			return;
	/*	printf("mview->tempAtoms[%d] = %d\n", i - natoms, mview->tempAtoms[i - natoms]); */
		if (mview->tempAtoms[i1 - natoms] >= 0)
			return;  /*  Already drawn  */
		ap = NULL;
		an1 = 6;
		r1 = mview->tempAtomPos[i1 - natoms];
		label[0] = 0;
	} else if (i1 < 0) {
		ExAtom *ep = mview->mol->exatoms + (- i1 - 1);
		ap = ATOM_AT_INDEX(mview->mol->atoms, ep->index);
		an1 = ap->atomicNumber;
		r1 = ap->r;
		trp = &(mview->mol->syms[ep->symop]);
		if (/* !mview->mol->is_xtal_coord && */ mview->mol->cell != NULL) {
			TransformVec(&r1, mview->mol->cell->rtr, &r1);
			TransformVec(&r1, *trp, &r1);
			TransformVec(&r1, mview->mol->cell->tr, &r1);
		} else TransformVec(&r1, *trp, &r1);
		VecInc(r1, ep->dr);
		label[0] = 0;
	} else {
		ap = ATOM_AT_INDEX(mview->mol->atoms, i1);
		if (ap == NULL)
			return;
		an1 = ap->atomicNumber;
		r1 = ap->r;
		strncpy(label, ap->aname, 4);
		label[4] = 0;
		if (SYMOP_ALIVE(ap->symop))
			expanded = 1;
	}
	if (!mview->showHydrogens && an1 == 1)
		return;
	if (!mview->showDummyAtoms && an1 == 0)
		return;
	if (!mview->showExpandedAtoms && expanded)
		return;
	if (ap != NULL && (ap->exflags & kAtomHiddenFlag))
		return;
	dp = &(gElementParameters[an1]);
	if (dp == NULL)
		return;
	if (selected) {
		memcpy(rgba, sRedColor, sizeof(rgba));
	} else {
		rgba[0] = dp->r; rgba[1] = dp->g; rgba[2] = dp->b; rgba[3] = 1.0;
	}
	if (expanded || periodicOffset != NULL) {
		rgba[0] *= 0.5;
		rgba[1] *= 0.5;
		rgba[2] *= 0.5;
	}
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, rgba);
	if (periodicOffset != NULL)
		VecInc(r1, *periodicOffset);
/*	if (mview->mol->is_xtal_coord)
		TransformVec(&r1, mview->mol->cell->tr, &r1); */
	p[0] = r1.x; p[1] = r1.y; p[2] = r1.z;
	if (mview->draggingMode == kMainViewDraggingSelectedAtoms && selected) {
		p[0] += dragOffset->x;
		p[1] += dragOffset->y;
		p[2] += dragOffset->z;
	}
	if (mview->showEllipsoids) {
		if (ap != NULL && ap->aniso != NULL) {
		/*	Double xp[3][3]; */
		/*	int i, j; */
			GLfloat elip[9];
			Mat33 pmat2;
			int i;
			if (trp != NULL) {
				MatrixMul(pmat2, mview->mol->cell->rtr, ap->aniso->pmat);
				MatrixMul(pmat2, *((Mat33 *)trp), pmat2);
				MatrixMul(pmat2, mview->mol->cell->tr, pmat2);
				MatrixTranspose(pmat2, pmat2);
			} else {
				MatrixTranspose(pmat2, ap->aniso->pmat);
			}
			for (i = 0; i < 9; i++)
				elip[i] = pmat2[i] * mview->probabilityScale;
		/*	for (i = 0; i < 3; i++) {
				Double w = ap->aniso->val[i];
				Vector vv = ap->aniso->axis[i];
				if (w <= 0.0)
					w = 0.001;
				if (trp != NULL) {
					TransformVec(&vv, mview->mol->cell->rtr, &vv);
					MatrixVec(&vv, *((Mat33 *)trp), &vv);
					TransformVec(&vv, mview->mol->cell->tr, &vv);
				}
				w *= mview->probabilityScale;
				xp[i][0] = w * vv.x;
				xp[i][1] = w * vv.y;
				xp[i][2] = w * vv.z;
			} */
			
			drawEllipsoid(p, elip, elip+3, elip+6, 15);
		} else {
			Double rad;
			rad = biso2radius(ap->tempFactor);
			rad *= mview->probabilityScale;
			drawSphere(p, rad, 8);
		}
	} else {
		drawSphere(p, dp->radius * mview->atomRadius, 8);
	}
	if (MainView_convertObjectPositionToScreenPosition(mview, p, p + 3)) {
	/*	fprintf(stderr, "atom %d: {%f, %f, %f}\n", i1, p[3], p[4], p[5]); */
		float fp[3];
		fp[0] = p[3]; fp[1] = p[4]; fp[2] = p[5];
		MainViewCallback_drawLabel(mview, fp, label);
	}
}

static void
drawBond(MainView *mview, int i1, int i2, int selected, int selected2, int draft, const Vector *dragOffset, const Vector *periodicOffset)
{
	const ElementPar *dp;
	int i, in;
	int an[2];
	int expanded[2];
	Vector r[2];
	GLfloat p[6];
	GLfloat rgba[4];
	int natoms = mview->mol->natoms;
	expanded[0] = expanded[1] = 0;

	for (i = 0; i < 2; i++) {
		const Atom *ap;
		in = (i == 0 ? i1 : i2);
		if (in >= natoms && in < natoms + 2) {
			if (mview->tempAtoms[in - natoms] >= 0) {
				ap = ATOM_AT_INDEX(mview->mol->atoms, mview->tempAtoms[in - natoms]);
				an[i] = ap->atomicNumber;
				r[i] = ap->r;
			} else {
				ap = NULL;
				r[i] = mview->tempAtomPos[in - natoms];
				an[i] = 6;
			}
		} else if (in < 0) {
			ExAtom *ep = mview->mol->exatoms + (- in - 1);
			if (!mview->showExpandedAtoms)
				return;
			ap = ATOM_AT_INDEX(mview->mol->atoms, ep->index);
			an[i] = ap->atomicNumber;
			r[i] = ap->r;
			TransformVec(&r[i], mview->mol->syms[ep->symop], &r[i]);
			VecInc(r[i], ep->dr);
		} else {
			ap = ATOM_AT_INDEX(mview->mol->atoms, in);
			an[i] = ap->atomicNumber;
			r[i] = ap->r;
			if (SYMOP_ALIVE(ap->symop))
				expanded[i] = 1;
		}
		if (!mview->showHydrogens && an[i] == 1)
			return;
		if (!mview->showDummyAtoms && an[i] == 0)
			return;
		if (!mview->showExpandedAtoms && expanded[i])
			return;
		if (ap != NULL && (ap->exflags & kAtomHiddenFlag))
			return;
	}

	if (periodicOffset != NULL) {
		VecInc(r[0], *periodicOffset);
		VecInc(r[1], *periodicOffset);
	}
/*	if (mview->mol->is_xtal_coord) {
		TransformVec(&r[0], mview->mol->cell->tr, &r[0]);
		TransformVec(&r[1], mview->mol->cell->tr, &r[1]);
	} */

	dp = &(gElementParameters[an[0]]);
	if (dp == NULL)
		return;
	if (selected && selected2) {
		memcpy(rgba, sRedColor, sizeof(rgba));
	} else {
		rgba[0] = dp->r; rgba[1] = dp->g; rgba[2] = dp->b; rgba[3] = 1.0;
	}
	if (expanded[0] || periodicOffset != NULL) {
		rgba[0] *= 0.5;
		rgba[1] *= 0.5;
		rgba[2] *= 0.5;
	}		
	if (mview->draggingMode == kMainViewDraggingSelectedAtoms) {
		if (selected)
			VecInc(r[0], *dragOffset);
		if (selected2)
			VecInc(r[1], *dragOffset);
	}
	p[0] = r[0].x; p[1] = r[0].y; p[2] = r[0].z;
	p[3] = (r[1].x + p[0]) * 0.5;
	p[4] = (r[1].y + p[1]) * 0.5;
	p[5] = (r[1].z + p[2]) * 0.5;
	if (draft) {
		glColor3f(rgba[0], rgba[1], rgba[2]);
		glBegin(GL_LINES);
		glVertex3f(p[0], p[1], p[2]);
		glVertex3f(p[3], p[4], p[5]);
		glEnd();
	} else {
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, rgba);
		drawCylinder(p, p + 3, mview->bondRadius, 6);
	}
	dp = &(gElementParameters[an[1]]);
	if (dp == NULL)
		return;
	if (!selected || !selected2) {
		rgba[0] = dp->r; rgba[1] = dp->g; rgba[2] = dp->b; rgba[3] = 1.0;
	}
	if (expanded[1] || periodicOffset != NULL) {
		rgba[0] *= 0.5;
		rgba[1] *= 0.5;
		rgba[2] *= 0.5;
	}		
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, rgba);
	p[0] = r[1].x; p[1] = r[1].y; p[2] = r[1].z;
	if (draft) {
		glColor3f(rgba[0], rgba[1], rgba[2]);
		glBegin(GL_LINES);
		glVertex3f(p[0], p[1], p[2]);
		glVertex3f(p[3], p[4], p[5]);
		glEnd();
	} else {
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, rgba);
		drawCylinder(p, p + 3, mview->bondRadius, 6);
	}
}

/*  Calculate drag offset during moving the selection.  */
static void
calcDragOffset(MainView *mview, Vector *outVector)
{
	GLfloat p[6];
	if (mview->draggingMode == kMainViewDraggingSelectedAtoms
	&& MainView_convertScreenPositionToObjectPosition(mview, mview->dragStartPos, p)
	&& MainView_convertScreenPositionToObjectPosition(mview, mview->dragEndPos, p + 3)) {
		outVector->x = p[3] - p[0];
		outVector->y = p[4] - p[1];
		outVector->z = p[5] - p[2];
	} else {
		outVector->x = outVector->y = outVector->z = 0;
	}

}

static void
drawModel(MainView *mview)
{
	Molecule *mol;
    int i, natoms, nbonds;
	int amin, amax, bmin, bmax, cmin, cmax, da, db, dc;
	Byte original;
	Double atomRadius, bondRadius;
	Vector dragOffset;
	Vector periodicOffset;
	Vector *axes;
	int selected, selected2;
	char *selectFlags;
	int draft = mview->lineMode;
/*    static Double gray[] = {0.8, 0.8, 0.8, 1}; */
	
	atomRadius = mview->atomRadius;
	bondRadius = mview->bondRadius;
	mol = mview->mol;
	natoms = mol->natoms;
/*    if (natoms == 0) {
        return;
    }
*/
	if (mview->draggingMode == kMainViewSelectingRegion)
		selectFlags = temporarySelection(mview, mview->modifierFlags, 0, 1);
	else selectFlags = NULL;
	
	if (mview->draggingMode == kMainViewDraggingSelectedAtoms)
		calcDragOffset(mview, &dragOffset);
	else dragOffset.x = dragOffset.y = dragOffset.z = 0;
	
	if (mview->showGraphite) {
		drawGraphite(mview);
	}
	
	amin = amax = bmin = bmax = cmin = cmax = 0;
	if (mview->showExpandedAtoms && mol->cell != NULL) {
		if (mol->cell->flags[0]) {
			amin = mview->showPeriodicImage[0];
			amax = mview->showPeriodicImage[1];
		}
		if (mol->cell->flags[1]) {
			bmin = mview->showPeriodicImage[2];
			bmax = mview->showPeriodicImage[3];
		}
		if (mol->cell->flags[2]) {
			cmin = mview->showPeriodicImage[4];
			cmax = mview->showPeriodicImage[5];
		}
		axes = mol->cell->axes;
	} else {
		axes = NULL;
	}
	
	if (draft == 0) {
		for (da = amin; da <= amax; da++) {
			for (db = bmin; db <= bmax; db++) {
				for (dc = cmin; dc <= cmax; dc++) {
					original = (da == 0 && db == 0 && dc == 0);
					if (!original) {
						VecScale(periodicOffset, axes[0], da);
						VecScaleInc(periodicOffset, axes[1], db);
						VecScaleInc(periodicOffset, axes[2], dc);
					}
					for (i = 0; i < natoms; i++) {
						if (mview->draggingMode != 0 && i % 50 == 0 && MainViewCallback_mouseCheck(mview)) {
							/*  Mouse event is detected  */
							draft = 1;
							goto skip;
						}
						if (mview->draggingMode == kMainViewCreatingBond && (i == mview->tempAtoms[0] || i == mview->tempAtoms[1] || i >= natoms))
							selected = 1;  /*  extra atoms  */
						else if (selectFlags != NULL)
							selected = selectFlags[i];
						else
							selected = MoleculeIsAtomSelected(mview->mol, i);
						drawAtom(mview, i, selected, &dragOffset, (original ? NULL : &periodicOffset));
					}
				}
	
				if (draft == 0) {
					if (original) {
						/*  Extra atoms  */
						drawAtom(mview, natoms, 1, &dragOffset, NULL);
						drawAtom(mview, natoms + 1, 1, &dragOffset, NULL);
					}
					/*  Expanded atoms  */
					if (mview->showExpandedAtoms) {
						for (i = 0; i < mview->mol->nexatoms; i++) {
							if (mview->draggingMode != 0 && i % 50 == 0 && MainViewCallback_mouseCheck(mview)) {
								/*  Mouse event is detected  */
								draft = 1;
								break;
							/*	goto cleanup;  */
							}
							drawAtom(mview, -i-1, (selectFlags == NULL ? 0 : selectFlags[i]), &dragOffset, (original ? NULL : &periodicOffset));
						}
					}
				}
			}
		}
	}
	
skip:
	nbonds = mol->nbonds;
	if (draft)
		glDisable(GL_LIGHTING);	
	for (da = amin; da <= amax; da++) {
		for (db = bmin; db <= bmax; db++) {
			for (dc = cmin; dc <= cmax; dc++) {
				original = (da == 0 && db == 0 && dc == 0);
				if (!original) {
					VecScale(periodicOffset, axes[0], da);
					VecScaleInc(periodicOffset, axes[1], db);
					VecScaleInc(periodicOffset, axes[2], dc);
				}
				
				for (i = 0; i < nbonds; i++) {
					int n1, n2;
					if (draft == 0 && mview->draggingMode != 0 && i % 50 == 0 && MainViewCallback_mouseCheck(mview)) {
						/*  Mouse event is detected  */
						draft = 1;
						glDisable(GL_LIGHTING);
					/*	goto cleanup;  */
					}
					n1 = mview->mol->bonds[i * 2];
					n2 = mview->mol->bonds[i * 2 + 1];
					if (selectFlags == NULL) {
						selected = MoleculeIsAtomSelected(mview->mol, n1);
						selected2 = MoleculeIsAtomSelected(mview->mol, n2);
					} else {
						selected = selectFlags[n1];
						selected2 = selectFlags[n2];
					}
					drawBond(mview, n1, n2, selected, selected2, draft, &dragOffset, (original ? NULL : &periodicOffset));
				}
				
				/*  Extra bond  */
				if (original && mview->draggingMode == kMainViewCreatingBond) {
					drawBond(mview, natoms, natoms + 1, 1, 1, draft, &dragOffset, NULL);
				}
				
				/*  Expanded bonds  */
				for (i = 0; i < mview->mol->nexbonds; i++) {
					int n1, n2;
					if (draft == 0 && mview->draggingMode != 0 && i % 50 == 0 && MainViewCallback_mouseCheck(mview)) {
						/*  Mouse event is detected  */
						draft = 1;
						glDisable(GL_LIGHTING);
					/*	goto cleanup;  */
					}
					n1 = mview->mol->exbonds[i * 2];
					n2 = mview->mol->exbonds[i * 2 + 1];
					if (n1 < 0)
						n1 = -n1 - 1;
					if (n2 < 0)
						n2 = -n2 - 1;
					if (selectFlags == NULL) {
						selected = MoleculeIsAtomSelected(mview->mol, n1);
						selected2 = MoleculeIsAtomSelected(mview->mol, n2);
					} else {
						selected = selectFlags[n1];
						selected2 = selectFlags[n2];
					}
					drawBond(mview, mview->mol->exbonds[i * 2], mview->mol->exbonds[i * 2 + 1], selected, selected2, draft, &dragOffset, (original ? NULL : &periodicOffset));
				}
			}
		}
	}
	
/*  cleanup: */
	if (draft)
		glEnable(GL_LIGHTING);
	if (selectFlags != NULL)
		free(selectFlags);
}

static void
drawUnitCell(MainView *mview)
{
	GLfloat a[3], b[3], c[3], ab[3], bc[3], ca[3], abc[3];
	XtalCell *cp;
	static GLfloat sGrayColor[] = {0.5, 0.5, 0.5, 1};
	GLfloat origin[3];
	int i;
	glDisable(GL_LIGHTING);
	for (i = 0; i < 1; i++) {
		if (i == 0) {
			if (!mview->showUnitCell || (cp = mview->mol->cell) == NULL)
				continue;
			a[0] = cp->tr[0];
			a[1] = cp->tr[3];
			a[2] = cp->tr[6];
			b[0] = cp->tr[1];
			b[1] = cp->tr[4];
			b[2] = cp->tr[7];
			c[0] = cp->tr[2];
			c[1] = cp->tr[5];
			c[2] = cp->tr[8];
			origin[0] = cp->origin.x;
			origin[1] = cp->origin.y;
			origin[2] = cp->origin.z;
			glColor3f(0.75, 0.2, 0.0);
/*		} else {
			if (!mview->showPeriodicBox || (bp = mview->mol->box) == NULL)
				continue;
			origin[0] = bp->origin.x;
			origin[1] = bp->origin.y;
			origin[2] = bp->origin.z;
			a[0] = bp->axes[0].x + origin[0];
			a[1] = bp->axes[0].y + origin[1];
			a[2] = bp->axes[0].z + origin[2];
			b[0] = bp->axes[1].x + origin[0];
			b[1] = bp->axes[1].y + origin[1];
			b[2] = bp->axes[1].z + origin[2];
			c[0] = bp->axes[2].x + origin[0];
			c[1] = bp->axes[2].y + origin[1];
			c[2] = bp->axes[2].z + origin[2];
			glColor3f(0.0, 0.75, 0.2); */
		}
	
		ab[0] = a[0] + b[0]; ab[1] = a[1] + b[1]; ab[2] = a[2] + b[2];
		bc[0] = b[0] + c[0]; bc[1] = b[1] + c[1]; bc[2] = b[2] + c[2];
		ca[0] = c[0] + a[0]; ca[1] = c[1] + a[1]; ca[2] = c[2] + a[2];
		abc[0] = a[0] + bc[0]; abc[1] = a[1] + bc[1]; abc[2] = a[2] + bc[2];

	/*	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, sGrayColor);
		drawCylinder(sOrigin, a, 0.05, 4);
		drawCylinder(sOrigin, b, 0.05, 4);
		drawCylinder(sOrigin, c, 0.05, 4);
		drawCylinder(a, ab, 0.05, 4);
		drawCylinder(a, ca, 0.05, 4);
		drawCylinder(b, ab, 0.05, 4);
		drawCylinder(b, bc, 0.05, 4);
		drawCylinder(c, ca, 0.05, 4);
		drawCylinder(c, bc, 0.05, 4);
		drawCylinder(ab, abc, 0.05, 4);
		drawCylinder(bc, abc, 0.05, 4);
		drawCylinder(ca, abc, 0.05, 4); */

		glBegin(GL_LINES);
		glVertex3fv(origin);
		glVertex3fv(a);
		glVertex3fv(origin);
		glVertex3fv(b);
		glVertex3fv(origin);
		glVertex3fv(c);
		glVertex3fv(a);
		glVertex3fv(ab);
		glVertex3fv(a);
		glVertex3fv(ca);
		glVertex3fv(b);
		glVertex3fv(ab);
		glVertex3fv(b);
		glVertex3fv(bc);
		glVertex3fv(c);
		glVertex3fv(ca);
		glVertex3fv(c);
		glVertex3fv(bc);
		glVertex3fv(ab);
		glVertex3fv(abc);
		glVertex3fv(bc);
		glVertex3fv(abc);
		glVertex3fv(ca);
		glVertex3fv(abc);
		glEnd();
	}
	glEnable(GL_LIGHTING);
}

static void
drawRotationCenter(MainView *mview)
{
	GLfloat ps[3], pe[3], col[4];
	float fd[2];  /*  Fovy and distance  */
	float tr[3];  /*  Translation  */
	float r, rr;
	if (mview == NULL || !mview->showRotationCenter)
		return;
	TrackballGetTranslate(mview->track, tr);
	TrackballGetPerspective(mview->track, fd);
	tr[0] *= -mview->dimension;
	tr[1] *= -mview->dimension;
	tr[2] *= -mview->dimension;
	r = fd[1] * mview->dimension * tan(fd[0] * 0.5 * kDeg2Rad) * 0.1;
	rr = r * 0.1;
	ps[0] = tr[0];
	ps[1] = tr[1];
	ps[2] = tr[2];
	col[0] = col[1] = col[2] = 0.5; col[3] = 1.0;
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
	drawSphere(ps, rr, 8);
	ps[0] = tr[0] - r;
	pe[0] = tr[0] + r;
	pe[1] = tr[1];
	pe[2] = tr[2];
	col[0] = 1.0; col[1] = col[2] = 0.0;
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
	drawSphere(ps, rr, 8);
	drawSphere(pe, rr, 8);
	drawCylinder(ps, pe, rr, 8);
	ps[0] = tr[0];
	ps[1] = tr[1] - r;
	pe[0] = tr[0];
	pe[1] = tr[1] + r;
	col[1] = 1.0; col[0] = col[2] = 0.0;
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
	drawSphere(ps, rr, 8);
	drawSphere(pe, rr, 8);
	drawCylinder(ps, pe, rr, 8);
	ps[1] = tr[1];
	ps[2] = tr[2] - r;
	pe[1] = tr[1];
	pe[2] = tr[2] + r;
	col[2] = 1.0; col[0] = col[1] = 0.0;
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
	drawSphere(ps, rr, 8);
	drawSphere(pe, rr, 8);
	drawCylinder(ps, pe, rr, 8);
}

static int
compareLabelByDepth(const void *ap, const void *bp)
{
	Double dz = ((const LabelRecord *)bp)->pos.z - ((const LabelRecord *)ap)->pos.z;
	if (dz > 0)
		return 1;
	else if (dz < 0)
		return -1;
	else return 0;
}

static void
drawLabels(MainView *mview)
{
	Transform *trp;
	Atom *ap;
	LabelRecord *lp;
	int i, nlabels;

	if (mview->nlabels == 0)
		return;
	
/*	mview->sortedLabels = (LabelRecord **)calloc(sizeof(LabelRecord *), mview->nlabels);
	if (mview->sortedLabels == NULL)
		return; */
	
/*	if (mview->mol->is_xtal_coord)
		trp = &(mview->mol->cell->tr);
	else trp = NULL; */

	/*  Get the screen coordinates of the labels  */
	nlabels = 0;
	for (i = 0; i < mview->mol->natoms; i++) {
		float scrp[3], f[3];
		ap = ATOM_AT_INDEX(mview->mol->atoms, i);
		if (ap->exflags & kAtomHiddenFlag)
			continue;
		if (!mview->showHydrogens && ap->atomicNumber == 1)
			continue;
		if (!mview->showDummyAtoms && ap->atomicNumber == 0)
			continue;
		if (ap->labelid <= 0 || ap->labelid > mview->nlabels)
			continue;
		lp = mview->labels + (ap->labelid - 1);
		MainView_screenCenterPointOfAtom(mview, lp->idx1, scrp);
	/*	if (lp->idx2 >= 0) {
			Atom *bp;
			r1 = ap->r;
			if (trp != NULL)
				TransformVec(&r1, *trp, &r1);
			bp = ATOM_AT_INDEX(mview->mol->atoms, lp->idx2);
			if (bp->exflags & kAtomHiddenFlag)
				continue;
			r2 = bp->r;
			if (trp != NULL)
				TransformVec(&r2, *trp, &r2);
			r1.x = (r1.x + r2.x) * 0.5;
			r1.y = (r1.y + r2.y) * 0.5;
			r1.z = (r1.z + r2.z) * 0.5;
			objp[0] = r1.x;
			objp[1] = r1.y;
			objp[2] = r1.z;
			MainView_convertObjectPositionToScreenPosition(mview, objp, scrp);
		} */
		lp->pos.x = scrp[0];
		lp->pos.y = scrp[1];
		lp->pos.z = scrp[2];
		MainViewCallback_labelSize(lp->label, f);
		f[0] = floor(lp->pos.x - f[0] * 0.5);
		f[1] = -floor(lp->pos.y + f[1] * 0.5);
		f[2] = lp->pos.z;
	/*	fprintf(stderr, "label position (%d) = {%f, %f, %f}\n", i, f[0], f[1], f[2]); */
		MainViewCallback_drawLabelAtPoint(lp->label, f);
	/*	mview->sortedLabels[nlabels++] = lp;
		if (nlabels >= mview->nlabels)
			break; */
	}
	
/*	//  Sort labels by z coordinates (descending order) 
	qsort(mview->sortedLabels, nlabels, sizeof(LabelRecord *), compareLabelByDepth);
	
	//  Draw labels 
	for (i = 0; i < nlabels; i++) {
		LabelRecord *lp = mview->sortedLabels[i];
		Double f[3];
		MainViewCallback_labelSize(lp->label, f);
		f[0] = floor(lp->pos.x - f[0] * 0.5);
		f[1] = -floor(lp->pos.y + f[1] * 0.5);
		f[2] = lp->pos.z;
		MainViewCallback_drawLabelAtPoint(lp->label, f);
	}
*/	
}

void
MainView_drawModel(MainView *mview)
{
    float w[4], dimension, distance;
	float frame[4], width, height;
	GLdouble *pp;
	Transform mtr;
	
/*    int i;  */
	
	if (mview == NULL)
		return;

    if (!mview->isInitialized) {
        MainView_initializeOpenGLView(mview);
    }
    
    /*  Clear the buffer  */
    glClearColor (0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT);

	if (mview->mol == NULL)
		return;

	dimension = mview->dimension;
/*    dimension = [model dimension];  */

	MainViewCallback_frame(mview, frame);
	width = frame[2] - frame[0];
	height = frame[3] - frame[1];

    glViewport(0, 0, width, height);
    
    /*  Set up the projection  */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	TrackballGetPerspective(mview->track, w);
    distance = w[1] * dimension;
    mview->perspective_vector[0] = w[0];
    mview->perspective_vector[1] = width / height;
    mview->perspective_vector[2] = dimension;
    mview->perspective_vector[3] = distance + 200.0 * dimension;
    gluPerspective(mview->perspective_vector[0], mview->perspective_vector[1], mview->perspective_vector[2], mview->perspective_vector[3]);

    /*  Set up the model view  */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -distance);
	TrackballGetRotate(mview->track, w);
    glRotatef(w[0], w[1], w[2], w[3]);
	TrackballGetTranslate(mview->track, w);
	w[0] *= dimension;
	w[1] *= dimension;
	w[2] *= dimension;
    glTranslatef(w[0], w[1], w[2]);
	mview->lookat.x = -w[0];
	mview->lookat.y = -w[1];
	mview->lookat.z = -w[2];
	
	MainViewCallback_clearLabels(mview);
    drawModel(mview);
	drawUnitCell(mview);
	drawRotationCenter(mview);

	/*  Get important matrices and vectors  */
    glGetDoublev(GL_MODELVIEW_MATRIX, mview->modelview_matrix);
    glGetDoublev(GL_PROJECTION_MATRIX, mview->projection_matrix);
	pp = mview->modelview_matrix;
	mtr[0] = pp[0]; mtr[1] = pp[4]; mtr[2] = pp[8];
	mtr[3] = pp[1]; mtr[4] = pp[5]; mtr[5] = pp[9];
	mtr[6] = pp[2]; mtr[7] = pp[6]; mtr[8] = pp[10];
	mtr[9] = pp[12]; mtr[10] = pp[13]; mtr[11] = pp[14];
	TransformInvert(mtr, mtr);
	mview->camera.x = mtr[9]; mview->camera.y = mtr[10]; mview->camera.z = mtr[11];
	mview->lookto.x = mtr[2]; mview->lookto.y = mtr[5]; mview->lookto.z = mtr[8];
	mview->up.x = mtr[1]; mview->up.y = mtr[4]; mview->up.z = mtr[7];

	/*  Draw labels  */
	glDisable(GL_LIGHTING);
//	glDisable (GL_DEPTH_TEST);
//	glEnable (GL_BLEND);
//	glBlendFunc (GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
#if __WXMAC__
	glEnable (GL_TEXTURE_RECTANGLE_EXT);
#endif
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	glOrtho(0, width, 0, -height, 0.0,-1.0);  /*  non-flipped view  */
	drawLabels(mview);
#if __WXMAC__
	glDisable (GL_TEXTURE_RECTANGLE_EXT);
#endif
//	glDisable(GL_BLEND);
//	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

    if (mview->draggingMode == kMainViewSelectingRegion) {
		/*  Draw selection rectangle  */
		glDisable(GL_LIGHTING);
		glDisable(GL_DEPTH_TEST);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, width, 0, height, -1.0, 1.0);
        glColor3f(1.0, 1.0, 0.0);
        glBegin(GL_LINE_STRIP);
		glVertex2f(mview->dragStartPos[0], mview->dragStartPos[1]);
		glVertex2f(mview->dragStartPos[0], mview->dragEndPos[1]);
		glVertex2f(mview->dragEndPos[0], mview->dragEndPos[1]);
		glVertex2f(mview->dragEndPos[0], mview->dragStartPos[1]);
		glVertex2f(mview->dragStartPos[0], mview->dragStartPos[1]);
        glEnd();
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
    }

    glFinish();
        
}

#pragma mark ====== Labels ======

void
MainView_attachLabelToAtom(MainView *mview, int index)
{
	LabelRecord rec;
	char buf[24], *p;
	Atom *ap;
/*	const ElementPar *dp; */
	static float foreColor[] = {1, 1, 0, 1};
	static float backColor[] = {1, 1, 1, 0};
	ap = ATOM_AT_INDEX(mview->mol->atoms, index);
	if (ap->resSeq == 0)
		snprintf(buf, sizeof buf, "%-.4s", ap->aname);
	else
		snprintf(buf, sizeof buf, "%d:%-.4s", ap->resSeq, ap->aname);
	for (p = buf; *p; p++) {
		if (isspace(*p)) {
			*p = 0;
			break;
		}
	}
/*	dp = &(gBuiltinParameters->atomPars[ap->atomicNumber]);
	foreColor[0] = 1.0 - dp->r;
	foreColor[1] = 1.0 - dp->g;
	foreColor[2] = 1.0 - dp->b; */
	rec.label = MainViewCallback_newLabel(mview, buf, 14, foreColor, backColor);
	rec.idx1 = index;
	rec.idx2 = kInvalidIndex;
	rec.labelid = mview->nlabels + 1;
	ap->labelid = rec.labelid;
	AssignArray(&mview->labels, &mview->nlabels, sizeof(LabelRecord), mview->nlabels, &rec);
}

void
MainView_detachLabelFromAtom(MainView *mview, int index)
{
	Atom *ap;
	if (index >= 0 && index < mview->mol->natoms) {
		ap = ATOM_AT_INDEX(mview->mol->atoms, index);
		if (ap->labelid > 0 && ap->labelid <= mview->nlabels) {
			mview->labels[ap->labelid - 1].idx1 = kInvalidIndex;
		}
		ap->labelid = 0;
	}
}

void
MainView_purgeUnusedLabels(MainView *mview)
{
	int i, n;
	Atom *ap;
	Int *tempid;
	if (mview == NULL || mview->nlabels == 0 || mview->mol->natoms == 0)
		return;
	tempid = (Int *)calloc(sizeof(Int), mview->nlabels);
	if (tempid == NULL)
		return;

	/*  Mark the labels in use  */
	for (i = 0, ap = mview->mol->atoms; i < mview->mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->labelid > 0 && ap->labelid <= mview->nlabels)
			tempid[ap->labelid - 1] = 1;
	}

	/*  Release the labels not in use, and assign new IDs  */
	n = 1;
	for (i = 0; i < mview->nlabels; i++) {
		if (tempid[i] == 0) {
			MainViewCallback_releaseLabel(mview->labels[i].label);
			mview->labels[i].label = NULL;
		} else {
			mview->labels[i].labelid = tempid[i] = n++;
		}
	}
	
	/*  Purge the unused entries  */
	for (i = mview->nlabels - 1; i >= 0; i--) {
		if (tempid[i] == 0) {
			memmove(mview->labels + i + 1, mview->labels + i, (mview->nlabels - 1 - i) * sizeof(LabelRecord));
			mview->nlabels--;
		}
	}
	if (mview->nlabels == 0) {
		free(mview->labels);
		mview->labels = NULL;
	}
	
	/*  Renumber  */
	for (i = 0, ap = mview->mol->atoms; i < mview->mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->labelid > 0 && ap->labelid <= mview->nlabels) {
			ap->labelid = tempid[ap->labelid - 1];
		} else {
			ap->labelid = 0;
		}
	}
}

#pragma mark ====== Mode ======

void
MainView_setMode(MainView *mview, int mode)
{
	if (mview != NULL)
		mview->mode = mode;
}

int
MainView_getMode(const MainView *mview)
{
	if (mview != NULL)
		return mview->mode;
	else return 0;
}

#pragma mark ====== Mouse operations ======

static void
mousePosToTrackballPos(MainView *mview, const float *mousePos, float *pos)
{
	float frame[4], width, height, radius;
	NULL_CHECK(mview, "mousePosToTrackballPos");
	MainViewCallback_frame(mview, frame);
	width = frame[2] - frame[0];
	height = frame[3] - frame[1];
	radius = (width > height ? height * 0.5 : width * 0.5);
	pos[0] = (mousePos[0] - frame[0] - width * 0.5) / radius;
	pos[1] = (mousePos[1] - frame[1] - height * 0.5) / radius;
}

static void
showAtomsInInfoText(MainView *mview, int n1, int n2)
{
	char buf[64];
	if (n1 >= 0) {
		MoleculeGetAtomName(mview->mol, n1, buf, sizeof buf - 1);
		if (n2 >= 0) {
			int nn = strlen(buf);
			buf[nn++] = '-';
			MoleculeGetAtomName(mview->mol, n2, buf + nn, sizeof buf - nn);
		}
		MainViewCallback_drawInfoText(mview, buf);
	} else MainViewCallback_drawInfoText(mview, NULL);
}

void
MainView_mouseDown(MainView *mview, const float *mousePos, int flags)
{
	float p[3];
	float screenPos[3], objectPos[3];
	int n1, n2, found;
	Atom *ap;
	mview->dragStartPos[0] = mousePos[0];
	mview->dragStartPos[1] = mousePos[1];
	mview->dragStartPos[2] = 0.5;
	found = MainView_findObjectAtPoint(mview, mousePos, &n1, &n2, 0, 0);
	mview->clickedAtoms[0] = n1;
	mview->clickedAtoms[1] = n2;
	if (found) {
		/*  Estimate the screen z-coordinate of the mouse position  */
		Vector r1, r2;
		ap = ATOM_AT_INDEX(mview->mol->atoms, n1);
		r1 = ap->r;
	/*	if (mview->mol->is_xtal_coord)
				TransformVec(&r1, mview->mol->cell->tr, &r1); */
		objectPos[0] = r1.x;
		objectPos[1] = r1.y;
		objectPos[2] = r1.z;
		if (MainView_convertObjectPositionToScreenPosition(mview, objectPos, screenPos)) {
			if (n2 >= 0) {
				ap = ATOM_AT_INDEX(mview->mol->atoms, n2);
				r2 = ap->r;
			/*	if (mview->mol->is_xtal_coord)
						TransformVec(&r2, mview->mol->cell->tr, &r2); */
				objectPos[0] = r2.x;
				objectPos[1] = r2.y;
				objectPos[2] = r2.z;
				if (MainView_convertObjectPositionToScreenPosition(mview, objectPos, p)) {
					Double w;
					r1.x = p[0] - screenPos[0];
					r1.y = p[1] - screenPos[1];
					r1.z = p[2] - screenPos[2];
					NormalizeVec(&r1, &r1);
					r2.x = mousePos[0] - screenPos[0];
					r2.y = mousePos[1] - screenPos[1];
					r2.z = 0.5;
					w = VecDot(r1, r2);
					VecScale(r2, r1, w);
					screenPos[2] += r2.z;
				}
			}
		}
		mview->dragStartPos[2] = screenPos[2];
	} else {
		/*  Set the 'depth' value of the molecule center (i.e. trackball->translate * dimension)
		    as the z-coordinate of the drag start position  */
		TrackballGetTranslate(mview->track, p);
		objectPos[0] = -mview->dimension * p[0];
		objectPos[1] = -mview->dimension * p[1];
		objectPos[2] = -mview->dimension * p[2];
		if (MainView_convertObjectPositionToScreenPosition(mview, objectPos, screenPos))
			mview->dragStartPos[2] = screenPos[2];
	}

	mview->isDragging = 0;

	switch (mview->mode) {
		case kTrackballRotateMode:
		case kTrackballTranslateMode:
		case kTrackballScaleMode:
			if (!found) {
				mousePosToTrackballPos(mview, mousePos, p);
				TrackballStartDragging(mview->track, p, mview->mode);
				mview->draggingMode = kMainViewMovingTrackball;
				break;
			} else {
				/*  No 'break' intentional; drop to next case  */
			}
		case kTrackballSelectionMode: {
			/*  Select or move around */
			if (found) {
				if (flags & kShiftKeyMask) {
					/*  Shift-key pressed: toggle selection  */
					MoleculeToggleSelectionOfAtom(mview->mol, n1);
					if (n2 >= 0)
						MoleculeToggleSelectionOfAtom(mview->mol, n2);
				} else {
					if (n2 < 0) {
						if (!MoleculeIsAtomSelected(mview->mol, n1)) {
							/*  If this atom is not selected, select this atom and unselect others  */
							MoleculeSelectAtom(mview->mol, n1, 0);
						}
					} else {
						if (!MoleculeIsBondSelected(mview->mol, n1, n2)) {
							/*  If this bond is not selected, select this bond and unselect others */
							MoleculeSelectAtom(mview->mol, n1, 0);
							MoleculeSelectAtom(mview->mol, n2, 1);
						}
					}
				}
				if (MoleculeIsAtomSelected(mview->mol, n1))
					mview->draggingMode = kMainViewDraggingSelectedAtoms;
				else mview->draggingMode = 0;
			} else mview->draggingMode = kMainViewSelectingRegion;
			mview->dragEndPos[0] = mousePos[0];
			mview->dragEndPos[1] = mousePos[1];
			mview->dragEndPos[2] = mview->dragStartPos[2];
			mview->modifierFlags = flags;
			break;
		}
		case kTrackballCreateMode: {
			if (md_is_running(mview->mol->arena)) {
				MoleculeCallback_cannotModifyMoleculeDuringMDError(mview->mol);
				mview->draggingMode = 0;
				break;
			}
			/*  Draw a new bond  */
			MoleculeSetSelection(mview->mol, NULL);
			if (found && n2 < 0) {
				/*  An atom under mouse: create a new atom and a bond  */
				ap = ATOM_AT_INDEX(mview->mol->atoms, n1);
				mview->tempAtoms[0] = n1;
				mview->tempAtomPos[0] = ap->r;
			} else {
				/*  No atom under mouse: create two new atoms and a bond  */
				mview->tempAtoms[0] = -1;
				screenPos[0] = mousePos[0];
				screenPos[1] = mousePos[1];
				screenPos[2] = mview->dragStartPos[2];
				if (MainView_convertScreenPositionToObjectPosition(mview, screenPos, objectPos) == 0) {
					mview->tempAtoms[0] = -1;  /*  Cannot create  */
					mview->draggingMode = 0;
				} else {
					mview->tempAtomPos[0].x = objectPos[0];
					mview->tempAtomPos[0].y = objectPos[1];
					mview->tempAtomPos[0].z = objectPos[2];
				/*	if (mview->mol->is_xtal_coord)
						TransformVec(&mview->tempAtomPos[0], mview->mol->cell->rtr, &mview->tempAtomPos[0]); */
				}
			}
			mview->tempAtoms[1] = -1;
			mview->tempAtomPos[1] = mview->tempAtomPos[0];
			mview->draggingMode = kMainViewCreatingBond;
			break;
		}
		default:
			mview->draggingMode = 0;
			break;
	}
	if (found)
		showAtomsInInfoText(mview, n1, n2);
	else
		showAtomsInInfoText(mview, -1, -1);
/*	MainViewCallback_setNeedsDisplay(mview, 1); */
}

void
MainView_mouseDragged(MainView *mview, const float *mousePos, int flags)
{
	float p[2];
	if (mview->isDragging == 0) {
		if (abs(mousePos[0] - mview->dragStartPos[0]) >= 3 || abs(mousePos[1] - mview->dragStartPos[1]) >= 3)
			mview->isDragging = 1;
		else return;
	}
	mousePosToTrackballPos(mview, mousePos, p);
	switch (mview->draggingMode) {
		case kMainViewMovingTrackball:
			TrackballDrag(mview->track, p);
			break;
		case kMainViewSelectingRegion:
		case kMainViewDraggingSelectedAtoms:
			mview->dragEndPos[0] = mousePos[0];
			mview->dragEndPos[1] = mousePos[1];
			mview->dragEndPos[2] = mview->dragStartPos[2];
			mview->modifierFlags = flags;
			MainViewCallback_display(mview);
			break;
		case kMainViewCreatingBond: {
			int n1, n2;
			if (MainView_findObjectAtPoint(mview, mousePos, &n1, &n2, 0, 0) && n2 < 0) {
				/*  An atom under mouse  */
				Atom *ap = ATOM_AT_INDEX(mview->mol->atoms, n1);
			/*	printf("n1=%d, n2=%d\n", n1, n2); */
				mview->tempAtoms[1] = n1;
				mview->tempAtomPos[1] = ap->r;
			} else {
				float screenPos[3], objectPos[3];
				Vector r1;
				mview->tempAtoms[1] = -1;
				/*  Convert the position of temporary atom 0 to screen position */
				r1 = mview->tempAtomPos[0];
			/*	if (mview->mol->is_xtal_coord)
					TransformVec(&r1, mview->mol->cell->tr, &r1); */
				objectPos[0] = r1.x;
				objectPos[1] = r1.y;
				objectPos[2] = r1.z;
				if (MainView_convertObjectPositionToScreenPosition(mview, objectPos, screenPos) == 0)
					break;  /*  Do nothing  */
				/*  Convert the mouse position to object position, while using the same Z depth calculated from the temporary atom 0 position  */
				screenPos[0] = mousePos[0];
				screenPos[1] = mousePos[1];
				if (MainView_convertScreenPositionToObjectPosition(mview, screenPos, objectPos) == 0)
					break;  /*  Do nothing  */
				r1.x = objectPos[0];
				r1.y = objectPos[1];
				r1.z = objectPos[2];
			/*	if (mview->mol->is_xtal_coord)
					TransformVec(&r1, mview->mol->cell->rtr, &r1); */
				mview->tempAtomPos[1] = r1;
			}
			if (mview->tempAtoms[0] < 0)
				showAtomsInInfoText(mview, mview->tempAtoms[1], -1);
			else
				showAtomsInInfoText(mview, mview->tempAtoms[0], mview->tempAtoms[1]);
			MainViewCallback_display(mview);
			break;
		}
		default:
			return;
	}
	MainViewCallback_display(mview);
}

void
MainView_mouseMoved(MainView *mview, const float *mousePos, int flags)
{
	int n1, n2, found;
	if (mview->isDragging)
		return;
	found = MainView_findObjectAtPoint(mview, mousePos, &n1, &n2, 1, 0);
	if (found)
		showAtomsInInfoText(mview, n1, n2);
	else
		showAtomsInInfoText(mview, -1, -1);
}

static int
sCreateNewAtom(Molecule *mol, Vector pos)
{
	Int i, j;
	Atom a;
	char name[6];
	memset(&a, 0, sizeof(Atom));
	a.occupancy = 1.0;
	for (i = 0; i < 1000; i++) {
		snprintf(name, sizeof name, "C%03d", i);
		for (j = 0; j < mol->natoms; j++) {
			if (strncmp(ATOM_AT_INDEX(mol->atoms, j)->aname, name, 4) == 0)
				break;
		}
		if (j == mol->natoms)
			break;
	}
	strncpy(a.aname, name, 4);
	a.atomicNumber = 6;
	strcpy(a.element, "C");
	a.type = AtomTypeEncodeToUInt("c3");
	a.weight = WeightForAtomicNumber(a.atomicNumber);
	a.r = pos;
	if (MolActionCreateAndPerform(mol, gMolActionAddAnAtom, &a, -1, &i) == 0)
		return mol->natoms - 1;
	else return -1;
/*	return MoleculeAddAtom(mol, &a); */
}

void
MainView_mouseUp(MainView *mview, const float *mousePos, int flags, int clickCount)
{
	float p[2];
	char buf[1024];

	mousePosToTrackballPos(mview, mousePos, p);

	if (clickCount == 2 && mview->isDragging == 0) {
		/*  Create a new molecular fragment  */
		buf[0] = 0;
		if (MyAppCallback_getTextWithPrompt("Enter formula (e.g. CH2OCH3)", buf, sizeof buf) == 0)
			return;
		MolActionCreateAndPerform(mview->mol, SCRIPT_ACTION("s"), "dock_formula", buf);
		goto exit;
	}

	if (mview->mode == kTrackballEraseMode) {
		int n1, n2, found;
		found = MainView_findObjectAtPoint(mview, mousePos, &n1, &n2, 0, 0);
		if (found && n1 == mview->clickedAtoms[0] && n2 == mview->clickedAtoms[1]) {
			if (n2 < 0) {
			/*	IntGroup *sel, *newsel;
				sel = MoleculeGetSelection(mview->mol);
				if (sel != NULL) {
					newsel = MoleculeModifySelectionByRemovingAtoms(mview->mol, sel, IntGroupNewWithPoints(n1, 1, -1));
					if (!IntGroupIsEqual(sel, newsel))
						MolActionCreateAndPerform(mview->mol, gMolActionSetSelection, newsel);
					if (newsel != NULL)
						IntGroupRelease(newsel);
				} */
				MolActionCreateAndPerform(mview->mol, gMolActionDeleteAnAtom, (Int)n1);
				
			} else {
				Int nn[2];
				nn[0] = n1;
				nn[1] = n2;
				MolActionCreateAndPerform(mview->mol, gMolActionDeleteBonds, 2, nn);
			}
		}
		goto exit;
	}
	
	switch (mview->draggingMode) {
		case kMainViewMovingTrackball:
			TrackballEndDragging(mview->track, p);
			break;
		case kMainViewDraggingSelectedAtoms: {
			Vector offset;
			if (mview->isDragging) {
				calcDragOffset(mview, &offset);
			/*	if (mview->mol->is_xtal_coord)
					TransformVec(&offset, mview->mol->cell->rtr, &offset); */
				MolActionCreateAndPerform(mview->mol, gMolActionTranslateAtoms, &offset, mview->mol->selection);
		/*	} else if (clickCount == 2) {
				buf[0] = 0;
				if (MyAppCallback_getTextWithPrompt("Enter formula (e.g. CH2OCH3)", buf, sizeof buf) == 0)
					return;
				MolActionCreateAndPerform(mview->mol, SCRIPT_ACTION("s"), "dock_formula", buf);
		*/
			}
			break;
		}
		case kMainViewSelectingRegion: {
			char *selectFlags = temporarySelection(mview, mview->modifierFlags, (mview->isDragging ? 0 : 1), 1);
			if (selectFlags != NULL) {
				IntGroup *ig = IntGroupNew();
				if (ig != NULL) {
					int i, natoms;
					natoms = mview->mol->natoms;
					for (i = 0; i < natoms; i++) {
						if (selectFlags[i])
							IntGroupAdd(ig, i, 1);
					}
					MoleculeSetSelection(mview->mol, ig);
				/*	printf("current selection = {");
					for (i = j = 0; i < natoms; i++) {
						if (selectFlags[i]) {
							printf("%s%.4s", (j == 0 ? "" : " "), mview->mol->atoms[i].name);
							j++;
						}
					}
					printf("}\n"); */
					IntGroupRelease(ig);
				}
				free(selectFlags);
			}
			break;
		}
		case kMainViewCreatingBond: {
			Int b[3];
			int n1, n2;
			if (mview->tempAtoms[0] >= 0 && mview->tempAtoms[1] >= 0) {
				n1 = mview->tempAtoms[0];
				n2 = mview->tempAtoms[1];
			} else {
				if (mview->tempAtoms[0] < 0) {
					/*  Create the first atom  */
					n1 = sCreateNewAtom(mview->mol, mview->tempAtomPos[0]);
				} else n1 = mview->tempAtoms[0];
				if (mview->tempAtoms[1] < 0) {
					/*  Create the second atom, if not too close  */
					Vector dr = ATOM_AT_INDEX(mview->mol->atoms, n1)->r;
					VecDec(dr, mview->tempAtomPos[1]);
				/*	if (mview->mol->is_xtal_coord)
						TransformVec(&dr, mview->mol->cell->tr, &dr); */
					if (VecLength2(dr) > 0.01)
						n2 = sCreateNewAtom(mview->mol, mview->tempAtomPos[1]);
					else n2 = -1;
				} else n2 = mview->tempAtoms[1];
			}
			if (n1 >= 0 && n2 >= 0 && n1 != n2) {
				b[0] = n1;
				b[1] = n2;
				b[2] = kInvalidIndex;
				MolActionCreateAndPerform(mview->mol, gMolActionAddBonds, 2, b);
			/*	MoleculeAddBonds(mview->mol, b, NULL); */
			}
			break;
		}
	}
  exit:
	mview->draggingMode = 0;
	mview->isDragging = 0;
	MainViewCallback_setNeedsDisplay(mview, 1);
	MainViewCallback_setKeyboardFocus(mview);
}

void
MainView_rotateBySlider(MainView *mview, float angle, int mode, int mouseStatus, int modifierFlags)
{
	int i, n1, n2;

	if (mouseStatus == 1) {  /*  mouseDown  */
	
		if (mode == kSliderRotateBondMode) {
			if (MoleculeIsFragmentRotatable(mview->mol, mview->mol->selection, &n1, &n2, &(mview->rotateFragment))) {
				mview->rotateCenter = ATOM_AT_INDEX(mview->mol->atoms, n1)->r;
				mview->rotateAxis = ATOM_AT_INDEX(mview->mol->atoms, n2)->r;
				VecDec(mview->rotateAxis, mview->rotateCenter);
				if (VecLength2(mview->rotateAxis) < 1e-10)
					return;
				NormalizeVec(&(mview->rotateAxis), &(mview->rotateAxis));
				if (modifierFlags & kAltKeyMask) {
					/*  Option key: reverse the fragment  */
					IntGroupRelease(mview->rotateFragment);
					mview->rotateFragment = MoleculeFragmentExcludingAtoms(mview->mol, n2, 1, &n1);
					VecScaleSelf(mview->rotateAxis, -1);
				}
			} else return;
		} else if ((modifierFlags & kAltKeyMask) != 0) {
			/*  Rotate selection  */
			float tr[3];
			if (mview->mol->selection == NULL)
				return;
			mview->rotateFragment = IntGroupNewFromIntGroup(mview->mol->selection);
			TrackballGetTranslate(mview->track, tr);
			mview->rotateCenter.x = -mview->dimension * tr[0];
			mview->rotateCenter.y = -mview->dimension * tr[1];
			mview->rotateCenter.z = -mview->dimension * tr[2];
			if (mode == kSliderRotateXMode) {
				VecCross(mview->rotateAxis, mview->up, mview->lookto);
			} else {
				mview->rotateAxis = mview->up;
				VecScaleSelf(mview->rotateAxis, -1);
			}
		} else {
			/*  Rotate along x/y axis (no coordinate transform)  */
			TrackballStartDragging(mview->track, NULL, kTrackballRotateMode);
			mview->rotateFragment = NULL;  /*  This is probably not necessary  */
		}
		if (mview->rotateFragment != NULL) {
			/*  Save the original position  */
			n1 = IntGroupGetCount(mview->rotateFragment);
			mview->rotateFragmentOldPos = (Vector *)calloc(sizeof(Vector), n1);
			if (mview->rotateFragmentOldPos == NULL) {
				IntGroupRelease(mview->rotateFragment);
				mview->rotateFragment = NULL;
				return;
			}
			for (i = 0; i < n1; i++) {
				n2 = IntGroupGetNthPoint(mview->rotateFragment, i);
				mview->rotateFragmentOldPos[i] = (ATOM_AT_INDEX(mview->mol->atoms, n2))->r;
			}			
		}
	
	} else {  /*  mouseDragged or mouseUp */
		
		if (mview->rotateFragment != NULL) {
		
			/*  Restore original positions  */
			n1 = IntGroupGetCount(mview->rotateFragment);
			for (i = 0; i < n1; i++) {
				n2 = IntGroupGetNthPoint(mview->rotateFragment, i);
				(ATOM_AT_INDEX(mview->mol->atoms, n2))->r = mview->rotateFragmentOldPos[i];
			}
			/*  Rotate  */
			if (mouseStatus == 2) { /* dragged */
				MoleculeRotate(mview->mol, &(mview->rotateAxis), angle, &(mview->rotateCenter), mview->rotateFragment);
			} else {  /* mouse up */
				IntGroup *ig = IntGroupNewFromIntGroup(mview->rotateFragment);
				MolActionCreateAndPerform(mview->mol, gMolActionRotateAtoms, &(mview->rotateAxis), (double)angle, &(mview->rotateCenter), ig);
				IntGroupRelease(ig);
				IntGroupRelease(mview->rotateFragment);
				mview->rotateFragment = NULL;
				free(mview->rotateFragmentOldPos);
				mview->rotateFragmentOldPos = NULL;
			}
			
		} else {
			float quat[4];
			float cs = cos(angle / 2);
			float sn = sin(angle / 2);
			if (mouseStatus == 0) { /* mouseUp */
				TrackballEndDragging(mview->track, NULL);
			} else { /* mouseDragged */
				quat[0] = cs;
				if (mode == kSliderRotateXMode) {
					/*  Rotation along x axis  */
					quat[1] = -sn;
					quat[2] = quat[3] = 0;
				} else {
					quat[2] = sn;
					quat[1] = quat[3] = 0;
				}
				TrackballSetTemporaryRotation(mview->track, quat);
			}
		}
	}	
	
	MainViewCallback_setNeedsDisplay(mview, 1);
	MainViewCallback_setKeyboardFocus(mview);
}

#pragma mark ====== Menu Commands ======

void
MainView_selectAll(MainView *mview)
{
	IntGroup *ig;
	if (mview == NULL || mview->mol == NULL)
		return;
	ig = IntGroupNew();
	if (ig == NULL)
		return;
	IntGroupAdd(ig, 0, mview->mol->natoms);
	MoleculeSetSelection(mview->mol, ig);
	IntGroupRelease(ig);
	MainViewCallback_setNeedsDisplay(mview, 1);
}

void
MainView_selectFragment(MainView *mview)
{
	IntGroup *ig;
	if (mview == NULL || mview->mol == NULL)
		return;
	ig = MoleculeFragmentWithAtomGroups(mview->mol, MoleculeGetSelection(mview->mol), NULL);
	if (ig == NULL)
		return;
	MoleculeSetSelection(mview->mol, ig);
	IntGroupRelease(ig);
	MainViewCallback_setNeedsDisplay(mview, 1);
}

void
MainView_selectReverse(MainView *mview)
{
	IntGroup *ig;
	if (mview == NULL || mview->mol == NULL)
		return;
	ig = IntGroupNewFromIntGroup(MoleculeGetSelection(mview->mol));
	IntGroupReverse(ig, 0, mview->mol->natoms);
	MoleculeSetSelection(mview->mol, ig);
	IntGroupRelease(ig);
	MainViewCallback_setNeedsDisplay(mview, 1);
}

void
MainView_centerSelection(MainView *mview)
{
	IntGroup *ig;
	Vector c;
	float tr[3];
	if (mview == NULL || mview->mol == NULL)
		return;
	ig = MoleculeGetSelection(mview->mol);
	MoleculeCenterOfMass(mview->mol, &c, ig);
	tr[0] = -c.x / mview->dimension;
	tr[1] = -c.y / mview->dimension;
	tr[2] = -c.z / mview->dimension;
	TrackballSetTranslate(mview->track, tr);
	MainViewCallback_setNeedsDisplay(mview, 1);
}

#pragma mark ====== Pasteboard Support ======

int
MainView_copy(MainView *mview)
{
	Molecule *mol, *mol2;
	IntGroup *sel;
	Int len, result, time;
	char *p;
	if (mview == NULL || (mol = mview->mol) == NULL || (sel = MoleculeGetSelection(mol)) == NULL)
		return 1;
	if (MoleculeExtract(mol, &mol2, MoleculeGetSelection(mol), 1) != 0
	|| mol2 == NULL
	|| (p = MoleculeSerialize(mol2, &len, &time)) == NULL)
		return 2;
	result = MoleculeCallback_writeToPasteboard(kMoleculePasteboardType, p, len);
	if (result != 0)
		return result;
	if (mol2 != NULL)
		MoleculeRelease(mol2);
	mview->pasteTimeStamp = time;
	mview->pasteCount = 0;
	return 0;
}

int
MainView_delete(MainView *mview)
{
	int result;
	Molecule *mol = mview->mol;
	if (md_is_running(mview->mol->arena)) {
		MoleculeCallback_cannotModifyMoleculeDuringMDError(mview->mol);
		return -1;
	}
	result = MolActionCreateAndPerform(mol, gMolActionUnmergeMolecule, MoleculeGetSelection(mol));
	if (result == 0)
		result = MolActionCreateAndPerform(mol, gMolActionSetSelection, NULL);
	return result;
}

int
MainView_cut(MainView *mview)
{
	int result;
	if (md_is_running(mview->mol->arena)) {
		MoleculeCallback_cannotModifyMoleculeDuringMDError(mview->mol);
		return -1;
	}
	result = MainView_copy(mview);
	if (result == 0)
		result = MainView_delete(mview);
	return result;
}

int
MainView_isPastable(MainView *mview)
{
	if (MoleculeCallback_isDataInPasteboard(kMoleculePasteboardType))
		return 1;
	else return 0;
}

int
MainView_paste(MainView *mview)
{
	void *p;
	Int len, result, time;
	Molecule *mol2;
	IntGroup *sel;

	if (md_is_running(mview->mol->arena)) {
		MoleculeCallback_cannotModifyMoleculeDuringMDError(mview->mol);
		return -1;
	}	
	if (!MainView_isPastable(mview))
		return 1;
	if (MoleculeCallback_readFromPasteboard(kMoleculePasteboardType, &p, &len) != 0)
		return 2;
	if ((mol2 = MoleculeDeserialize(p, len, &time)) == NULL) {
		free(p);
		return 3;
	}
	free(p);

	if (time == mview->pasteTimeStamp) {
		/*  Offset the pasted fragment by something  */
		Vector v;
		mview->pasteCount++;
		v.x = v.y = v.z = mview->pasteCount * 0.5 + sin((double)mview->pasteCount) * 0.1; /* sin(..) is to avoid accidental coincidence  */
		MoleculeTranslate(mol2, &v, NULL);
	} else {
		mview->pasteTimeStamp = time;
		mview->pasteCount = 0;
	}
	
	sel = MoleculeGetSelection(mview->mol);
	if (sel == NULL || IntGroupGetCount(sel) == 0) {
		/*  Remove dummy atoms from mol2  */
		Atom *ap;
		int i;
		sel = IntGroupNew();
		for (i = 0, ap = mol2->atoms; i < mol2->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->atomicNumber == 0 && ap->aname[0] == '_') {
				/*  Dummy atom  */
				IntGroupAdd(sel, i, 1);
			}
		}
		if (IntGroupGetCount(sel) > 0)
			MoleculeUnmerge(mol2, NULL, sel, 0);
		IntGroupRelease(sel);
	}
	
	result = MolActionCreateAndPerform(mview->mol, SCRIPT_ACTION("M"), "dock_fragment", mol2);
	MoleculeRelease(mol2);
	return result;
}

static int
sMainView_AppendParameterToSerialBuffer(void **pb_ptr, Int *pb_len, Int parType, Int count, const UnionPar *up)
{
	int len = sizeof(Int) * 2 + sizeof(UnionPar) * count;
	char *p;
	if (*pb_len == 0) {
		*pb_ptr = malloc(len + 1);
		if (*pb_ptr == NULL)
			return -1;
		p = (char *)(*pb_ptr);
		*pb_len = len + 1;
	} else {
		*pb_ptr = realloc(*pb_ptr, *pb_len + len);
		if (*pb_ptr == NULL)
			return -1;
		p = (char *)(*pb_ptr) + *pb_len - 1;  /*  *pb_len includes the end mark  */
		*pb_len += len;
	}
	*((Int *)p) = 0;
	*p = parType;
	p += sizeof(Int);
	*((Int *)p) = count;
	p += sizeof(Int);
	memmove(p, up, sizeof(UnionPar) * count);
	p += sizeof(UnionPar) * count;
	*p = 0;  /*  End mark  */
	return 0;
}

int
MainView_pasteParameters(MainView *mview)
{
	char *p, *pp;
	Int len, i;
	IntGroup *newsel;
	if (mview == NULL || mview->mol == NULL || mview->mol->par == NULL)
		return -1;
	if (!MoleculeCallback_isDataInPasteboard(kParameterPasteboardType))
		return 1;
	if (MoleculeCallback_readFromPasteboard(kParameterPasteboardType, (void **)&p, &len) != 0)
		return 2;
	pp = p;
	newsel = IntGroupNew();
	while (*p != 0) {
		int count, c;
		UnionPar *up;
		IntGroup *ig;
		int parType = *p;
		if (parType < kFirstParType || parType > kLastParType)
			break;
		p += sizeof(Int);
		count = *((Int *)p);
		p += sizeof(Int);
		up = (UnionPar *)calloc(sizeof(UnionPar), count);
		memmove(up, p, sizeof(UnionPar) * count);

		/*  The global parameters become local when pasted  */
		for (i = 0; i < count; i++) {
			if (up[i].bond.src > 0)
				up[i].bond.src = 0;
		}
		
		c = ParameterGetCountForType(mview->mol->par, parType);
		ig = IntGroupNewWithPoints(c, count, -1);
		MolActionCreateAndPerform(mview->mol, gMolActionAddParameters, parType, ig, count, up);
		free(up);
		IntGroupRelease(ig);
		p += sizeof(UnionPar) * count;
		c = ParameterTableGetRowFromTypeAndIndex(mview->mol->par, parType, c);
		IntGroupAdd(newsel, c, count);
	}
	free(pp);
	
	/*  Select newly pasted parameters  */
	MainViewCallback_setTableSelection(mview, newsel);
	IntGroupRelease(newsel);

	/*  Request MD rebuild  */
	mview->mol->needsMDRebuild = 1;

	/*  Suppress clear of parameter table selection  */
	mview->mol->parameterTableSelectionNeedsClear = 0;

	MoleculeCallback_notifyModification(mview->mol, 0);
	return 0;
}

/*  flags: 1, delete; 2, copy; 3, cut  */
int
MainView_copyOrCutParameters(MainView *mview, int flags)
{
	IntGroup *ig = MainViewCallback_getTableSelection(mview);
	IntGroup *ig2;
	int i, n, idx, type = -1, t1;
	Parameter *par = mview->mol->par;
	void *pb_ptr = NULL;
	Int pb_len = 0;
	i = 0;
	ig2 = NULL;
	while (1) {
		n = IntGroupGetNthPoint(ig, i);
		if (n >= 0)
			idx = ParameterTableGetItemIndex(par, n, &t1);
		else {
			idx = -1;
			t1 = -2;
		}
		if (t1 != type) {
			/*  Process Parameters for the last group  */
			if (type >= kFirstParType && ig2 != NULL && (n = IntGroupGetCount(ig2)) > 0) {
				UnionPar *up = (UnionPar *)calloc(sizeof(UnionPar), n);
				if (flags & 1) {
					MolAction *act;
					if (ParameterDelete(par, type, up, ig2) < 0)
						return -1;
					act = MolActionNew(gMolActionAddParameters, type, ig2, n, up);
					MolActionCallback_registerUndo(mview->mol, act);
					MolActionRelease(act);
				} else {
					if (ParameterCopy(par, type, up, ig2) < 0)
						return -1;
				}
				if (flags & 2) {
					if (sMainView_AppendParameterToSerialBuffer(&pb_ptr, &pb_len, type, n, up) != 0)
						return -1;
				}
				free(up);
				IntGroupRelease(ig2);
				ig2 = NULL;
			}
			if (t1 == -2)
				break;
			type = t1;
		}
		if (idx >= 0) {
			if (ig2 == NULL)
				ig2 = IntGroupNew();
			IntGroupAdd(ig2, idx, 1);
		}
		i++;
	}
	if (ig2 != NULL)
		IntGroupRelease(ig2);
	IntGroupRelease(ig);
	
	if (flags & 2) {
		n = MoleculeCallback_writeToPasteboard(kParameterPasteboardType, pb_ptr, pb_len);
		if (n != 0)
			return n;
	}
	if (flags & 1) {
		/*  Clear selection  */
		MainViewCallback_setTableSelection(mview, NULL);
		/*  Request MD rebuild  */
		mview->mol->needsMDRebuild = 1;
	}
	MoleculeCallback_notifyModification(mview->mol, 0);
	return 0;
}

#pragma mark ====== Table View (also for gBuiltinParameters) ======

/*  As a special case, the table view functions will handle gBuiltinParameters when mview == NULL.  */

typedef struct ColumnInfoRecord {
	char *name;
	int width;
	int editable;
} ColumnInfoRecord;

static ColumnInfoRecord sAtomColumns[] = {
{"atom", 4, 0}, {"name", 4, 1}, {"type", 4, 1}, {"element", 4, 1}, {"residue", 6, 1},
{"x", 6, 1}, {"y", 6, 1}, {"z", 6, 1}, {"charge", 6, 1}, {NULL}
};
static ColumnInfoRecord sBondColumns[] = {
{"atoms", 9, 0}, {"names", 9, 0}, {"type", 9, 0}, {"length", 8, 0}, {"r0", 8, 0}, {"force", 8, 0}, {NULL}
};
static ColumnInfoRecord sAngleColumns[] = {
{"atoms", 12, 0}, {"names", 12, 0}, {"type", 12, 0}, {"angle", 8, 0}, {"a0", 8, 0}, {"force", 8, 0}, {NULL}
};
static ColumnInfoRecord sDihedralColumns[] = {
{"atoms", 15, 0}, {"names", 15, 0}, {"type", 15, 0}, {"dihedral", 8, 0}, {"force", 8, 0}, {"period", 4, 0}, {"phi0", 8, 0}, {NULL}
};
static ColumnInfoRecord sImproperColumns[] = {
{"atoms", 15, 0}, {"names", 15, 0}, {"type", 15, 0}, {"improper", 8, 0}, {"force", 8, 0}, {"period", 4, 0}, {"phi0", 8, 0}, {NULL}
};
static ColumnInfoRecord sParameterColumns[] = {
{"class", 5, 0}, {"type", 9, 0}, {"", 6, 0}, {"", 6, 0}, {"", 6, 0}, {"", 6, 0}, {"", 6, 0}, {"", 6, 0}, {"src", 8, 0}, {"comment", 25, 0}, {NULL}
};
static ColumnInfoRecord sMOInfoColumns[] = {
{"MO", 10, 0}, {"energy", 20, 0}, {NULL}
};
static ColumnInfoRecord *sColumnInfo[] = {
sAtomColumns, sBondColumns, sAngleColumns, sDihedralColumns, sImproperColumns, sParameterColumns, sMOInfoColumns
};
static char *sTableTitles[] = {
	"atom", "bond", "angle", "dihedral", "improper", "parameter", "MO info"
};

void
MainView_tableTitleForIndex(MainView *mview, int idx, char *buf, int bufsize)
{
	if (mview == NULL)
		idx = kMainViewParameterTableIndex;
	if (idx < 0 || idx >= sizeof(sColumnInfo) / sizeof(sColumnInfo[0])) {
		buf[0] = 0;
		return;
	}
	snprintf(buf, bufsize, "%s", sTableTitles[idx]);
}

void
MainView_createColumnsForTableAtIndex(MainView *mview, int idx)
{
	int i;

	if (mview == NULL)
		idx = kMainViewParameterTableIndex;
	if (idx < 0 || idx >= sizeof(sColumnInfo) / sizeof(sColumnInfo[0]))
		return;

	/*  Remove all existing columns  */
	while (MainViewCallback_removeTableColumnAtIndex(mview, 0) > 0);
	
	if (mview != NULL)
		mview->tableIndex = idx;
	
	/*  Create columns  */
	for (i = 0; ; i++) {
		int width;
		ColumnInfoRecord *recp = &(sColumnInfo[idx][i]);
		if (recp->name == NULL)
			break;
		width = recp->width;
		if (mview == 0 && strcmp(recp->name, "comment") == 0)
			width = 80;
		MainViewCallback_addTableColumn(mview, recp->name, width, recp->editable);
	}
}

void
MainView_refreshTable(MainView *mview)
{
	/*  Reload data  */
	if (mview != NULL && mview->mol != NULL)
		MainView_refreshCachedInfo(mview);
	
	MainViewCallback_reloadTableData(mview);

	if (mview != NULL && mview->mol != NULL && mview->tableIndex >= kMainViewAtomTableIndex && mview->tableIndex <= kMainViewImproperTableIndex)
		MainViewCallback_setTableSelection(mview, mview->tableSelection);
}

int
MainView_numberOfRowsInTable(MainView *mview)
{
	if (mview == NULL)
		return ParameterTableNumberOfRows(gBuiltinParameters);
	if (mview->mol == NULL)
		return 0;
	if (mview->tableIndex == kMainViewParameterTableIndex)
		return ParameterTableNumberOfRows(mview->mol->par);
	if (mview->tableCache == NULL)
		MainView_refreshCachedInfo(mview);
	return IntGroupGetCount(mview->tableCache);
}

static char *
sAtomDescription(Atom *ap, char *buf, int bufsize)
{
	snprintf(buf, bufsize, "%d:%.4s", ap->resSeq, ap->aname);
	return buf;
}

static UnionPar *
sParameterOfTypeAtIndex(Molecule *mol, int parType, int idx)
{
	int i;
	if (mol->arena == NULL || mol->needsMDRebuild || mol->arena->par == NULL)
		return NULL;
	switch (parType) {
		case kVdwParType:
			i = mol->arena->vdw_par_i[idx];
			if (i >= 0 && i < mol->arena->par->nvdwPars)
				return (UnionPar *)(mol->arena->par->vdwPars + i);
			else return NULL;
		case kBondParType:
			i = mol->arena->bond_par_i[idx];
			if (i >= 0 && i < mol->arena->par->nbondPars)
				return (UnionPar *)(mol->arena->par->bondPars + i);
			else return NULL;
		case kAngleParType:
			i = mol->arena->angle_par_i[idx];
			if (i >= 0 && i < mol->arena->par->nanglePars)
				return (UnionPar *)(mol->arena->par->anglePars + i);
			else return NULL;
		case kDihedralParType:
			i = mol->arena->dihedral_par_i[idx];
			if (i >= 0 && i < mol->arena->par->ndihedralPars)
				return (UnionPar *)(mol->arena->par->dihedralPars + i);
			else return NULL;
		case kImproperParType:
			i = mol->arena->improper_par_i[idx];
			if (i >= 0 && i < mol->arena->par->nimproperPars)
				return (UnionPar *)(mol->arena->par->improperPars + i);
			else return NULL;
	}
	return NULL;
}

void
MainView_valueForTable(MainView *mview, int column, int row, char *buf, int bufsize)
{
	int idx;
	Int *ip;
	Atom *ap[4];
	char descbuf[4][20], typebuf[4][8];
	Molecule *mol;

	if (mview == NULL)
		return ParameterTableGetItemText(gBuiltinParameters, column, row, buf, bufsize);

	if (mview == NULL || (mol = mview->mol) == NULL) {
		buf[0] = 0;
		return;
	}
	
	if (mview->tableIndex == kMainViewParameterTableIndex)
		return ParameterTableGetItemText(mview->mol->par, column, row, buf, bufsize);

	if (mview->tableIndex == kMainViewAtomTableIndex) { /* Atoms */
		idx = IntGroupGetNthPoint(mview->tableCache, row);
		ap[0] = ATOM_AT_INDEX(mol->atoms, idx);
		switch (column) {
			case 0: snprintf(buf, bufsize, "%d", idx); break;
			case 1: snprintf(buf, bufsize, "%.4s", ap[0]->aname); break;
			case 2: snprintf(buf, bufsize, "%.6s", AtomTypeDecodeToString(ap[0]->type, NULL)); break;
			case 3: snprintf(buf, bufsize, "%.2s", ap[0]->element); break;
			case 4: 
				if (ap[0]->resSeq == 0)
					buf[0] = 0;
				else
					snprintf(buf, bufsize, "%.4s.%d", ap[0]->resName, ap[0]->resSeq);
				break;
			case 5: snprintf(buf, bufsize, "%.3f", ap[0]->r.x); break;
			case 6: snprintf(buf, bufsize, "%.3f", ap[0]->r.y); break;
			case 7: snprintf(buf, bufsize, "%.3f", ap[0]->r.z); break;
			case 8: snprintf(buf, bufsize, "%.3f", ap[0]->charge); break;
			default: buf[0] = 0; break;
		}
	} else if (mview->tableIndex == kMainViewBondTableIndex) { /* Bonds */
		idx = IntGroupGetNthPoint(mview->tableCache, row);
		ip = mol->bonds + idx * 2;
		ap[0] = ATOM_AT_INDEX(mol->atoms, ip[0]);
		ap[1] = ATOM_AT_INDEX(mol->atoms, ip[1]);
		switch (column) {
			case 0: snprintf(buf, bufsize, "%d-%d", ip[0], ip[1]); break;
			case 1: snprintf(buf, bufsize, "%s-%s", sAtomDescription(ap[0], descbuf[0], 20), sAtomDescription(ap[1], descbuf[1], 20)); break;
			case 2: snprintf(buf, bufsize, "%.6s-%.6s", AtomTypeDecodeToString(ap[0]->type, typebuf[0]), AtomTypeDecodeToString(ap[1]->type, typebuf[1])); break;
			case 3:
				snprintf(buf, bufsize, "%.3f", MoleculeMeasureBond(mview->mol, &(ap[0]->r), &(ap[1]->r)));
				break;
			case 4:
			case 5: {
				BondPar *bp = (BondPar *)sParameterOfTypeAtIndex(mol, kBondParType, idx);
				if (bp == NULL)
					buf[0] = 0;
				else snprintf(buf, bufsize, "%.3f", (column == 4 ? bp->r0 : bp->k * INTERNAL2KCAL));
				break;
			}
			default: buf[0] = 0; break;
		}
	} else if (mview->tableIndex == kMainViewAngleTableIndex) { /* Angles */
		idx = IntGroupGetNthPoint(mview->tableCache, row);
		ip = mol->angles + idx * 3;
		ap[0] = ATOM_AT_INDEX(mol->atoms, ip[0]);
		ap[1] = ATOM_AT_INDEX(mol->atoms, ip[1]);
		ap[2] = ATOM_AT_INDEX(mol->atoms, ip[2]);
		switch (column) {
			case 0: snprintf(buf, bufsize, "%d-%d-%d", ip[0], ip[1], ip[2]); break;
			case 1: snprintf(buf, bufsize, "%s-%s-%s", sAtomDescription(ap[0], descbuf[0], 20), sAtomDescription(ap[1], descbuf[1], 20), sAtomDescription(ap[2], descbuf[2], 20)); break;
			case 2: snprintf(buf, bufsize, "%.6s-%.6s-%.6s", AtomTypeDecodeToString(ap[0]->type, typebuf[0]), AtomTypeDecodeToString(ap[1]->type, typebuf[1]), AtomTypeDecodeToString(ap[2]->type, typebuf[2])); break;
			case 3:
				snprintf(buf, bufsize, "%.3f", MoleculeMeasureAngle(mview->mol, &(ap[0]->r), &(ap[1]->r), &(ap[2]->r)));
				break;
			case 4:
			case 5: {
				AnglePar *anp = (AnglePar *)sParameterOfTypeAtIndex(mol, kAngleParType, idx);
				if (anp == NULL)
					buf[0] = 0;
				else snprintf(buf, bufsize, "%.3f", (column == 4 ? anp->a0 * kRad2Deg : anp->k * INTERNAL2KCAL));
				break;
			}
			default: buf[0] = 0; break;
		}		
	} else if (mview->tableIndex == kMainViewDihedralTableIndex || mview->tableIndex == kMainViewImproperTableIndex) { /* Dihedrals, Impropers */
		int f = (mview->tableIndex == kMainViewDihedralTableIndex);
		idx = IntGroupGetNthPoint(mview->tableCache, row);
		ip = (f ? mview->mol->dihedrals : mview->mol->impropers) + idx * 4;
		ap[0] = ATOM_AT_INDEX(mview->mol->atoms, ip[0]);
		ap[1] = ATOM_AT_INDEX(mview->mol->atoms, ip[1]);
		ap[2] = ATOM_AT_INDEX(mview->mol->atoms, ip[2]);
		ap[3] = ATOM_AT_INDEX(mview->mol->atoms, ip[3]);
		switch (column) {
			case 0: snprintf(buf, bufsize, "%d-%d-%d-%d", ip[0], ip[1], ip[2], ip[3]); break;
			case 1: snprintf(buf, bufsize, "%s-%s-%s-%s", sAtomDescription(ap[0], descbuf[0], 20), sAtomDescription(ap[1], descbuf[1], 20), sAtomDescription(ap[2], descbuf[2], 20), sAtomDescription(ap[3], descbuf[3], 20)); break;
			case 2: snprintf(buf, bufsize, "%.6s-%.6s-%.6s-%.6s", AtomTypeDecodeToString(ap[0]->type, typebuf[0]), AtomTypeDecodeToString(ap[1]->type, typebuf[1]), AtomTypeDecodeToString(ap[2]->type, typebuf[2]), AtomTypeDecodeToString(ap[3]->type, typebuf[3])); break;
			case 3:
				snprintf(buf, bufsize, "%.3f", MoleculeMeasureDihedral(mview->mol, &(ap[0]->r), &(ap[1]->r), &(ap[2]->r), &(ap[3]->r)));
				break;
			case 4:
			case 5:
			case 6: {
				TorsionPar *tp = (TorsionPar *)sParameterOfTypeAtIndex(mol, (f ? kDihedralParType : kImproperParType), idx);
				if (tp == NULL)
					buf[0] = 0;
				else if (column == 5)
					snprintf(buf, bufsize, "%d", tp->period[0]);
				else snprintf(buf, bufsize, "%.3f", (column == 4 ? tp->k[0] * INTERNAL2KCAL : tp->phi0[0] * kRad2Deg));
				break;
			}
			default: buf[0] = 0; break;
		}		
	} else if (mview->tableIndex == kMainViewMOTableIndex) { /* MO info */
		BasisSet *bset = mview->mol->bset;
		idx = row;
		buf[0] = 0;
		if (bset != NULL && idx >= 0 && idx < bset->nmos) {
			char *s = "";
			if (idx * 2 + 2 == bset->nelectrons)
				s = " (HOMO)";
			else if (idx * 2 + 1 == bset->nelectrons)
				s = " (SOMO)";
			else if ((bset->nelectrons + 1) / 2 == idx)
				s = " (LUMO)";
			switch (column) {
				case 0: snprintf(buf, bufsize, "%d%s", idx + 1, s); break;
				case 1: snprintf(buf, bufsize, "%.8f", bset->moenergies[idx]); break;
			}
		}
	}
}

/*  Set color for the locally defined or undefined MM parameters  */
int
MainView_setColorForTable(MainView *mview, int column, int row, float *fg, float *bg)
{
	int parType = -1;
	int idx;
	UnionPar *up;
	switch (mview->tableIndex) {
		case kMainViewBondTableIndex:
			parType = kBondParType;
			break;
		case kMainViewAngleTableIndex:
			parType = kAngleParType;
			break;
		case kMainViewDihedralTableIndex:
			parType = kDihedralParType;
			break;
		case kMainViewImproperTableIndex:
			parType = kImproperParType;
			break;
		default:
			return 0;
	}
	if (column < 3)
		return 0;

	idx = IntGroupGetNthPoint(mview->tableCache, row);
	up = sParameterOfTypeAtIndex(mview->mol, parType, idx);
	if (up == NULL)
		return 0;

	if (column == 3) {
		/*  Value column; warn if the value is "abnormal"  */
		int f;
		switch (parType) {
			case kBondParType: f = md_check_abnormal_bond(mview->mol->arena, mview->mol, idx); break;
			case kAngleParType: f = md_check_abnormal_angle(mview->mol->arena, mview->mol, idx); break;
			case kDihedralParType: f = md_check_abnormal_dihedral(mview->mol->arena, mview->mol, idx); break;
			case kImproperParType: f = md_check_abnormal_improper(mview->mol->arena, mview->mol, idx); break;
			default: return 0;
		}
		if (f == 1) {
			bg[0] = 1.0;
			bg[1] = bg[2] = 0.5;
			return 2;
		} else return 0;
	} else {
		if (up->bond.src == 0) {
			bg[0] = bg[1] = 1.0;
			bg[2] = 0.6;
			return 2;
		} else if (up->bond.src == -1) {
			bg[0] = 1.0;
			bg[1] = bg[2] = 0.2;
			return 2;
		}
		return 0;
	}
}

void
MainView_setSelectionFromTable(MainView *mview)
{
	IntGroup *ig, *sel;
	if (mview == NULL || mview->mol == NULL)
		return;
	if (mview->tableIndex == kMainViewMOTableIndex || mview->tableIndex == kMainViewParameterTableIndex)
		return;  /*  Do nothing  */
	ig = MainViewCallback_getTableSelection(mview);
	sel = IntGroupNew();
	if (ig != NULL) {
		int i, i1, i2;
		Int *ip;
		for (i = 0; (i1 = IntGroupGetNthPoint(ig, i)) >= 0; i++) {
			i2 = IntGroupGetNthPoint(mview->tableCache, i1);
			if (i2 < 0)
				continue;
			if (mview->tableIndex == kMainViewAtomTableIndex) {  /* Atoms */
				IntGroupAdd(sel, i2, 1);
			} else if (mview->tableIndex == kMainViewBondTableIndex) {  /* Bonds */
				ip = mview->mol->bonds + i2 * 2;
				IntGroupAdd(sel, ip[0], 1);
				IntGroupAdd(sel, ip[1], 1);
			} else if (mview->tableIndex == kMainViewAngleTableIndex) {  /* Angles */
				ip = mview->mol->angles + i2 * 3;
				IntGroupAdd(sel, ip[0], 1);
				IntGroupAdd(sel, ip[1], 1);
				IntGroupAdd(sel, ip[2], 1);
			} else if (mview->tableIndex == kMainViewDihedralTableIndex || mview->tableIndex == kMainViewImproperTableIndex) {  /* Dihedrals, impropers */
				ip = (mview->tableIndex == kMainViewDihedralTableIndex ? mview->mol->dihedrals : mview->mol->impropers) + i2 * 4;
				IntGroupAdd(sel, ip[0], 1);
				IntGroupAdd(sel, ip[1], 1);
				IntGroupAdd(sel, ip[2], 1);
				IntGroupAdd(sel, ip[3], 1);
			}
		}
		IntGroupRelease(ig);
	}
	MoleculeSetSelection(mview->mol, sel);
	IntGroupRelease(sel);
}

static void
sMainView_ParameterTableSetItemText(Molecule *mol, int column, int row, const char *s)
{
	Parameter *par;
	int type, idx;
	const char *kstr = NULL;
	UnionPar *up = NULL;
	if (mol == NULL || (par = mol->par) == NULL || row < 0)
		return;
	up = ParameterTableGetItemPtr(par, row, &type);
	if (up == NULL)
		return;
	switch (type) {
		case kVdwParType:
			switch (column) {
				case 1: kstr = "atom_type"; break;
				case 2: kstr = "eps"; break;
				case 3: kstr = "r_eq"; break;
				case 4: kstr = "eps14"; break;
				case 5: kstr = "r_eq14"; break;
				case 6: kstr = "atomic_number"; break;
				case 7: kstr = "weight"; break;
			}
			break;
		case kBondParType:
			switch (column) {
				case 1: kstr = "atom_types"; break;
				case 2: kstr = "k"; break;
				case 3: kstr = "r0"; break;
			}
			break;
		case kAngleParType:
			switch (column) {
				case 1: kstr = "atom_types"; break;
				case 2: kstr = "k"; break;
				case 3: kstr = "a0"; break;
			}
			break;
		case kDihedralParType:
		case kImproperParType:
			switch (column) {
				case 1: kstr = "atom_types"; break;
				case 2: kstr = "k"; break;
				case 3: kstr = "period"; break;
				case 4: kstr = "phi0"; break;
			}
			break;
		case kVdwPairParType:
			switch (column) {
				case 1: kstr = "atom_types"; break;
				case 2: kstr = "eps"; break;
				case 3: kstr = "r_eq"; break;
				case 4: kstr = "eps14"; break;
				case 5: kstr = "r_eq14"; break;
			}
			break;
	}
	if (column == 9)
		kstr = "comment";
	if (kstr == NULL)
		return;
	idx = ParameterTableGetItemIndex(par, row, &type);
	MolActionCreateAndPerform(mol, SCRIPT_ACTION("iissi"), "set_parameter_attr", type, idx, kstr, s, 0);
}

void
MainView_setValueForTable(MainView *mview, int column, int row, const char *buf)
{
	int idx;
	char *key;
	Molecule *mol;
	char temp[256];
	if (mview == NULL || (mol = mview->mol) == NULL)
		return;
	MainView_valueForTable(mview, column, row, temp, sizeof temp);
	if (strcmp(buf, temp) == 0 || buf[0] == 0)
		return;  /*  No change  */
	idx = IntGroupGetNthPoint(mview->tableCache, row);
	if (mview->tableIndex == kMainViewAtomTableIndex) { /* Atoms */
		switch (column) {
			case 1: key = "name"; break;
			case 2: key = "atom_type"; break;
			case 3: key = "element"; break;
			case 4: {  /*  Residue  */
				IntGroup *ig = IntGroupNewWithPoints(idx, 1, -1);
				MolActionCreateAndPerform(mol, SCRIPT_ACTION("Gs"), "assign_residue", ig, buf);
				IntGroupRelease(ig);
				return;
			}
			case 5: key = "x"; break;
			case 6: key = "y"; break;
			case 7: key = "z"; break;
			case 8: key = "charge"; break;
			default: return;
		}
		MolActionCreateAndPerform(mol, SCRIPT_ACTION("iss"), "set_atom_attr", idx, key, buf);
	} else if (mview->tableIndex == kMainViewParameterTableIndex) { /* Parameters */
		sMainView_ParameterTableSetItemText(mview->mol, column, row, buf);
	}
}

int
MainView_isTableItemEditable(MainView *mview, int column, int row)
{
	if (mview == NULL)
		return ParameterTableIsItemEditable(gBuiltinParameters, column, row);
	if (mview->mol == NULL)
		return 0;
	if (mview->tableIndex == kMainViewParameterTableIndex)
		return ParameterTableIsItemEditable(mview->mol->par, column, row);
	if (mview->tableIndex >= 0 && mview->tableIndex < sizeof(sColumnInfo) / sizeof(sColumnInfo[0]))
		return sColumnInfo[mview->tableIndex][column].editable != 0;
	else return 0;
}

int
MainView_tableType(MainView *mview)
{
	if (mview == NULL)
		return kMainViewParameterTableIndex;
	if (mview->mol == NULL)
		return -1;
	return mview->tableIndex;
}

void
MainView_dragTableSelectionToRow(MainView *mview, int row)
{
	Int *new2old, i, n, count, natoms, start_row;
	IntGroup *sel;
	if (mview == NULL || mview->mol == NULL || mview->tableIndex != 0 || row < 0 || row > (natoms = mview->mol->natoms))
		return;
	if (md_is_running(mview->mol->arena)) {
		MoleculeCallback_cannotModifyMoleculeDuringMDError(mview->mol);
		return;
	}	
	sel = MoleculeGetSelection(mview->mol);
	if (sel == NULL || (count = IntGroupGetCount(sel)) == 0)
		return;
	new2old = (Int *)calloc(sizeof(Int), natoms);
	if (new2old == NULL)
		return;

	//  Place the atoms above the target position
	for (i = n = 0; i < row; i++) {
		if (IntGroupLookupPoint(sel, i) < 0)
			new2old[n++] = i;
	}
	start_row = n;
	//  Place the atoms within the selection
	for (i = 0; i < count; i++) {
		new2old[n++] = IntGroupGetNthPoint(sel, i);
	}
	//  Place the remaining atoms
	for (i = row; i < natoms; i++) {
		if (IntGroupLookupPoint(sel, i) < 0)
			new2old[n++] = i;
	}
	MolActionCreateAndPerform(mview->mol, gMolActionRenumberAtoms, n, new2old);
	
	//  Change selection
	sel = IntGroupNewWithPoints(start_row, count, -1);
	MolActionCreateAndPerform(mview->mol, gMolActionSetSelection, sel);
	IntGroupRelease(sel);
}

IntGroup *
MainView_selectedMO(MainView *mview)
{
	if (mview == NULL || mview->mol == NULL || mview->tableIndex != kMainViewMOTableIndex)
		return NULL;
	return MainViewCallback_getTableSelection(mview);  /*  Note: the indices are 0 based  */
}

	