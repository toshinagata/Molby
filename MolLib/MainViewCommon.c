/*
 *  MainViewCommon.c
 *
 *  Created by Toshi Nagata on 12/09/01.
 *  Copyright 2012 Toshi Nagata. All rights reserved.
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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#pragma mark ====== Drawing Settings ======

MainView *
MainView_new(void)
{
	MainView *mview = (MainView *)malloc(sizeof(MainView));
	if (mview != NULL) {
		/*  Initialize common members (that are used both in GUI and CMD versions  */
		memset(mview, 0, sizeof(MainView));
		mview->track = TrackballNew();
		mview->atomRadius = 0.4;
		mview->bondRadius = 0.1;
		mview->probabilityScale = 1.5382;
		mview->dimension = 10.0;
		mview->showHydrogens = mview->showDummyAtoms = mview->showExpandedAtoms = 1;
		mview->showPeriodicBox = 1;
		mview->showGraphite = 5;
	}
	return mview;
}

void
MainView_release(MainView *mview)
{
	if (mview != NULL) {
		if (mview->ref != NULL) {
			fprintf(stderr, "%s:%d:Memory leak warning: mview->ref is not NULL\n", __FILE__, __LINE__);
			return;
		}
	//	MainView_setMolecule(mview, NULL);
		TrackballRelease(mview->track);
#if !defined(__CMDMAC__)
		IntGroupRelease(mview->tableCache);
		IntGroupRelease(mview->tableSelection);
		if (mview->nlabels > 0) {
			int i;
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
#endif
	}
	free(mview);
}

/*
void
MainView_setMolecule(MainView *mview, struct Molecule *mol)
{
	if (mview == NULL || mview->mol == mol)
		return;
	if (mview->mol != NULL) {
		mview->mol->mview = NULL;  //  No need to release
		MoleculeRelease(mview->mol);
	}
	mview->mol = mol;
	if (mol != NULL) {
		MoleculeRetain(mol);
		mol->mview = mview;  //  No retain
		MainViewCallback_moleculeReplaced(mview, mol);
		MoleculeCallback_notifyModification(mol, 0);
	}
}
*/

void
MainView_setBackgroundColor(MainView *mview, float red, float green, float blue)
{
	if (mview != NULL) {
		mview->background_color[0] = red;
		mview->background_color[1] = green;
		mview->background_color[2] = blue;
		MoleculeCallback_notifyModification(mview->mol, 0);
	}
}

void
MainView_getBackgroundColor(const MainView *mview, float *rgb)
{
	if (mview != NULL) {
		rgb[0] = mview->background_color[0];
		rgb[1] = mview->background_color[1];
		rgb[2] = mview->background_color[2];
	}
}

#pragma mark ====== Graphics ======

int
MainView_insertGraphic(MainView *mview, int index, const MainViewGraphic *graphic)
{
	if (index < 0 || index >= mview->ngraphics)
		index = mview->ngraphics;
	InsertArray(&mview->graphics, &mview->ngraphics, sizeof(MainViewGraphic), index, 1, graphic);
	MoleculeCallback_notifyModification(mview->mol, 0);
	return index;
}

int
MainView_removeGraphic(MainView *mview, int index)
{
	MainViewGraphic *g;
	if (index < 0 || index >= mview->ngraphics)
		return -1;
	g = &mview->graphics[index];
	if (g->points != NULL)
		free(g->points);
	if (g->normals != NULL)
		free(g->normals);
	DeleteArray(&mview->graphics, &mview->ngraphics, sizeof(MainViewGraphic), index, 1, NULL);
	MoleculeCallback_notifyModification(mview->mol, 0);
	return index;
}

