/*
 *  MainView.h
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

#ifndef __mainview_h__
#define __mainview_h__

#include "MolLib.h"

/*  The OpenGL header location may be different for different platform  */
#if defined(__WXMAC__) || defined(__CMDMAC__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/vvector.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
	
enum {
	kShiftKeyMask = 1,
	kAltKeyMask = 2
};

struct Molecule;
struct Trackball;
struct IntGroup;

struct Label;  /*  A customized data record for drawing text in OpenGL views. Should be defined somewhere in one of the machine-dependent source files. */

typedef struct LabelRecord {
	Int labelid;      /*  LabelID is 1-based! (0 means "no label")  */
	struct Label *label;
	Int idx1, idx2;
	Vector pos;  /*  Screen position with depth  */
} LabelRecord;

enum MainViewDraggingMode {
	kMainViewMovingTrackball = 1,  /* Rotate or scale by trackball */
	kMainViewSelectingRegion = 2,  /* Selecting a region */
	kMainViewDraggingSelectedAtoms = 3, /* Dragging selection */
	kMainViewCreatingBond = 4,     /* Creating a bond */
};

enum MainViewSliderMode {
	kSliderRotateBondMode = 1,
	kSliderRotateXMode = 2,
	kSliderRotateYMode = 3
};

enum {
	kMainViewAtomTableIndex = 0,
	kMainViewBondTableIndex = 1,
	kMainViewAngleTableIndex = 2,
	kMainViewDihedralTableIndex = 3,
	kMainViewImproperTableIndex = 4,
	kMainViewParameterTableIndex = 5,
	kMainViewMOTableIndex = 6
};
	
/*  The custom python classes for MyMolView GUI  */
typedef struct MainView {
//	PyObject_HEAD
	struct Molecule *mol;
	void *ref;  /*  A platform-dependent pointer to the main view object (or the main window **controller** object)  */
	void *tableRef;  /*  The table view object  */

	unsigned char isInitialized;
	int mode;
	struct Trackball *track;
    GLdouble modelview_matrix[16];
    GLdouble projection_matrix[16];
    GLdouble perspective_vector[4];
	
	/*  Camera position and direction in object (Cartesian) coordinate  */
	Vector camera;  /*  The camera position  */
	Vector lookat;  /*  The center of the screen  */
	Vector lookto;  /*  The direction from the camera to the screen center; this is easily derived by normalizing (lookat - camera), but provided for convenience  */
	Vector up;      /*  The direction up in the screen  */

	float atomRadius; /* Scale the vdW radius by this value */
	float bondRadius; /* in angstrom */
	float probabilityScale;
	float dimension;
	
	Byte showUnitCell;
	Byte showPeriodicBox;
	Byte showExpandedAtoms;
	Byte showEllipsoids;
	Byte showHydrogens;
	Byte showDummyAtoms;
	Byte showRotationCenter;

	Byte showGraphiteFlag;
	Int  showGraphite;
	Byte showPeriodicImageFlag;
	Int  showPeriodicImage[6];  /* amin, amax, bmin, bmax, cmin, cmax  */

	Byte *visibleFlags;     /*  This is used only as internal cache;
	                            The attribute of "hidden" atom is set as (ap->exflags & kAtomHiddenFlag).  */
	Int countHidden;
	Byte freezeScreen;
	
	unsigned char lineMode;     /*  Draw the model with lines  */
	unsigned char draggingMode; /*  MainViewDraggingMode  */
	unsigned char isDragging;   /*  Becomes true if mouse moved by a certain amount  */
	int clickedAtoms[2]; /*  The object under the mouse on mouseDown event. [-1,-1]: nothing, [n,-1]: atom n, [n,m]: bond between atoms n and m  */
	int modifierFlags;  /* The modifier flags during the dragging operation */
	float dragStartPos[3];  /* If starting position is on some object, then dragStartPos[2] is the screen z-coordinates of that object. Otherwise, dragStartPos[2] = 0.5.  */
	float dragEndPos[3];    /* dragEndPos[2] is always == dragStartpos[2]  */
	Vector tempAtomPos[2];  /* The positions of the atoms forming the temporary bond */
	int tempAtoms[2];       /* The atoms forming the temporary bond */

	Int pasteCount;         /* Used to offset the pasted fragment when the same fragment is pasted multiple times */
	Int pasteTimeStamp;     /* A time stamp for the last cut/copied/pasted fragment */

	/*  Rotate fragment  */
//	int rotateBond[2];      /* The bond along which the fragment is to be rotated  */
	struct IntGroup *rotateFragment;  /*  The fragment to rotate  */
	Vector rotateCenter;
	Vector rotateAxis;
	Vector *rotateFragmentOldPos; /*  The original positions (malloc'ed pointer)  */
	
	/*  Labels  */
	Int nlabels;
	LabelRecord *labels;
	LabelRecord **sortedLabels; /* (LabelRecord *)[nlabels], internally used in drawLabels(). Should be updated when nlabels changes  */

	/*  Table view  */
	int tableIndex;  /* kMainViewAtomTableIndex, etc.  */

	/*  Caches for the table view; recalculated in MainView_refreshCachedInfo  */
	struct IntGroup *tableCache;     /* Indices of atoms etc. that are shown in the table */
	struct IntGroup *tableSelection; /* Selected rows in the table  */

} MainView;
//extern PyTypeObject MainViewType;

/*  Public functions  */
MainView *MainView_newMainView(void *ref);
void MainView_initializeOpenGLView(MainView *mview);
void MainView_release(MainView *mview);
void MainView_setMolecule(MainView *mview, struct Molecule *mol);
void MainView_refreshCachedInfo(MainView *mview);
int MainView_isAtomHidden(MainView *mview, int index);
void MainView_getCamera(MainView *mview, Vector *outCamera, Vector *outLookAt, Vector *outUp);
void MainView_resizeToFit(MainView *mview);
void MainView_drawModel(MainView *view);
void MainView_invalidateLabels(MainView *mview);
void MainView_mouseDown(MainView *view, const float *p, int eventMask);
void MainView_mouseUp(MainView *view, const float *p, int eventMask, int clickCount);
void MainView_mouseDragged(MainView *view, const float *p, int eventMask);
void MainView_mouseMoved(MainView *view, const float *p, int eventMask);
void MainView_setMode(MainView *mview, int mode);
int MainView_getMode(const MainView *mview);
void MainView_attachLabelToAtom(MainView *mview, int index);
void MainView_detachLabelFromAtom(MainView *mview, int index);
void MainView_purgeUnusedLabels(MainView *mview);
void MainView_rotateBySlider(MainView *mview, float angle, int mode, int mouseStatus, int modifierFlags);
void MainView_selectAll(MainView *mview);
void MainView_selectFragment(MainView *mview);
void MainView_selectReverse(MainView *mview);
void MainView_centerSelection(MainView *mview);
int MainView_copy(MainView *mview);
int MainView_cut(MainView *mview);
int MainView_delete(MainView *mview);
int MainView_paste(MainView *mview);
int MainView_pasteParameters(MainView *mview);
int MainView_copyOrCutParameters(MainView *mview, int flags);

/*  Table view  */
void MainView_tableTitleForIndex(MainView *mview, int idx, char *buf, int bufsize);
void MainView_createColumnsForTableAtIndex(MainView *mview, int idx);
void MainView_refreshTable(MainView *mview);
int MainView_numberOfRowsInTable(MainView *mview);
int MainView_indexToTableRow(MainView *mview, int idx);
int MainView_tableRowToIndex(MainView *mview, int row);
void MainView_valueForTable(MainView *mview, int column, int row, char *buf, int bufsize);
void MainView_setValueForTable(MainView *mview, int column, int row, const char *buf);
int MainView_setColorForTable(MainView *mview, int column, int row, float *fg, float *bg);
void MainView_setSelectionFromTable(MainView *mview);
int MainView_isTableItemEditable(MainView *mview, int column, int row);
int MainView_tableType(MainView *mview);
void MainView_dragTableSelectionToRow(MainView *mview, int row);
IntGroup *MainView_selectedMO(MainView *mview);

/*  Stubs  */
/*  These operations should work for the main *view* object  */
STUB int MainViewCallback_modifierFlags(void *eventRef);
STUB int MainViewCallback_clickCount(void *eventRef);
STUB void MainViewCallback_lockFocus(MainView *mview);
STUB void MainViewCallback_unlockFocus(MainView *mview);
STUB void MainViewCallback_frame(MainView *mview, float *rect);
STUB void MainViewCallback_display(MainView *mview);
STUB void MainViewCallback_setNeedsDisplay(MainView *mview, int flag);
STUB void MainViewCallback_setKeyboardFocus(MainView *mview);
STUB int MainViewCallback_mouseCheck(MainView *mview);
STUB void MainViewCallback_clearLabels(MainView *mview);
STUB void MainViewCallback_drawLabel(MainView *mview, const float *pos, const char *label);
STUB void MainViewCallback_drawInfoText(MainView *mview, const char *label);
STUB void MainViewCallback_selectMatrixCellForMode(MainView *mview, int mode);
//STUB int MainViewCallback_getTag(MainView *mview);
STUB MainView *MainViewCallback_viewWithTag(int tag);
STUB MainView *MainViewCallback_activeView(void);
STUB MainView *MainViewCallback_newFromFile(const char *fname);
STUB int MainViewCallback_importFromFile(MainView *mview, const char *fname);
STUB void MainViewCallback_getFilename(MainView *mview, char *buf, int bufsize);
STUB void MainViewCallback_moleculeReplaced(MainView *mview, struct Molecule *mol);

STUB struct Label *MainViewCallback_newLabel(MainView *mview, const char *message, float fontsize, const float *forecolor, const float *backcolor); /* colors are rgba */
STUB void MainViewCallback_releaseLabel(struct Label *label);
STUB void MainViewCallback_drawLabelAtPoint(struct Label *label, const float *pos);
STUB void MainViewCallback_labelSize(struct Label *label, float *outSize);
STUB int MainViewCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize);

/*  Table view  */
STUB void MainViewCallback_selectTable(MainView *mview, int idx);
STUB int MainViewCallback_numberOfTableColumns(MainView *mview);
STUB int MainViewCallback_addTableColumn(MainView *mview, const char *name, int width, int editable);
STUB int MainViewCallback_removeTableColumnAtIndex(MainView *mview, int idx);
STUB void MainViewCallback_reloadTableData(MainView *mview);
STUB void MainViewCallback_setTableSelection(MainView *mview, IntGroup *selection);
STUB IntGroup *MainViewCallback_getTableSelection(MainView *mview);
STUB void MainViewCallback_showTable(MainView *mview);
STUB void MainViewCallback_hideTable(MainView *mview);
STUB void MainViewCallback_ensureVisible(MainView *mview, int row);

/*  Register the type definition  */
//extern void MainView_register(PyObject *module);

#ifdef __cplusplus
}
#endif

#endif /* __mainview_h__ */
