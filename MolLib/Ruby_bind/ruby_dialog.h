/*
 *  RubyDialog.h
 *
 *  Created by Toshi Nagata on 08/04/13.
 *  Copyright 2008 Toshi Nagata. All rights reserved.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
*/

#ifndef __ruby_dialog_h__
#define __ruby_dialog_h__

#include "Molby_extern.h"

#ifdef __cplusplus
extern "C" {
#endif
	
/*  RubyDialog class  */
/* extern VALUE rb_cDialog; */

/*  Style of the dialog frame  */
enum {
	rd_Resizable = 1,
	rd_HasCloseBox = 2,
};

/*  True if y-coordinate grows from bottom to top (like Cocoa)  */
extern int gRubyDialogIsFlipped;
	
/*  Opaque structures (actually they are used only as pointers)  */
typedef struct RubyDialog RubyDialog;
typedef struct RDItem RDItem;

/*  Non-opaque structures  */
typedef struct RDPoint { float x, y; } RDPoint;
typedef struct RDSize  { float width, height; } RDSize;
typedef struct RDRect  { RDPoint origin; RDSize size; } RDRect;

extern const RDPoint gZeroPoint;
extern const RDSize gZeroSize;
extern const RDRect gZeroRect;
	
/*  Utility function  */
extern int RubyDialog_validateItemContent(RubyValue self, RDItem *ip, const char *s);
extern void RubyDialog_doItemAction(RubyValue self, RDItem *ip);
extern void RubyDialog_doTimerAction(RubyValue self);
extern void RubyDialog_doKeyAction(RubyValue self, int keyCode);
extern int RubyDialog_getFlexFlags(RubyValue self, RDItem *ip);
extern void RubyDialog_doCloseWindow(RubyValue self, int isModal);

extern int RubyDialog_GetTableItemCount(RubyValue self, RDItem *ip);
extern void RubyDialog_GetTableItemText(RubyValue self, RDItem *ip, int row, int column, char *buf, int blen);
extern int RubyDialog_SetTableItemText(RubyValue self, RDItem *ip, int row, int column, const char *str);
extern void RubyDialog_DragTableSelectionToRow(RubyValue self, RDItem *ip, int row);
extern int RubyDialog_IsTableItemEditable(RubyValue self, RDItem *ip, int row, int column);
extern int RubyDialog_IsTableDragAndDropEnabled(RubyValue self, RDItem *ip);
extern void RubyDialog_OnTableSelectionChanged(RubyValue self, RDItem *ip);
extern int RubyDialog_SetTableItemColor(RubyValue self, RDItem *ip, int row, int column, float *fg, float *bg);
extern int RubyDialog_HasPopUpMenu(RubyValue self, RDItem *ip, int row, int column, char ***menu_titles);
extern void RubyDialog_OnPopUpMenuSelected(RubyValue self, RDItem *ip, int row, int column, int selected_index);
	
extern void RubyDialogInitClass(void);

/*  Stub routines  */
STUB RubyDialog *RubyDialogCallback_new(int style);
STUB void RubyDialogCallback_release(RubyDialog *dref);
STUB void RubyDialogCallback_setRubyObject(RubyDialog *dref, RubyValue val);
STUB void RubyDialogCallback_setWindowTitle(RubyDialog *dref, const char *title);
STUB int RubyDialogCallback_runModal(RubyDialog *dref);
STUB void RubyDialogCallback_endModal(RubyDialog *dref, int status);
STUB void RubyDialogCallback_close(RubyDialog *dref);
STUB void RubyDialogCallback_show(RubyDialog *dref);
STUB void RubyDialogCallback_hide(RubyDialog *dref);
STUB int RubyDialogCallback_startIntervalTimer(RubyDialog *dref, float interval);
STUB void RubyDialogCallback_stopIntervalTimer(RubyDialog *dref);
STUB void RubyDialogCallback_enableOnKeyHandler(RubyDialog *dref, int flag);

STUB RDSize RubyDialogCallback_windowMinSize(RubyDialog *dref);
STUB void RubyDialogCallback_setWindowMinSize(RubyDialog *dref, RDSize size);
STUB RDSize RubyDialogCallback_windowSize(RubyDialog *dref);
STUB void RubyDialogCallback_setWindowSize(RubyDialog *dref, RDSize size);

STUB void RubyDialogCallback_setAutoResizeEnabled(RubyDialog *dref, int flag);
STUB int RubyDialogCallback_isAutoResizeEnabled(RubyDialog *dref);
STUB int RubyDialogCallback_Listen(RubyDialog *dref, void *obj, const char *objtype, const char *msg, RubyValue oval, RubyValue pval);

	
STUB void RubyDialogCallback_createStandardButtons(RubyDialog *dref, const char *oktitle, const char *canceltitle);
STUB RDItem *RubyDialogCallback_createItem(RubyDialog *dref, const char *type, const char *title, RDRect frame);
STUB RDItem *RubyDialogCallback_dialogItemAtIndex(RubyDialog *dref, int idx);
//STUB RDItem *RubyDialogCallback_itemWithTag(RubyDialog *dref, int tag);
STUB int RubyDialogCallback_indexOfItem(RubyDialog *dref, RDItem *item);
//STUB int RubyDialogCallback_tagOfItem(RDItem *item);
STUB void RubyDialogCallback_moveItemUnderView(RDItem *item, RDItem *superView, RDPoint origin);
STUB RDItem *RubyDialogCallback_superview(RDItem *item);
STUB char *RubyDialogCallback_titleOfItem(RDItem *item);
STUB RDRect RubyDialogCallback_frameOfItem(RDItem *item);
STUB void RubyDialogCallback_setFrameOfItem(RDItem *item, RDRect rect);
STUB void RubyDialogCallback_setStringToItem(RDItem *item, const char *s);
STUB void RubyDialogCallback_getStringFromItem(RDItem *item, char *buf, int bufsize);
STUB char *RubyDialogCallback_getStringPtrFromItem(RDItem *item);
STUB void RubyDialogCallback_setTitleToItem(RDItem *item, const char *s);
STUB void RubyDialogCallback_setEnabledForItem(RDItem *item, int flag);
STUB int RubyDialogCallback_isItemEnabled(RDItem *item);
STUB void RubyDialogCallback_setEditableForItem(RDItem *item, int flag);
STUB int RubyDialogCallback_isItemEditable(RDItem *item);
STUB void RubyDialogCallback_setStateForItem(RDItem *item, int state);
STUB int RubyDialogCallback_getStateForItem(RDItem *item);
STUB void RubyDialogCallback_setHiddenForItem(RDItem *item, int flag);
STUB int RubyDialogCallback_isItemHidden(RDItem *item);
STUB void RubyDialogCallback_setNeedsDisplay(RDItem *item, int flag);
STUB void RubyDialogCallback_setFontForItem(RDItem *item, int size, int family, int style, int weight);
STUB int RubyDialogCallback_getFontForItem(RDItem *item, int *size, int *family, int *style, int *weight);
STUB void RubyDialogCallback_setForegroundColorForItem(RDItem *item, const double *col);
STUB void RubyDialogCallback_setBackgroundColorForItem(RDItem *item, const double *col);
STUB void RubyDialogCallback_getForegroundColorForItem(RDItem *item, double *col);
STUB void RubyDialogCallback_getBackgroundColorForItem(RDItem *item, double *col);
STUB int RubyDialogCallback_appendString(RDItem *item, const char *str);

STUB int RubyDialogCallback_countSubItems(RDItem *item);
STUB int RubyDialogCallback_appendSubItem(RDItem *item, const char *s);
STUB int RubyDialogCallback_insertSubItem(RDItem *item, const char *s, int pos);
STUB int RubyDialogCallback_deleteSubItem(RDItem *item, int pos);
STUB char *RubyDialogCallback_titleOfSubItem(RDItem *item, int pos);
STUB void RubyDialogCallback_setSelectedSubItem(RDItem *item, int pos);	
STUB int RubyDialogCallback_selectedSubItem(RDItem *item);
	
STUB RDSize RubyDialogCallback_sizeOfString(RDItem *item, const char *s);
STUB RDSize RubyDialogCallback_resizeToBest(RDItem *item);

STUB char RubyDialogCallback_insertTableColumn(RDItem *item, int col, const char *heading, int format, int width);
STUB char RubyDialogCallback_deleteTableColumn(RDItem *item, int col);
STUB int RubyDialogCallback_countTableColumn(RDItem *item);
STUB char RubyDialogCallback_isTableRowSelected(RDItem *item, int row);
STUB char RubyDialogCallback_setTableRowSelected(RDItem *item, int row, int flag);
STUB void RubyDialogCallback_refreshTable(RDItem *item);

STUB int RubyDialogCallback_savePanel(const char *title, const char *dirname, const char *wildcard, char *buf, int bufsize);
STUB int RubyDialogCallback_openPanel(const char *title, const char *dirname, const char *wildcard, char ***array, int for_directories, int multiple_selection);

#ifdef __cplusplus
}
#endif
		
#endif /* __ruby_dialog_h__ */
