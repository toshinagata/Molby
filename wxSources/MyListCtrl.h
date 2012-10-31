/*
 *  MyListCtrl.h
 *
 *  Created by Toshi Nagata on 08/12/09.
 *  Copyright 2008-2009 Toshi Nagata. All rights reserved.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 */

#ifndef __MyListCtrl_h__
#define __MyListCtrl_h__

#include "wx/listctrl.h"
#include "wx/string.h"
#include "wx/textctrl.h"
#include "wx/dc.h"

#if __WXMSW__
#include <wx/generic/listctrl.h>
#endif

class MyListCtrl;
extern const wxEventType MyListCtrlEvent;
enum {
	MyListCtrlEvent_tableSelectionChanged,
	MyListCtrlEvent_enableTableSelectionNotification
};

/*  Data source protocol  */
class MyListCtrlDataSource {
public:
	virtual int GetItemCount(MyListCtrl *ctrl) { return 0; }
	virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const { return _T(""); }
	virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value) { return 0; }
	virtual void DragSelectionToRow(MyListCtrl *ctrl, long row) {}
	virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column) { return false; }
	virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl) { return false; }
	virtual void OnSelectionChanged(MyListCtrl *ctrl) {}
	
	//  Return 1 if foreground color should be modified, 2 if background color should be modified, 3 if both
	virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg) { return 0; }
};

class MyListCtrl: public wxGenericListCtrl {

public:
	MyListCtrlDataSource *dataSource;
	wxTextCtrl *editText;
	int editRow, editColumn;
	wxListItemAttr *subTitleRowAttr;  /*  Used in SetItemColor()  */
	int dragTargetRow;

	MyListCtrl();
	virtual ~MyListCtrl();
	
	bool Create(wxWindow* parent, wxWindowID wid, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize);

	void SetDataSource(MyListCtrlDataSource *source);

	void RefreshTable();

	virtual wxString OnGetItemText(long item, long column) const;
	virtual wxListItemAttr *OnGetItemAttr(long item) const;

	void StartEditing(long item, long column);
	void EndEditing(long item, long column);

	void GetScrollPixelsPerUnit(int *xunit, int *yunit);
	bool GetItemRectForRowAndColumn(wxRect &rect, int row, int column);
	bool FindItemAtPosition(const wxPoint &pos, int *row, int *column);
	
	void SetItemTextForColumn(long item, long column, const wxString &text);

	void StartEditText(int row, int column);
	void EndEditTextAndRestart(bool setValueFlag, int newRow, int newColumn);
	void EndEditText(bool setValueFlag = true);
	void OnKeyDownOnEditText(wxKeyEvent &event);
	void OnKillFocusOnEditText(wxFocusEvent &event);
	void OnIdle(wxIdleEvent &event);
	
	void OnPaintCallback(wxDC *dc);

	/*  Override the wxListCtrl functions to take care of the internal text editor  */
	bool DeleteColumn(int col);
	bool InsertColumn(long col, const wxString &heading, int format = wxLIST_FORMAT_LEFT, int width = -1);
	
	void OnItemSelectionChanged(wxListEvent &event);
	void OnTableSelectionChanged(wxCommandEvent &event);
	void EnableSelectionChangeNotification(bool flag) { selectionChangeNotificationEnabled = flag; }
	void OnEnableTableSelectionNotification(wxCommandEvent &event);

	void OnBeginLabelEdit(wxListEvent &event);
	void OnEndLabelEdit(wxListEvent &event);
	void OnItemActivated(wxListEvent &event);
	void OnBeginDrag(wxListEvent &event);
	
	void OnChar(wxKeyEvent &event);
	void OnMouseDown(wxMouseEvent &event);
	void OnLeftDClick(wxMouseEvent &event);
	
	void PostSelectionChangeNotification();
	
	bool selectionChangeNotificationSent;
	bool selectionChangeNotificationEnabled;
	
private:
	DECLARE_DYNAMIC_CLASS(MyListCtrl)
	DECLARE_EVENT_TABLE()
};

#endif /* __MyListCtrl_h__ */
