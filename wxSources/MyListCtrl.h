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

#include "wx/scrolwin.h"

#include "wx/string.h"
#include "wx/textctrl.h"
#include "wx/dc.h"

#include <vector>

class MyListCtrl;
extern const wxEventType MyListCtrlEvent;
enum {
	MyListCtrlEvent_tableSelectionChanged,
	MyListCtrlEvent_enableTableSelectionNotification
};

enum MyListColumnFormat {
    MyLIST_FORMAT_LEFT,
    MyLIST_FORMAT_RIGHT,
    MyLIST_FORMAT_CENTRE
};

/*  Data source protocol  */
class MyListCtrlDataSource {
public:
    virtual int GetNumberOfRows(MyListCtrl *ctrl) { return 0; }
    virtual int GetNumberOfColumns(MyListCtrl *ctrl) { return 0; }
    virtual int GetColumnWidth(MyListCtrl *ctrl, int index) { return 0; }
    virtual int GetRowHeight(MyListCtrl *ctrl) { return 0; }
    virtual int GetHeaderHeight(MyListCtrl *ctrl) { return 0; }
    virtual wxString GetCellText(MyListCtrl *ctrl, int row, int col) { return _T(""); }
    virtual bool GetCellAttr(MyListCtrl *ctrl, int row, int col, wxTextAttr &attr) { return false; }
    virtual bool IsDragEnabled(MyListCtrl *ctrl, int row = -1) { return true; }
    virtual bool OnSelectionChanged(MyListCtrl *ctrl, std::vector<int> &oldsel, std::vector<int> &newsel) { return true; }
    virtual bool IsEditable(MyListCtrl *ctrl, int row = -1, int col = -1) { return true; }
    virtual bool SetCellText(MyListCtrl *ctrl, int row, int col, const wxString &str) { return true; }
    virtual bool DragSelectionToRow(MyListCtrl *ctrl, std::vector<int> &sel, int row) { return true; }

	//  If a popup menu is attached to the cell, then returns a positive integer, and *menu_titles should
	//  contain a malloc()'ed array of char* pointers (that are also malloc()'ed or strdup()'ed)
	virtual int HasPopUpMenu(MyListCtrl *ctrl, int row, int col, char ***menu_titles) { return 0; }
	virtual void OnPopUpMenuSelected(MyListCtrl *ctrl, int row, int col, int selected_index) {}

    //virtual int GetItemCount(MyListCtrl *ctrl) { return 0; }
    //virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const { return _T(""); }
    //virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value) { return 0; }
    //virtual void DragSelectionToRow(MyListCtrl *ctrl, long row) {}
    //virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column) { return false; }
    //virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl) { return false; }
    //virtual void OnSelectionChanged(MyListCtrl *ctrl) {}
    
	//  Return 1 if foreground color should be modified, 2 if background color should be modified, 3 if both
	//virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg) { return 0; }
};

class MyListCtrl: public wxWindow {

public:
	MyListCtrlDataSource *dataSource;
	wxTextCtrl *editText;
    wxWindow *header;
    wxScrolledWindow *scroll;
    
    int ncols, nrows;
    std::vector<int> colWidths;
    std::vector<int> colFormats;
    wxArrayString colNames;
    
    int headerHeight;
    int rowHeight;
    int pageWidth, pageHeight;
    std::vector<int> selection;
    
	int editRow, editColumn;
	int dragTargetRow;
    int mouseMode;
    int mouseRow;
    int lastMouseRow;
    bool draggingRows;
    std::vector<int> oldSelection;
    
	MyListCtrl();
	virtual ~MyListCtrl();
	
	bool Create(wxWindow* parent, wxWindowID wid, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize);

	void SetDataSource(MyListCtrlDataSource *source);

	void RefreshTable();

    void PrepareSelectionChangeNotification();
    bool IsRowSelected(int row);
    bool SelectRow(int row);
    bool UnselectRow(int row);
    void UnselectAllRows();
    
    void DragRows(int x, int y);
    
//	virtual wxString OnGetItemText(long item, long column) const;

//	void StartEditing(long item, long column);
//	void EndEditing(long item, long column);

//	void GetScrollPixelsPerUnit(int *xunit, int *yunit);
//	bool GetItemRectForRowAndColumn(wxRect &rect, int row, int column);
	bool FindItemAtPosition(const wxPoint &pos, int *col, int *row);
	
//	void SetItemTextForColumn(long item, long column, const wxString &text);

	void StartEditText(int col, int row);
	void EndEditTextAndRestart(bool setValueFlag, int newCol, int newRow);
	void EndEditText(bool setValueFlag = true);
    void FinalizeEdit();
    
//	void OnKeyDownOnEditText(wxKeyEvent &event);
//	void OnKillFocusOnEditText(wxFocusEvent &event);
//	void OnIdle(wxIdleEvent &event);
	
//	void OnPaintCallback(wxDC *dc);

//	bool DeleteColumn(int col);
//	bool InsertColumn(int col, const wxString &heading, int width = -1, int format = MyLIST_FORMAT_LEFT);
//  int GetNumberOfColumns();
    
    void OnPaintHeader(wxPaintEvent &event);
    void OnPaint(wxPaintEvent &event);
    void OnLeftDown(wxMouseEvent &event);
    void OnLeftUp(wxMouseEvent &event);
    void OnLeftDClick(wxMouseEvent &event);
    void OnMotion(wxMouseEvent &event);
    void OnScrollWin(wxScrollWinEvent &event);
    void OnCharInText(wxKeyEvent &event);
    
//	void OnItemSelectionChanged(wxListEvent &event);
//	void OnTableSelectionChanged(wxCommandEvent &event);
	void EnableSelectionChangeNotification(bool flag) { selectionChangeNotificationEnabled = flag; }
//	void OnEnableTableSelectionNotification(wxCommandEvent &event);

//	void OnBeginLabelEdit(wxListEvent &event);
//	void OnEndLabelEdit(wxListEvent &event);
//	void OnItemActivated(wxListEvent &event);
//	void OnBeginDrag(wxListEvent &event);
	
//	void OnChar(wxKeyEvent &event);
//	void OnMouseDown(wxMouseEvent &event);

    void OnPopUpMenuSelected(wxCommandEvent &event);
	
	void PostSelectionChangeNotification();
	
	bool selectionChangeNotificationRequired;
	bool selectionChangeNotificationEnabled;
	int lastPopUpColumn, lastPopUpRow;

private:
	DECLARE_DYNAMIC_CLASS(MyListCtrl)
	DECLARE_EVENT_TABLE()
};

#endif /* __MyListCtrl_h__ */
