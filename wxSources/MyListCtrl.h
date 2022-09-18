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
#include "wx/listctrl.h" //  some constants are included for compatibility

#include "wx/string.h"
#include "wx/textctrl.h"
#include "wx/dc.h"

#undef max
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
    virtual int GetItemCount(MyListCtrl *ctrl) { return 0; }
    virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const { return _T(""); }
    virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value) { return 0; }
    virtual void DragSelectionToRow(MyListCtrl *ctrl, long row) {}
    virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column) { return false; }
    virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl, long row = -1) { return false; }
    virtual void OnSelectionChanged(MyListCtrl *ctrl) {}
    virtual bool IsRowSelectable(MyListCtrl *ctrl, long row) { return true; }

    //  If a popup menu is attached to the cell, then returns a positive integer, and *menu_titles should
    //  contain a malloc()'ed array of char* pointers (that are also malloc()'ed or strdup()'ed)
    virtual int HasPopUpMenu(MyListCtrl *ctrl, long row, long column, char ***menu_titles) { return 0; }
    virtual void OnPopUpMenuSelected(MyListCtrl *ctrl, long row, long column, int selected_index) {}
    
    //  Return 1 if foreground color should be modified, 2 if background color should be modified, 3 if both
    virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg) { return 0; }

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
    MyListCtrlDataSource *GetDataSource() { return dataSource; }
    wxWindow *GetHeaderWindow() { return header; }
    wxScrolledWindow *GetScrolledWindow() { return scroll; }

    void SetNeedsReload(bool flag = true);
	void RefreshTable(bool refreshWindow = true);

    void PrepareSelectionChangeNotification();
    
    // wxListCtrl compatibility: only wxLIST_STATE_SELECTED is implemented
    int GetItemCount() { return (dataSource ? dataSource->GetItemCount(this) : 0); }
    int GetItemState(int item, int stateMask);
    bool SetItemState(int item, int state, int stateMask);

    bool IsRowSelected(int row);
    bool SelectRow(int row);
    bool UnselectRow(int row);
    void UnselectAllRows();
    
    void DragRows(int x, int y);
    
	bool FindItemAtPosition(const wxPoint &pos, int *col, int *row);
	
    void GetScrollPosition(int *xpos, int *ypos);
    bool SetScrollPosition(int xpos, int ypos);
    
    bool GetItemRectForRowAndColumn(wxRect &rect, int row, int column);
    bool EnsureVisible(int row, int col = -1);

	void StartEditText(int col, int row);
	void EndEditTextAndRestart(bool setValueFlag, int newCol, int newRow);
	void EndEditText(bool setValueFlag = true);
    void FinalizeEdit();
    
    int GetColumnCount();
	bool DeleteColumn(int col);
	bool InsertColumn(int col, const wxString &heading, int format = MyLIST_FORMAT_LEFT, int width = -1);
    void SetHeaderHeight(int headerHeight);
    int GetHeaderHeight() { return headerHeight; }
    void SetColumnWidth(int col, int width);

    void OnPaintHeader(wxPaintEvent &event);
    void OnPaint(wxPaintEvent &event);
    void OnLeftDown(wxMouseEvent &event);
    void OnLeftUp(wxMouseEvent &event);
    void OnLeftDClick(wxMouseEvent &event);
    void OnMotion(wxMouseEvent &event);
    void OnScrollWin(wxScrollWinEvent &event);
    void OnCharInText(wxKeyEvent &event);
    
	void EnableSelectionChangeNotification(bool flag) { selectionChangeNotificationEnabled = flag; }

    void OnPopUpMenuSelected(wxCommandEvent &event);
	
	void PostSelectionChangeNotification();
	
	bool selectionChangeNotificationRequired;
	bool selectionChangeNotificationEnabled;
    bool needsReload;
	int lastPopUpColumn, lastPopUpRow;
    wxFont cellFont;
    wxFont headerFont;

private:
	DECLARE_DYNAMIC_CLASS(MyListCtrl)
	DECLARE_EVENT_TABLE()
};

#endif /* __MyListCtrl_h__ */
