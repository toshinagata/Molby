/*
 *  MyListCtrl.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 2022/09/11.
 *  Copyright 2022 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "MyListCtrl.h"
#include "MyMBConv.h"

#include "wx/dcclient.h"
#include "wx/scrolwin.h"
#include "wx/glcanvas.h"
#include "wx/menu.h"
#include "wx/sizer.h"

const wxEventType MyListCtrlEvent = wxNewEventType();

IMPLEMENT_DYNAMIC_CLASS(MyListCtrl, wxWindow)

BEGIN_EVENT_TABLE(MyListCtrl, wxWindow)
//EVT_COMMAND(MyListCtrlEvent_tableSelectionChanged, MyListCtrlEvent, MyListCtrl::OnTableSelectionChanged)
//EVT_COMMAND(MyListCtrlEvent_enableTableSelectionNotification, MyListCtrlEvent, MyListCtrl::OnEnableTableSelectionNotification)
END_EVENT_TABLE()

MyListCtrl::MyListCtrl()
{
	editText = NULL;
    dataSource = NULL;
    header = NULL;
    scroll = NULL;
    ncols = nrows = 0;
    headerHeight = rowHeight = 0;
    pageWidth = pageHeight = 0;
    editRow = editColumn = -1;
    dragTargetRow = -1;
    mouseMode = 0;
    mouseRow = -1;
    lastMouseRow = -1;
    draggingRows = false;
    selectionChangeNotificationRequired = false;
    selectionChangeNotificationEnabled = true;
    lastPopUpColumn = lastPopUpRow = -1;

#if defined(__WXMAC__)
	//  On OSX, the default font seems to be 14-point, which is too big.
//	wxFont font = this->GetFont();
//	font.SetPointSize(12);
//	this->SetFont(font);
#endif
}

MyListCtrl::~MyListCtrl()
{
	if (editText != NULL) {
	/*	editText->Destroy(); */  /*  May be unnecessary  */
		editText = NULL;
	}
}

bool
MyListCtrl::Create(wxWindow* parent, wxWindowID wid, const wxPoint& pos, const wxSize& size)
{
	this->wxWindow::Create(parent, wid, pos, size);

    header = new wxWindow(this, 1001, wxPoint(0, 0), wxSize(size.x, 12));
    scroll = new wxScrolledWindow(this, 1002, wxPoint(0, 12), wxSize(size.x, (size.y <= 12 ? 100 : size.y) - 12));
    
    //  Set sizer
    wxBoxSizer *vsizer = new wxBoxSizer(wxVERTICAL);
    vsizer->Add(header, wxSizerFlags(0).Expand().Border(wxALL, 0));
    vsizer->Add(scroll, wxSizerFlags(1).Expand().Border(wxALL, 0));
    this->SetSizer(vsizer);
    
    //  Connect events
    header->Bind(wxEVT_PAINT, &MyListCtrl::OnPaintHeader, this);
    scroll->Bind(wxEVT_PAINT, &MyListCtrl::OnPaint, this);
    scroll->Bind(wxEVT_LEFT_DOWN, &MyListCtrl::OnLeftDown, this);
    scroll->Bind(wxEVT_LEFT_UP, &MyListCtrl::OnLeftUp, this);
    scroll->Bind(wxEVT_LEFT_DCLICK, &MyListCtrl::OnLeftDClick, this);
    scroll->Bind(wxEVT_SCROLLWIN_BOTTOM, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_TOP, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_LINEDOWN, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_LINEUP, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_PAGEDOWN, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_PAGEUP, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_THUMBRELEASE, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_THUMBTRACK, &MyListCtrl::OnScrollWin, this);

    headerHeight = 12;
    rowHeight = 12;
    
	selectionChangeNotificationRequired = false;
	selectionChangeNotificationEnabled = true;
	return true;
}

void
MyListCtrl::SetDataSource(MyListCtrlDataSource *source)
{
	dataSource = source;
	RefreshTable();
}

void
MyListCtrl::RefreshTable()
{
    if (dataSource == NULL) {
        nrows = ncols = 0;
    } else {
        int i;
        wxSize sz = header->GetSize();
        nrows = dataSource->GetNumberOfRows(this);
        ncols = dataSource->GetNumberOfColumns(this);
        colWidths.resize(ncols, 10);
        pageWidth = 0;
        for (i = 0; i < ncols; i++) {
            colWidths[i] = dataSource->GetColumnWidth(this, i);
            pageWidth += colWidths[i];
        }
        rowHeight = dataSource->GetRowHeight(this);
        headerHeight = dataSource->GetHeaderHeight(this);
        //  "+4" is for drawing marker during cell dragging
        pageHeight = rowHeight * nrows + 4;
        //  Set the subwindow infos
        sz.y = headerHeight;
        header->SetMinSize(sz);
        scroll->SetScrollbars(rowHeight, rowHeight, floor((pageWidth + rowHeight - 1) / rowHeight), nrows);
        scroll->SetVirtualSize(pageWidth, pageHeight);
        Fit();
        Layout();
    }
    Refresh();
}

static wxBrush lightBlueBrush(wxColour(128, 128, 255));

void
MyListCtrl::OnPaint(wxPaintEvent &event)
{
    if (dataSource == NULL)
        return;

    wxPaintDC dc(scroll);
    scroll->DoPrepareDC(dc);
    int ox, oy;
    scroll->CalcUnscrolledPosition(0, 0, &ox, &oy);
    int col = -1, row, basex;
    wxSize sz = scroll->GetClientSize();
    bool showDragTarget = (draggingRows && (dragTargetRow != mouseRow && dragTargetRow != mouseRow));
    //  Draw background
    int i, j;
    basex = 0;
    for (i = 0; i < ncols; i++) {
        basex += colWidths[i];
        if (basex > ox) {
            col = i;
            basex -= colWidths[i];
            break;
        }
    }
    if (col >= 0) {
        wxTextAttr attr;
        wxString str;
        int x, y;
        x = basex;
        row = floor(oy / rowHeight);
        for (i = col; i < ncols; i++) {
            for (j = row; j < nrows; j++) {
                y = j * rowHeight;
                if (showDragTarget && j >= dragTargetRow) {
                    y += 5;
                }
                if (y > sz.y + oy)
                    break;
                str = dataSource->GetCellText(this, j, i);
                if (dataSource->GetCellAttr(this, j, i, attr)) {
                    // TODO: attribute change
                }
                if (IsRowSelected(j)) {
                    if (showDragTarget) {
                        dc.SetBrush(lightBlueBrush);
                    } else {
                        dc.SetBrush(*wxBLUE_BRUSH);
                    }
                    dc.SetPen(wxNullPen);
                    dc.SetTextForeground(*wxWHITE);
                } else {
                    dc.SetBrush(*wxWHITE_BRUSH);
                    dc.SetPen(wxNullPen);
                    dc.SetTextForeground(*wxBLACK);
                }
                dc.DrawRectangle(x, y, colWidths[i], rowHeight - 1);
                dc.DrawText(str, x, y);
            }
            if (showDragTarget) {
                y = dragTargetRow * rowHeight + 1;
                dc.SetBrush(*wxRED_BRUSH);
                dc.SetPen(wxNullPen);
                dc.DrawRectangle(x, y, colWidths[i], 3);
            }
            x += colWidths[i];
            if (x > sz.y + ox)
                break;
        }
    }
}

void
MyListCtrl::OnPaintHeader(wxPaintEvent &event)
{
    wxPaintDC dc(header);
    dc.SetPen(*wxGREY_PEN);
    dc.SetBrush(*wxLIGHT_GREY_BRUSH);
    wxSize sz = header->GetSize();
    dc.DrawRectangle(0, 0, sz.x, sz.y);
    if (dataSource) {
        int ox, oy;
        int x, x1;
        int i;
        wxTextAttr attr;
        scroll->CalcUnscrolledPosition(0, 0, &ox, &oy);
        x = -ox;
        for (i = 0; i < ncols; i++) {
            x1 = x + colWidths[i];
            if (x1 > 0) {
                wxString str = dataSource->GetCellText(this, -1, i);
                if (dataSource->GetCellAttr(this, -1, i, attr)) {
                    //  TODO: do attribute change
                }
                dc.DrawText(str, x, 0);
            }
        }
    }
}

void
MyListCtrl::PrepareSelectionChangeNotification()
{
    int i;
    if (selectionChangeNotificationRequired || !selectionChangeNotificationEnabled)
        return;  //  Do nothing
    oldSelection.resize(selection.size());
    for (i = 0; i < selection.size(); i++)
        oldSelection[i] = selection[i];
    selectionChangeNotificationRequired = true;
}

bool
MyListCtrl::IsRowSelected(int row)
{
    int i;
    for (i = 0; i < selection.size(); i++) {
        if (selection[i] == row)
            return true;
    }
    return false;
}

bool
MyListCtrl::SelectRow(int row)
{
    int i;
    for (i = 0; i < selection.size(); i++) {
        if (selection[i] == row)
            return false;
        if (selection[i] > row)
            break;
    }
    selection.insert(selection.begin() + i, row);
    Refresh();
    return true;
}

bool
MyListCtrl::UnselectRow(int row)
{
    int i;
    for (i = 0; i < selection.size(); i++) {
        if (selection[i] == row) {
            selection.erase(selection.begin() + i);
            Refresh();
            return true;
        }
    }
    return false;
}

void
MyListCtrl::UnselectAllRows()
{
    selection.clear();
    Refresh();
}

//  Mouse down
//  1. On a selected row with no modifier
//    If dragging starts, then start dragging the selected cells (if dragging is enabled)
//    If mouse up on another row, then try to reorder the selected cells (if dragging is enabled)
//      (the moved cells are left selected)
//    If mouse up on the same row, then deselect all other rows
//      (and leave only this row selected)
//  2. On a selected row with shift
//    Same as 1, except that if mouse up on the same row then do nothing
//  3. On a selected row with command (Mac) or control (Win)
//    Unselect this row, and do nothing on drag or mouse-up
//  4. On an unselected row with no modifier
//    Unselect all other rows and select this row (in this handler)
//    And do the same way as 2
//  5. On an unselected row with shift
//    Select all cells between the last-clicked row and this row (in this handler)
//    And do the same way as 2
//  6. On an unselected row with command or control key
//    Select this row (in this handler), and do the same way as 2
void
MyListCtrl::OnLeftDown(wxMouseEvent &event)
{
    int x, y, ux, uy, row, modifiers, i;
    bool isRowSelected, selectionChanged = false;
    if (editText)
        EndEditText();
    x = event.GetX();
    y = event.GetY();
    scroll->CalcUnscrolledPosition(x, y, &ux, &uy);
    row = floor(uy / rowHeight);
    isRowSelected = IsRowSelected(row);
    modifiers = event.GetModifiers();
    PrepareSelectionChangeNotification();
    if ((modifiers & wxMOD_SHIFT) != 0)
        mouseMode = (isRowSelected ? 2 : 5);
#if defined(__WXMAC__) || defined(__WXOSX__)
    else if ((modifiers & wxMOD_CMD) != 0)
        mouseMode = (isRowSelected ? 3 : 6);
#else
    else if ((modifiers & wxMOD_CONTROL) != 0)
        mouseMode = (isRowSelected ? 3 : 6);
#endif
    else
        mouseMode = (isRowSelected ? 1 : 4);
    mouseRow = row;
    if (mouseMode == 3) {
        UnselectRow(row);
        selectionChanged = true;
    } else if (mouseMode == 4) {
        UnselectAllRows();
        SelectRow(row);
        selectionChanged = true;
    } else if (mouseMode == 5) {
        if (lastMouseRow >= 0) {
            int rs, re;
            if (lastMouseRow < row) {
                rs = lastMouseRow;
                re = row;
            } else {
                rs = row;
                re = lastMouseRow;
            }
            for (i = rs; i <= re; i++) {
                SelectRow(i);
            }
            selectionChanged = true;
        }
    } else if (mouseMode == 6) {
        SelectRow(row);
        selectionChanged = true;
    }
    if (!selectionChanged) {
        //  Actually no change occurred
        selectionChangeNotificationRequired = false;
    }
    Refresh();
}

void
MyListCtrl::DragRows(int x, int y)
{
    wxSize sz = scroll->GetClientSize();
    int ux, uy;
    if (y < 0)
        y = 0;
    else if (y > sz.y)
        y = sz.y;
    scroll->CalcUnscrolledPosition(x, y, &ux, &uy);
    dragTargetRow = floor(uy / rowHeight + 0.5);
    if (dragTargetRow < 0)
        dragTargetRow = 0;
    else if (dragTargetRow > nrows)
        dragTargetRow = nrows;
    Refresh();
}

void
MyListCtrl::OnMotion(wxMouseEvent &event)
{
    if (event.LeftIsDown()) {
        //  Dragging
        if (mouseMode > 0 && !scroll->HasCapture()) {
            //  Start dragging
            scroll->CaptureMouse();
            if (mouseMode != 3 && dataSource->IsDragEnabled(this, mouseRow)) {
                draggingRows = true;
            }
        }
        if (draggingRows) {
            wxPoint cp = scroll->ScreenToClient(wxGetMousePosition());
            DragRows(cp.x, cp.y);
        }
    }
}

void
MyListCtrl::OnLeftDClick(wxMouseEvent &event)
{
    //  Start editing this cell
    int x, y, ux, uy, row, col, cx, i;
    x = event.GetX();
    y = event.GetY();
    scroll->CalcUnscrolledPosition(x, y, &ux, &uy);
    row = floor(uy / rowHeight);
    cx = 0;
    for (i = 0; i < ncols; i++) {
        if ((ux >= cx && ux < cx + colWidths[i]) || i == ncols - 1) {
            col = i;
            break;
        }
        cx += colWidths[i];
    }
    if (!dataSource->IsEditable(this, col, row))
        return;
    StartEditText(col, row);
}

void
MyListCtrl::OnLeftUp(wxMouseEvent &event)
{
    int x, y, ux, uy, row, i;
    bool dragged = false;
    bool selectionChanged = selectionChangeNotifiationRequired;
    x = event.GetX();
    y = event.GetY();
    scroll->CalcUnscrolledPosition(x, y, &ux, &uy);
    row = floor(uy / rowHeight);
    PrepareSelectionChangeNotification();
    if (scroll->HasCapture()) {
        scroll->ReleaseMouse();
        dragged = true;
        if (row != mouseRow) {
            if (draggingRows) {
                //  TODO: move cells
                selectionChanged = true;
            }
        }
    }
    if (!dragged) {
        if (mouseMode == 1 || mouseMode == 4) {
            if (selection.size() != 1 || !IsRowSelected(mouseRow)) {
                UnselectAllRows();
                SelectRow(mouseRow);
                selectionChanged = true;
            }
        }
        lastMouseRow = row;
    }
    if (!selectionChanged)
        selectionChangeNotificationRequired = false;
    mouseMode = 0;
    mouseRow = -1;
    draggingRows = false;
    Refresh();
}

void
MyListCtrl::OnScrollWin(wxScrollWinEvent &event)
{
    wxPoint cp = scroll->ScreenToClient(wxGetMousePosition());
    wxSize sz = scroll->GetClientSize();
    wxEventType etype = event.GetEventType();
    wxEventType etype_org = etype;
    int vx, vy;
    int step = rowHeight;
    scroll->CalcUnscrolledPosition(0, 0, &vx, &vy);
    wxSize vs = scroll->GetVirtualSize();
    int orient = event.GetOrientation();
    
    if (scroll->IsAutoScrolling()) {
        //  Autoscrolling
        if (mouseMode == 3) {
            return;  //  Mouse is captured by do not autoscroll
        }
        if (etype == wxEVT_SCROLLWIN_LINEUP || etype == wxEVT_SCROLLWIN_LINEDOWN) {
            if (orient == wxHORIZONTAL) {
                //  Is the mouse outside the client width?
                if (cp.x < 0) {
                    etype = wxEVT_SCROLLWIN_LINEUP;
                    if (vx < step)
                        etype = 0;
                } else if (cp.x > sz.x) {
                    etype = wxEVT_SCROLLWIN_LINEDOWN;
                    if (vx > vs.x - sz.x)
                        etype = 0;
                } else etype = 0;
            } else {
                //  Is the mouse outsize the client height?
                if (cp.y < 0) {
                    etype = wxEVT_SCROLLWIN_LINEUP;
                    if (vy < step)
                        etype = 0;
                } else if (cp.y > sz.y) {
                    etype = wxEVT_SCROLLWIN_LINEDOWN;
                    if (vy > vs.y - sz.y)
                        etype = 0;
                } else etype = 0;
            }
        }
        if (etype == 0)
            return;  //  Pause scrolling
        event.SetEventType(etype);
    }
    event.Skip();
    if (draggingRows) {
        //  Handle dragging rows
        DragRows(cp.x, cp.y);
    }
    header->Refresh();  //  Adjust the header display
}

void
MyListCtrl::PostSelectionChangeNotification()
{
	if (selectionChangeNotificationRequired && selectionChangeNotificationEnabled) {
		wxCommandEvent myEvent(MyListCtrlEvent, MyListCtrlEvent_tableSelectionChanged);
		wxPostEvent(this, myEvent);
        selectionChangeNotificationRequired = false;
	}
}

//  Find item on list control
//  Pos is the client coordinate in MyListCtrl (not scroll)
bool
MyListCtrl::FindItemAtPosition(const wxPoint &pos, int *row, int *col)
{
    int r, cx, i;
    wxPoint p = this->ClientToScreen(pos);
    p = scroll->ScreenToClient(p);
    p = scroll->CalcUnscrolledPosition(p);
    r = floor(p.y / rowHeight);
    if (r < 0 || r >= nrows)
        r = -1;
    cx = 0;
    for (i = 0; i < ncols; i++) {
        if (p.x >= cx && p.x < cx + colWidths[i])
            break;
    }
    if (i >= ncols)
        i = -1;
    if (row != NULL)
        *row = r;
    if (col != NULL)
        *col = i;
    return (r >= 0 && i >= 0);
}

void
MyListCtrl::StartEditText(int row, int col)
{
    int i, tx, ty;
    int cy = rowHeight * row;
    int cx = 0;
    for (i = 0; i < col; i++) {
        cx += colWidths[i];
    }
    scroll->CalcScrolledPosition(cx, cy, &tx, &ty);
    if (editText == NULL) {
        editText = new wxTextCtrl(scroll, -1, "", wxPoint(tx - 2, ty - 2), wxSize(colWidths[col] + 4, rowHeight + 4), wxTE_PROCESS_ENTER);
        editText->Bind(wxEVT_CHAR, &MyListCtrl::OnCharInText, this);
    } else {
        FinalizeEdit();
        editText->SetPosition(wxPoint(tx - 2, ty - 2));
        editText->SetSize(wxSize(colWidths[col] + 4, rowHeight + 4));
    }
    wxSize sz = scroll->GetClientSize();
    int scx = -1, scy = -1;
    if (tx < 0) {
        scx = floor(cx / rowHeight);
    } else if (tx + colWidths[col] > sz.x) {
        scx = ceil((colWidths[col] - sz.x) / rowHeight);
        if (scx < 0)
            scx = 0;
    }
    if (ty < 0) {
        scy = floor(cy / rowHeight);
    } else if (ty + rowHeight > sz.y) {
        scy = ceil((cy + rowHeight - sz.y) / rowHeight);
        if (scy < 0)
            scy = 0;
    }
    if (scx >= 0 || scy >= 0)
        scroll->Scroll(scx, scy);
    editText->Show();
    editText->SetFocus();
    editText->SelectAll();
    editRow = row;
    editColumn = col;
    if (selection.size() != 1 || !IsRowSelected(row)) {
        UnselectAllRows();
        SelectRow(row);
    }
}

void
MyListCtrl::EndEditText(bool setValueFlag)
{
    if (!editText)
        return;
    if (setValueFlag)
        FinalizeEdit();
    editText->Hide();
    editText->Destroy();
    editText = NULL;
}
void
MyListCtrl::EndEditTextAndRestart(bool setValueFlag, int newRow, int newCol)
{
    EndEditText(setValueFlag);
    StartEditText(newRow, newCol);
}

void
MyListCtrl::OnCharInText(wxKeyEvent &event)
{
    int kc = event.GetKeyCode();
    bool shiftDown = event.ShiftDown();
    bool altDown = event.AltDown();
    int row = editRow;
    int col = editColumn;
    if (kc == WXK_RETURN) {
        if (altDown) {
            EndEditText();
        } else {
            do {
                if (shiftDown) {
                    if (row <= 0)
                        row = -1;
                    else
                        row--;
                } else {
                    if (row >= nrows - 1)
                        row = -1;
                    else
                        row++;
                }
            } while (row >= 0 && !dataSource->IsEditable(this, row, col));
            if (row >= 0) {
                StartEditText(row, col);
            } else {
                EndEditText();
            }
        }
    } else if (kc == WXK_NUMPAD_ENTER) {
        EndEditText();
    } else if (kc == WXK_TAB) {
        do {
            if (shiftDown) {
                col--;
                if (col < 0) {
                    col = ncols - 1;
                    row--;
                    if (row < 0)
                        row = -1;
                }
            } else {
                col++;
                if (col >= ncols) {
                    col = 0;
                    row++;
                    if (row >= nrows)
                        row = -1;
                }
            }
        } while (row >= 0 && !dataSource->IsEditable(this, row, col));
        if (row >= 0)
            StartEditText(row, col);
        else
            EndEditText();
    } else
        event.Skip();
}

void
MyListCtrl::OnPopUpMenuSelected(wxCommandEvent &event)
{
	if (dataSource != NULL)
		dataSource->OnPopUpMenuSelected(this, lastPopUpRow, lastPopUpColumn, event.GetId() - 1);
}
