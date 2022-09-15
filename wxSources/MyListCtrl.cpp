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

    header = new wxWindow(this, 1001, wxPoint(0, 0), wxSize(size.x, 16));
    scroll = new wxScrolledWindow(this, 1002, wxPoint(0, 16), wxSize(size.x, (size.y <= 16 ? -1 : size.y - 16)));
    
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
    scroll->Bind(wxEVT_MOTION, &MyListCtrl::OnMotion, this);
    scroll->Bind(wxEVT_SCROLLWIN_BOTTOM, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_TOP, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_LINEDOWN, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_LINEUP, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_PAGEDOWN, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_PAGEUP, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_THUMBRELEASE, &MyListCtrl::OnScrollWin, this);
    scroll->Bind(wxEVT_SCROLLWIN_THUMBTRACK, &MyListCtrl::OnScrollWin, this);

    //  Set Fonts
    cellFont = GetFont();
    headerFont = cellFont.Smaller();
    header->SetFont(headerFont);
    {
        //  Measure line height
        wxClientDC dc(this);
        int w, h, descent, leading;
        dc.GetTextExtent(_T("M"), &w, &h, &descent, &leading, &cellFont);
        rowHeight = h;
        dc.GetTextExtent(_T("M"), &w, &h, &descent, &leading, &headerFont);
        headerHeight = h;
    }
    
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

int
MyListCtrl::GetColumnCount()
{
    return colNames.Count();
}

bool
MyListCtrl::DeleteColumn(int col)
{
    if (col < 0 || col >= GetColumnCount())
        return false;
    colNames.RemoveAt(col);
    colWidths.erase(colWidths.begin() + col);
    colFormats.erase(colFormats.begin() + col);
    SetNeedsReload();
    return true;
}

bool
MyListCtrl::InsertColumn(int col, const wxString &heading, int format, int width)
{
    if (col < 0 || col > GetColumnCount())
        return false;
    colNames.Insert(heading, col);
    colWidths.insert(colWidths.begin() + col, width);
    colFormats.insert(colFormats.begin() + col, format);
    SetNeedsReload();
    return true;
}

void
MyListCtrl::SetColumnWidth(int col, int width)
{
    if (col < 0 || col > GetColumnCount())
        return;
    colWidths[col] = width;
    SetNeedsReload();
}

void
MyListCtrl::RefreshTable(bool refreshWindow)
{
    if (dataSource == NULL) {
        nrows = ncols = 0;
    } else {
        int i;
        wxSize sz = header->GetSize();
        nrows = dataSource->GetItemCount(this);
        ncols = GetColumnCount();
        pageWidth = 0;
        for (i = 0; i < ncols; i++) {
            pageWidth += colWidths[i];
        }
        // rowHeight = dataSource->GetRowHeight(this);
        //  "+4" is for drawing marker during cell dragging
        pageHeight = rowHeight * nrows + 4;
        //  Set the subwindow infos
        sz.y = headerHeight;
        header->SetMinSize(sz);
        scroll->SetScrollbars(rowHeight, rowHeight, floor((pageWidth + rowHeight - 1) / rowHeight), nrows);
        scroll->SetVirtualSize(pageWidth, pageHeight);
        Layout();
    }
    needsReload = false;
    if (refreshWindow)
        Refresh();
}

void
MyListCtrl::SetNeedsReload(bool flag)
{
    needsReload = flag;
    if (needsReload)
        Refresh();
}

static wxBrush lightBlueBrush(wxColour(128, 128, 255));

void
MyListCtrl::OnPaint(wxPaintEvent &event)
{
    if (dataSource == NULL)
        return;

    if (needsReload)
        RefreshTable(false);

    wxColour colour;
    wxPaintDC dc(scroll);
    scroll->DoPrepareDC(dc);
    int ox, oy;
    scroll->CalcUnscrolledPosition(0, 0, &ox, &oy);
    int col = -1, row, basex;
    wxSize sz = scroll->GetClientSize();
    bool showDragTarget = (draggingRows && (dragTargetRow != mouseRow && dragTargetRow != mouseRow + 1));
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
        row = floor(oy / rowHeight);
        //  TODO: exchange the order of i and j, and get the "all-column" attribute
        //  from SetItemColor() with col = -1
        for (j = row; j < nrows; j++) {
            float fg0[4], bg0[4];
            int n0 = dataSource->SetItemColor(this, j, -1, fg0, bg0);
            y = j * rowHeight;
            if (showDragTarget && j >= dragTargetRow) {
                y += 5;
            }
            if (y > sz.y + oy)
                break;
            x = basex;
            for (i = col; i < ncols; i++) {
                float fg[4], bg[4];
                int n;
                str = dataSource->GetItemText(this, j, i);
                n = dataSource->SetItemColor(this, j, i, fg, bg);
                if (IsRowSelected(j)) {
                    if (showDragTarget) {
                        if (n & 2) {
                            bg[0] = bg[0] * 0.5 + 0.25;
                            bg[1] = bg[1] * 0.5 + 0.25;
                            bg[2] = bg[2] * 0.5 + 0.5;
                        } else if (n0 & 2) {
                            bg[0] = bg0[0] * 0.5 + 0.25;
                            bg[1] = bg0[1] * 0.5 + 0.25;
                            bg[2] = bg0[2] * 0.5 + 0.5;
                        } else {
                            bg[0] = bg[1] = 0.5;
                            bg[2] = 1.0;
                        }
                    } else {
                        if (n & 2) {
                            bg[0] = bg[0] * 0.5;
                            bg[1] = bg[1] * 0.5;
                            bg[2] = bg[2] * 0.5 + 0.5;
                        } else if (n0 & 2) {
                            bg[0] = bg0[0] * 0.5;
                            bg[1] = bg0[1] * 0.5;
                            bg[2] = bg0[2] * 0.5 + 0.5;
                        } else {
                            bg[0] = bg[1] = 0;
                            bg[2] = 1.0;
                        }
                    }
                    if (n & 1) {
                        //  Leave fg[] as they are
                    } else if (n0 & 1) {
                        fg[0] = fg0[0];
                        fg[1] = fg0[1];
                        fg[2] = fg0[2];
                    } else {
                        fg[0] = fg[1] = fg[2] = 1.0;
                    }
                } else {
                    if (n & 2) {
                        //  Leave bg[] as they are
                    } else if (n0 & 2) {
                        bg[0] = bg0[0];
                        bg[1] = bg0[1];
                        bg[2] = bg0[2];
                    } else {
                        bg[0] = bg[1] = bg[2] = 1.0;
                    }
                    if (n & 1) {
                        //  Leave fg[] as they are
                    } else if (n & 1) {
                        fg[0] = fg0[0];
                        fg[1] = fg0[1];
                        fg[2] = fg0[2];
                    } else {
                        fg[0] = fg[1] = fg[2] = 0.0;
                    }
                }
                colour.Set(bg[0] * 255, bg[1] * 255, bg[2] * 255);
                dc.SetBrush(*wxTheBrushList->FindOrCreateBrush(colour));
                dc.SetPen(*wxTRANSPARENT_PEN);
                colour.Set(fg[0] * 255, fg[1] * 255, fg[2] * 255);
                dc.SetTextForeground(colour);
                dc.DrawRectangle(x, y, colWidths[i], rowHeight - 1);
                dc.DrawText(str, x, y);
                x += colWidths[i];
                if (x > sz.y + ox)
                    break;
            }
            if (showDragTarget) {
                y = dragTargetRow * rowHeight + 1;
                dc.SetBrush(*wxRED_BRUSH);
                dc.SetPen(wxNullPen);
                dc.DrawRectangle(basex, y, x - basex, 3);
            }
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
                wxString str = colNames[i];
                dc.DrawText(str, x, 0);
            }
            x = x1;
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

int
MyListCtrl::GetItemState(int item, int stateMask)
{
    if (stateMask & wxLIST_STATE_SELECTED)
        return IsRowSelected(item) ? wxLIST_STATE_SELECTED : 0;
    else return 0;
}

bool
MyListCtrl::SetItemState(int item, int state, int stateMask)
{
    if (stateMask & wxLIST_STATE_SELECTED) {
        if (state & wxLIST_STATE_SELECTED)
            SelectRow(item);
        else
            UnselectRow(item);
    }
    return true;
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
            if (mouseMode != 3 && dataSource->IsDragAndDropEnabled(this, mouseRow)) {
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
    if (!dataSource->IsItemEditable(this, row, col))
        return;
    StartEditText(row, col);
}

void
MyListCtrl::OnLeftUp(wxMouseEvent &event)
{
    int x, y, ux, uy, row;
    bool dragged = false;
    bool selectionChanged = selectionChangeNotificationRequired;
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
    int vx, vy;
    int step = rowHeight;
    scroll->CalcUnscrolledPosition(0, 0, &vx, &vy);
    wxSize vs = scroll->GetVirtualSize();
    int orient = event.GetOrientation();
    
    if (scroll->IsAutoScrolling()) {
        //  Autoscrolling
        if (mouseMode == 3) {
            return;  //  Mouse is captured but do not autoscroll
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

//  The return rect is the client coordinate in MyListCtrl (not scroll)
bool
MyListCtrl::GetItemRectForRowAndColumn(wxRect &rect, int row, int col)
{
    int i, tx, ty, cx, cy;
    if (col < 0 || col >= ncols || row < 0 || row >= nrows)
        return false;
    cy = rowHeight * row;
    cx = 0;
    for (i = 0; i < col; i++) {
        cx += colWidths[i];
    }
    scroll->CalcScrolledPosition(cx, cy, &tx, &ty);
    rect.x = tx;
    rect.y = ty + headerHeight;
    rect.width = colWidths[col];
    rect.height = rowHeight;
    return true;
}


bool
MyListCtrl::EnsureVisible(int row, int col)
{
    wxRect r;
    if (!GetItemRectForRowAndColumn(r, row, (col == -1 ? 0 : col)))
        return false;
    r.y -= headerHeight;  //  Convert to client coord in scroll
    wxSize sz = scroll->GetClientSize();
    int scx = -1, scy = -1;
    int ux, uy;
    scroll->CalcUnscrolledPosition(r.x, r.y, &ux, &uy);
    if (col >= 0) {
        if (r.x < 0) {
            scx = floor(ux / rowHeight);
        } else if (r.x + r.width > sz.x) {
            scx = ceil((ux + r.width - sz.x) / rowHeight);
            if (scx < 0)
                scx = 0;
        }
    }  // If col is negative, then do not scroll horizontally
    if (r.y < 0) {
        scy = floor(uy / rowHeight);
    } else if (r.y + r.height > sz.y) {
        scy = ceil((uy + r.height - sz.y) / rowHeight);
        if (scy < 0)
            scy = 0;
    }
    if (scx >= 0 || scy >= 0) {
        scroll->Scroll(scx, scy);
        return true;
    } else return false;
}

void
MyListCtrl::StartEditText(int row, int col)
{
    wxRect r;
    if (!GetItemRectForRowAndColumn(r, row, col))
        return;
    r.y -= headerHeight;  //  Convert to client coord in scroll
    int i, tx, ty;
    int cy = rowHeight * row;
    int cx = 0;
    for (i = 0; i < col; i++) {
        cx += colWidths[i];
    }
    scroll->CalcScrolledPosition(cx, cy, &tx, &ty);
    if (editText == NULL) {
        editText = new wxTextCtrl(scroll, -1, "", wxPoint(r.x - 2, r.y - 2), wxSize(r.width + 4, r.height + 4), wxTE_PROCESS_ENTER);
        editText->Bind(wxEVT_CHAR, &MyListCtrl::OnCharInText, this);
    } else {
        FinalizeEdit();
        editText->SetPosition(wxPoint(r.x - 2, r.y - 2));
        editText->SetSize(wxSize(r.width + 4, r.height + 4));
    }
    EnsureVisible(row, col);
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
MyListCtrl::FinalizeEdit()
{
    if (editText != NULL) {
        wxString sval = editText->GetValue();
        if (dataSource)
            dataSource->SetItemText(this, editRow, editColumn, sval);
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
            } while (row >= 0 && !dataSource->IsItemEditable(this, row, col));
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
        } while (row >= 0 && !dataSource->IsItemEditable(this, row, col));
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
