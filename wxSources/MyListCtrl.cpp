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
#include "wx/toplevel.h"

#if wxCHECK_VERSION(3,1,0)
#define FromFrameDIP(frame, x) frame->FromDIP(x)
#define ToFrameDIP(frame, x) frame->ToDIP(x)
#else
#define FromFrameDIP(frame, x) (x)
#define ToFrameDIP(frame, x) (x)
#endif

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
//    wxFont font = this->GetFont();
//    font.SetPointSize(12);
//    this->SetFont(font);
#endif
}

MyListCtrl::~MyListCtrl()
{
    if (editText != NULL) {
    /*    editText->Destroy(); */  /*  May be unnecessary  */
        editText = NULL;
    }
}

bool
MyListCtrl::Create(wxWindow* parent, wxWindowID wid, const wxPoint& pos, const wxSize& size)
{
    this->wxWindow::Create(parent, wid, pos, size);

    header = new wxWindow(this, 1001, wxPoint(0, 0), wxSize(size.x, 16), wxWANTS_CHARS);
    scroll = new wxScrolledWindow(this, 1002, wxPoint(0, 16), wxSize(size.x, (size.y <= 16 ? -1 : size.y - 16)));
    
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
    scroll->Bind(wxEVT_CHAR, &MyListCtrl::OnCharInScroll, this);
    scroll->Bind(wxEVT_SET_FOCUS, &MyListCtrl::OnSetFocusInScroll, this);
    scroll->Bind(wxEVT_KILL_FOCUS, &MyListCtrl::OnKillFocusInScroll, this);

    //  Set Fonts
    cellFont = GetFont();
    headerFont = cellFont.Smaller();
    header->SetFont(headerFont);
    {
        //  Measure line height
        wxClientDC dc(this);
        int w, h, descent, leading;
        dc.GetTextExtent(_T("M"), &w, &h, &descent, &leading, &cellFont);
        rowHeight = ToFrameDIP(scroll, h) + 2;
        dc.GetTextExtent(_T("M"), &w, &h, &descent, &leading, &headerFont);
        headerHeight = ToFrameDIP(scroll, h) + 2;
        header->SetSize(wxSize(size.x, FromFrameDIP(scroll, headerHeight)));
        pageHeight = rowHeight;
        pageWidth = rowHeight;
        scroll->SetScrollbars(FromFrameDIP(scroll, rowHeight), FromFrameDIP(scroll, rowHeight), 1, 1, true);
    }

    //  Set sizer
    wxBoxSizer *vsizer = new wxBoxSizer(wxVERTICAL);
    vsizer->Add(header, wxSizerFlags(0).Expand().Border(wxALL, 0));
    vsizer->Add(scroll, wxSizerFlags(1).Expand().Border(wxALL, 0));
    this->SetSizer(vsizer);
    
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
    colWidths.insert(colWidths.begin() + col, ToFrameDIP(this, width));
    colFormats.insert(colFormats.begin() + col, format);
    SetNeedsReload();
    return true;
}

int
MyListCtrl::GetHeaderHeight()
{
    return FromFrameDIP(this, headerHeight);
}

void
MyListCtrl::SetColumnWidth(int col, int width)
{
    if (col < 0 || col > GetColumnCount())
        return;
    colWidths[col] = ToFrameDIP(this, width);
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
        header->SetMinSize(FromFrameDIP(this, sz));
        scroll->SetVirtualSize(FromFrameDIP(this, pageWidth), FromFrameDIP(this, pageHeight));
        int pageSize = FromFrameDIP(this, floor((pageWidth + rowHeight - 1) / rowHeight));
        if (scroll->GetScrollPageSize(wxHORIZONTAL) != pageSize)
            scroll->SetScrollPageSize(wxHORIZONTAL, pageSize);
        if (scroll->GetScrollPageSize(wxVERTICAL) != pageSize)
            scroll->SetScrollPageSize(wxVERTICAL, nrows);
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
    bool isActive;
    isActive = scroll->HasFocus();
    isPaintActive = isActive;
    wxPaintDC dc(scroll);
    scroll->DoPrepareDC(dc);
    int ox, oy;
    scroll->CalcUnscrolledPosition(0, 0, &ox, &oy);
    int col = -1, row, basex;
    wxSize sz = scroll->GetClientSize();
    wxSize szd = ToFrameDIP(this, sz);
    bool showDragTarget = (draggingRows && (dragTargetRow != mouseRow && dragTargetRow != mouseRow + 1));
    //  Draw background
    dc.SetPen(*wxTRANSPARENT_PEN);
    dc.SetBrush(*wxWHITE_BRUSH);
    dc.DrawRectangle(ox, oy, sz.x, sz.y);
    int i, j;
    int oxd, oyd;
    oxd = ToFrameDIP(this, ox);
    oyd = ToFrameDIP(this, oy);
    basex = 0;
    for (i = 0; i < ncols; i++) {
        basex += colWidths[i];
        if (basex > oxd) {
            col = i;
            basex -= colWidths[i];
            break;
        }
    }
    if (col >= 0) {
        wxTextAttr attr;
        wxString str;
        int x, y;
        int mg = 2;
        row = floor(oyd / rowHeight);
        for (j = row; j < nrows; j++) {
            float fg0[4], bg0[4];
            int n0 = dataSource->SetItemColor(this, j, -1, fg0, bg0);
            y = j * rowHeight;
            if (showDragTarget && j >= dragTargetRow) {
                y += 5;
            }
            if (y > szd.y + oyd)
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
                        float bgbase[3] = { 0, 0, 1 };
                        if (!isActive) {
                            bgbase[0] = 0.6;
                            bgbase[1] = 0.6;
                            bgbase[2] = 0.7;
                        }
                        if (n & 2) {
                            bg[0] = (bg[0] + bgbase[0]) * 0.5;
                            bg[1] = (bg[1] + bgbase[1]) * 0.5;
                            bg[2] = (bg[2] + bgbase[2]) * 0.5;
                        } else if (n0 & 2) {
                            bg[0] = (bg0[0] + bgbase[0]) * 0.5;
                            bg[1] = (bg0[1] + bgbase[1]) * 0.5;
                            bg[2] = (bg0[2] + bgbase[2]) * 0.5;
                        } else {
                            bg[0] = bgbase[0];
                            bg[1] = bgbase[1];
                            bg[2] = bgbase[2];
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
                dc.DrawRectangle(FromFrameDIP(this, x), FromFrameDIP(this, y), FromFrameDIP(this, colWidths[i]), FromFrameDIP(this, rowHeight - 1));
                dc.SetPen(*wxLIGHT_GREY_PEN);
                dc.DrawLine(FromFrameDIP(this, x), FromFrameDIP(this, y + rowHeight - 1), FromFrameDIP(this, x + colWidths[i]), FromFrameDIP(this, y + rowHeight - 1));
                if (i == ncols - 1) {
                    dc.DrawLine(FromFrameDIP(this, x + colWidths[i]), FromFrameDIP(this, y), FromFrameDIP(this, x + colWidths[i]), FromFrameDIP(this, y + rowHeight - 1));
                }
                dc.SetClippingRegion(FromFrameDIP(this, x + mg), FromFrameDIP(this, y), FromFrameDIP(this, colWidths[i] - mg * 2), FromFrameDIP(this, rowHeight - 1));
                dc.DrawText(str, FromFrameDIP(this, x + mg), FromFrameDIP(this, y));
                dc.DestroyClippingRegion();
                x += colWidths[i];
                if (x > oxd + szd.x)
                    break;
            }
            if (showDragTarget) {
                y = dragTargetRow * rowHeight + 1;
                dc.SetBrush(*wxRED_BRUSH);
                dc.SetPen(wxNullPen);
                dc.DrawRectangle(FromFrameDIP(this, basex), FromFrameDIP(this, y), FromFrameDIP(this, x - basex), FromFrameDIP(this, 3));
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
    wxSize szd = ToFrameDIP(this, sz);
    dc.DrawRectangle(0, 0, sz.x, sz.y);
    if (dataSource) {
        int ox, oy, oxd, oyd;
        int x, x1;
        int i;
        wxTextAttr attr;
        int mg = FromFrameDIP(this, 2);
        scroll->CalcUnscrolledPosition(0, 0, &ox, &oy);
        oxd = ToFrameDIP(this, ox);
        oyd = ToFrameDIP(this, oy);
        x = -oxd;
        for (i = 0; i < ncols; i++) {
            x1 = x + colWidths[i];
            if (x1 > 0) {
                wxString str = colNames[i];
                dc.DrawLine(FromFrameDIP(this, x + colWidths[i]), 0, FromFrameDIP(this, x + colWidths[i]), FromFrameDIP(this, szd.y - 1));
                dc.SetClippingRegion(FromFrameDIP(this, x + mg), 0, FromFrameDIP(this, colWidths[i] - mg * 2), FromFrameDIP(this, szd.y));
                dc.DrawText(str, FromFrameDIP(this, x + mg), 0);
                dc.DestroyClippingRegion();
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
    if (!dataSource->IsRowSelectable(this, row))
        return false;
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
    int ux, uy, row, col, modifiers, i, n;
    wxPoint pos;
    char **items;
    bool isRowSelected, selectionChanged = false;
    if (editText)
        EndEditText();
    if (!isPaintActive) {
        scroll->SetFocus();
    }
    pos = event.GetPosition();
    if (FindItemAtPosition(pos, &row, &col) && dataSource != NULL && (n = dataSource->HasPopUpMenu(this, row, col, &items)) > 0) {
        wxMenu mnu;
        for (i = 0; i < n; i++) {
            char *p = items[i];
            bool enabled = true;
            if (*p == '-') {
                if (p[1] == 0) {
                    //  Separator
                    mnu.AppendSeparator();
                    p = NULL;
                } else {
                    //  Disabled item
                    p++;
                    enabled = false;
                }
            }
            if (p != NULL) {
                wxString itemStr(p, WX_DEFAULT_CONV);
                mnu.Append(i + 1, itemStr);
                if (!enabled)
                    mnu.Enable(i + 1, false);
            }
            free(items[i]);
            items[i] = NULL;
        }
        free(items);
        lastPopUpColumn = col;
        lastPopUpRow = row;
        mnu.Bind(wxEVT_COMMAND_MENU_SELECTED, &MyListCtrl::OnPopUpMenuSelected, this);
        PopupMenu(&mnu);
        n = dataSource->GetItemCount(this);
        for (i = 0; i < n; i++)
            SetItemState(i, (i == row ? wxLIST_STATE_SELECTED : 0), wxLIST_STATE_SELECTED);
        PostSelectionChangeNotification();
        return;
    }
    
    scroll->CalcUnscrolledPosition(pos.x, pos.y, &ux, &uy);
    int uxd, uyd;
    uxd = ToFrameDIP(this, ux);
    uyd = ToFrameDIP(this, uy);
    row = floor(uyd / rowHeight);
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
    if (selectionChanged) {
        PostSelectionChangeNotification();
    } else {
        //  Actually no change occurred
        selectionChangeNotificationRequired = false;
    }
    Refresh();
}

void
MyListCtrl::DragRows(int x, int y)
{
    /*  (x, y) is physical coordinate (not DIP)  */
    wxSize sz = scroll->GetClientSize();
    int ux, uy;
    if (y < 0)
        y = 0;
    else if (y > sz.y)
        y = sz.y;
    scroll->CalcUnscrolledPosition(x, y, &ux, &uy);
    dragTargetRow = floor(ToFrameDIP(this, uy) / rowHeight + 0.5);
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
    int uxd = ToFrameDIP(this, ux);
    row = floor(ToFrameDIP(this, uy) / rowHeight);
    cx = 0;
    for (i = 0; i < ncols; i++) {
        if ((uxd >= cx && uxd < cx + colWidths[i]) || i == ncols - 1) {
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
    row = floor(ToFrameDIP(this, uy) / rowHeight);
    PrepareSelectionChangeNotification();
    if (scroll->HasCapture()) {
        scroll->ReleaseMouse();
        dragged = true;
        if (row != mouseRow) {
            if (draggingRows) {
                dataSource->DragSelectionToRow(this, dragTargetRow);
                selectionChanged = true;
            }
        }
        lastMouseRow = dragTargetRow;
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
    if (selectionChanged)
        PostSelectionChangeNotification();
    else
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
    int step = FromFrameDIP(this, rowHeight);
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
    dataSource->OnSelectionChanged(this);
    if (selectionChangeNotificationRequired && selectionChangeNotificationEnabled) {
        wxCommandEvent myEvent(MyListCtrlEvent, MyListCtrlEvent_tableSelectionChanged);
        wxPostEvent(this, myEvent);
        selectionChangeNotificationRequired = false;
    }
}

//  Find item on list control
//  Pos is the *client* coordinate in scroll (i.e. scrolled position)
bool
MyListCtrl::FindItemAtPosition(const wxPoint &pos, int *row, int *col)
{
    int r, cx, i;
    wxPoint p = scroll->CalcUnscrolledPosition(pos);
    wxPoint pd = ToFrameDIP(this, p);
    r = floor(pd.y / rowHeight);
    if (r < 0)
        r = -1;
    else if (r >= nrows)
        r = nrows;
    cx = 0;
    if (pd.x < 0)
        i = -1;
    else {
        for (i = 0; i < ncols; i++) {
            if (pd.x >= cx && pd.x < cx + colWidths[i])
                break;
            cx += colWidths[i];
        }
    }
    if (row != NULL)
        *row = r;
    if (col != NULL)
        *col = i;
    return (r >= 0 && r < nrows && i >= 0 && i < ncols);
}

//  The return rect is the *client* coordinate in scroll (i.e. scrolled position)
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
    scroll->CalcScrolledPosition(FromFrameDIP(this, cx), FromFrameDIP(this, cy), &tx, &ty);
    rect.x = tx;
    rect.y = ty;
    rect.width = FromFrameDIP(this, colWidths[col]);
    rect.height = FromFrameDIP(this, rowHeight);
    return true;
}

//  Get the left-top position in scroll unit (= rowHeight)
void
MyListCtrl::GetScrollPosition(int *xpos, int *ypos)
{
    *xpos = scroll->GetScrollPos(wxHORIZONTAL);
    *ypos = scroll->GetScrollPos(wxVERTICAL);
}

//  Scroll so that (xpos, ypos) position is left-top
//  Return false if the position is outside the scrolling limit
bool
MyListCtrl::SetScrollPosition(int xpos, int ypos)
{
    bool retval = true;
    int xlim = scroll->GetScrollLines(wxHORIZONTAL) - scroll->GetScrollPageSize(wxHORIZONTAL);
    int ylim = scroll->GetScrollLines(wxVERTICAL) - scroll->GetScrollPageSize(wxVERTICAL);
    if (xpos > xlim) {
        retval = false;
        xpos = xlim;
    }
    if (xpos < 0) {
        retval = false;
        xpos = 0;
    }
    if (ypos > ylim) {
        retval = false;
        ypos = ylim;
    }
    if (ypos < 0) {
        retval = false;
        ypos = 0;
    }
    //  TextCtrl may be moved during the scroll, so amend the position
    scroll->Scroll(xpos, ypos);
    if (editText) {
        wxRect r;
        int delta = FromFrameDIP(scroll, 2);
        if (GetItemRectForRowAndColumn(r, editRow, editColumn)) {
            editText->SetPosition(wxPoint(r.x - delta, r.y - delta));
        }
    }
    header->Refresh();
    return retval;
}

bool
MyListCtrl::EnsureVisible(int row, int col)
{
    wxRect r;
    if (!GetItemRectForRowAndColumn(r, row, (col == -1 ? 0 : col)))
        return false;
    wxSize sz = scroll->GetClientSize();
    wxSize szd = ToFrameDIP(this, sz);
    wxPoint rposd = wxPoint(ToFrameDIP(this, r.x), ToFrameDIP(this, r.y));
    wxSize rsized = wxSize(ToFrameDIP(this, r.width), ToFrameDIP(this, r.height));
    int scx = -1, scy = -1;
    int ux, uy;
    scroll->CalcUnscrolledPosition(r.x, r.y, &ux, &uy);
    if (col >= 0) {
        if (rposd.x < 0) {
            scx = floor(ToFrameDIP(this, ux) / rowHeight);
        } else if (rposd.x + rsized.x > szd.x) {
            scx = ceil((ToFrameDIP(this, ux) + rsized.x - szd.x) / rowHeight);
            if (scx < 0)
                scx = 0;
        }
    }  // If col is negative, then do not scroll horizontally
    if (rposd.y < 0) {
        scy = floor(ToFrameDIP(this, uy) / rowHeight);
    } else if (r.y + r.height > sz.y) {
        scy = ceil((ToFrameDIP(this, uy) + rsized.y - szd.y) / rowHeight);
        if (scy < 0)
            scy = 0;
    }
    if (scx >= 0 || scy >= 0) {
        scroll->Scroll(FromFrameDIP(this, scx), FromFrameDIP(this, scy));
        header->Refresh();
        return true;
    } else return false;
}

void
MyListCtrl::StartEditText(int row, int col)
{
    wxRect r;
    int delta = FromFrameDIP(this, 2);
    EnsureVisible(row, col);
    if (!GetItemRectForRowAndColumn(r, row, col))
        return;
    if (editText == NULL) {
        editText = new wxTextCtrl(scroll, -1, "", wxPoint(r.x - delta, r.y - delta), wxSize(r.width + delta * 2, r.height + delta * 2), wxTE_PROCESS_ENTER | wxTE_PROCESS_TAB | wxWANTS_CHARS);
        editText->Bind(wxEVT_CHAR, &MyListCtrl::OnCharInText, this);
        editText->Bind(wxEVT_CHAR_HOOK, &MyListCtrl::OnCharHookInText, this);
    } else {
        FinalizeEdit();
        editText->SetPosition(wxPoint(r.x - delta, r.y - delta));
        editText->SetSize(wxSize(r.width + delta * 2, r.height + delta * 2));
    }
    if (selection.size() != 1 || !IsRowSelected(row)) {
        UnselectAllRows();
        SelectRow(row);
        PostSelectionChangeNotification();
    }
    wxString str = dataSource->GetItemText(this, row, col);
    editText->SetValue(str);
    editText->Show();
    editText->SetFocus();
    editText->SelectAll();
    editRow = row;
    editColumn = col;
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
    scroll->SetFocus();
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
        if (row >= 0) {
            StartEditText(row, col);
        } else
            EndEditText();
    } else
        event.Skip();
}

void
MyListCtrl::OnCharHookInText(wxKeyEvent &event)
{
#if defined(__WXMAC__) || defined(__WXOSX__)
    //  On macOS, shift-TAB is consumed by wxWidgets even when wxTE_PROCESS_TAB and wxWANTS_CHARS
    //  are specified (why?); so, we intercept the TAB and shift-TAB here
    int kc = event.GetKeyCode();
    if (kc == WXK_TAB || kc == WXK_NUMPAD_TAB) {
        OnCharInText(event);
    } else event.Skip();
#else
    event.Skip();
#endif
}

void
MyListCtrl::OnCharInScroll(wxKeyEvent &event)
{
    int kc = event.GetKeyCode();
    if (kc == WXK_DOWN || kc == WXK_UP) {
        int row;
        if (selection.size() > 0) {
            if (selection.size() > 1) {
                if (lastMouseRow >= 0 && lastMouseRow < nrows)
                    row = lastMouseRow;
                else if (kc == WXK_DOWN)
                    row = selection[selection.size() - 1];
                else
                    row = selection[0];
            } else row = selection[0];
            if (kc == WXK_UP && row > 0)
                row--;
            else if (kc == WXK_DOWN && row < nrows - 1)
                row++;
            else return;  //  Ignore key
            UnselectAllRows();
            SelectRow(row);
            PostSelectionChangeNotification();
            Refresh();
            lastMouseRow = row;  //  Fake as if this row was clicked
        }
    } else event.Skip();
}

void
MyListCtrl::OnSetFocusInScroll(wxFocusEvent &event)
{
    Refresh();
}

void
MyListCtrl::OnKillFocusInScroll(wxFocusEvent &event)
{
    Refresh();
}

bool
MyListCtrl::HasFocusInScroll()
{
    return (scroll->HasFocus() || this->HasFocus());
}

void
MyListCtrl::OnPopUpMenuSelected(wxCommandEvent &event)
{
    if (dataSource != NULL)
        dataSource->OnPopUpMenuSelected(this, lastPopUpRow, lastPopUpColumn, event.GetId() - 1);
}
