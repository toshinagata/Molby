/*
 *  MyListCtrl.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 08/12/09.
 *  Copyright 2008 Toshi Nagata. All rights reserved.
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

const wxEventType MyListCtrlEvent = wxNewEventType();

IMPLEMENT_DYNAMIC_CLASS(MyListCtrl, wxGenericListCtrl)

BEGIN_EVENT_TABLE(MyListCtrl, wxGenericListCtrl)
EVT_LIST_ITEM_SELECTED(-1, MyListCtrl::OnItemSelectionChanged)
EVT_LIST_ITEM_DESELECTED(-1, MyListCtrl::OnItemSelectionChanged)
EVT_COMMAND(MyListCtrlEvent_tableSelectionChanged, MyListCtrlEvent, MyListCtrl::OnTableSelectionChanged)
EVT_COMMAND(MyListCtrlEvent_enableTableSelectionNotification, MyListCtrlEvent, MyListCtrl::OnEnableTableSelectionNotification)
EVT_LIST_BEGIN_DRAG(-1, MyListCtrl::OnBeginDrag)
EVT_LEFT_DCLICK(MyListCtrl::OnLeftDClick)
EVT_CHAR(MyListCtrl::OnChar)
EVT_LEFT_DOWN(MyListCtrl::OnMouseDown)
END_EVENT_TABLE()

MyListCtrl::MyListCtrl()
{
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
	this->wxGenericListCtrl::Create(parent, wid, pos, size, wxLC_REPORT | wxLC_VIRTUAL | wxBORDER_SIMPLE);
	dataSource = NULL;
	editText = NULL;
	selectionChangeNotificationSent = false;
	selectionChangeNotificationEnabled = true;
	subTitleRowAttr = new wxListItemAttr;
	dragTargetRow = -1;
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
	if (dataSource != NULL) {
		int nrows = dataSource->GetItemCount(this);
		SetItemCount(nrows);
		if (nrows > 0) {
			RefreshItems(0, nrows - 1);
		}
	}
}

// Define the repainting behaviour
void
MyListCtrl::OnPaintCallback(wxDC *dc)
{
	if (dragTargetRow >= 0) {
		wxRect r;
		wxPen pen = *wxCYAN_PEN;
		int dx, dy, y;
		/*  m_mainWin is a protected member in wxGenericListCtrl  */
		((wxScrolledWindow *)m_mainWin)->CalcScrolledPosition(0, 0, &dx, &dy);	
		pen.SetWidth(3);
		dc->SetPen(pen);
		if (dragTargetRow == GetItemCount()) {
			GetItemRect(dragTargetRow - 1, r);
			y = r.y - dy;
		} else {
			GetItemRect(dragTargetRow, r);
			y = r.y - dy - r.height;
		}
	/*	printf("dragTargetRow = %d, r.y = %d, y = %d\n", dragTargetRow, r.y, y); */
		dc->DrawLine(r.x - dx, y, r.x - dx + r.width, y);
		
	/*	r.Inflate(-1);
		dc->DrawLine(r.x, r.y, r.x + r.width, r.y);
		dc->DrawLine(r.x + r.width, r.y, r.x + r.width, r.y - r.height);
		dc->DrawLine(r.x + r.width, r.y - r.height, r.x, r.y - r.height);
		dc->DrawLine(r.x, r.y - r.height, r.x, r.y); */
	}
}

//  Callback function that is called from wxListMainWindow::OnPaint().
//  This is a very ugly hack, but I can think of no alternative...
void
wxListCtrl_onPaintCallback(wxGenericListCtrl *listctrl, wxDC *dc)
{
	MyListCtrl *ctrl = wxStaticCast(listctrl, MyListCtrl);
	if (ctrl != NULL && dc != NULL)
		ctrl->OnPaintCallback(dc);
}

wxString
MyListCtrl::OnGetItemText(long item, long column) const
{
	if (dataSource == NULL)
		return wxEmptyString;
	return dataSource->GetItemText((MyListCtrl *)this, item, column);
}

wxListItemAttr *
MyListCtrl::OnGetItemAttr(long item) const
{
	float fg[3], bg[3];
	long row, col;
	int ret;
	row = item % 1000000;
	col = item / 1000000 - 1;  //  -1 for all rows
	ret = dataSource->SetItemColor((MyListCtrl *)this, row, col, fg, bg);
	if (ret == 0)
		return NULL;
	if (ret & 1) {
		wxColour fgcol(fg[0] * 255, fg[1] * 255, fg[2] * 255);
		subTitleRowAttr->SetTextColour(fgcol);
	} else subTitleRowAttr->SetTextColour(*wxBLACK);
	if (ret & 2) {
		wxColour bgcol(bg[0] * 255, bg[1] * 255, bg[2] * 255);
		subTitleRowAttr->SetBackgroundColour(bgcol);
	} else subTitleRowAttr->SetBackgroundColour(*wxWHITE);
	return subTitleRowAttr;
}

void
MyListCtrl::PostSelectionChangeNotification()
{
	if (!selectionChangeNotificationSent && selectionChangeNotificationEnabled) {
		selectionChangeNotificationSent = true;
		wxCommandEvent myEvent(MyListCtrlEvent, MyListCtrlEvent_tableSelectionChanged);
		wxPostEvent(this, myEvent);
	}
}

void
MyListCtrl::OnBeginLabelEdit(wxListEvent &event)
{
//	printf("OnBeginLabelEdit: item index = %d\n", event.GetIndex());
}

void
MyListCtrl::OnEndLabelEdit(wxListEvent &event)
{
//	printf("OnEndLabelEdit: item index = %d\n", event.GetIndex());
}

void
MyListCtrl::OnItemActivated(wxListEvent &event)
{
//	printf("OnItemActivated: item index = %d, col = %d\n", event.GetIndex(), event.GetItem().m_col);
}

void
MyListCtrl::OnBeginDrag(wxListEvent &event)
{
	int count = GetItemCount();
	
	if (dataSource == NULL || !dataSource->IsDragAndDropEnabled(this))
		return;
	EndEditText(true);
	dragTargetRow = -1;
	while (1) {
		wxRect r;
		wxMouseState mstate = wxGetMouseState();
		wxPoint pt(mstate.GetX(), mstate.GetY());
		pt = ScreenToClient(pt);
		long newRow = FindItem(-1, pt, 0);
		if (newRow != dragTargetRow) {
			if (newRow >= 0)
				EnsureVisible(newRow);
			else {
				GetItemRect(0, r);
				if (pt.y < r.y)
					EnsureVisible(0);
				else {
					GetItemRect(count - 1, r);
					if (pt.y > r.y) {
						EnsureVisible(count - 1);
						if (pt.y < r.y + r.height)
							newRow = count;
					}
				}
			}
			if (newRow >= 0) {
				if (newRow == count) {
					GetItemRect(newRow - 1, r);
					r.y += r.height / 2;
				} else {
					GetItemRect(newRow, r);
					r.y -= r.height / 2;
				}
				RefreshRect(r);
			}
			if (dragTargetRow >= 0) {
				if (dragTargetRow == count) {
					GetItemRect(dragTargetRow - 1, r);
					r.y += r.height / 2;
				} else {
					GetItemRect(dragTargetRow, r);
					r.y -= r.height / 2;
				}
				RefreshRect(r);
			}
			dragTargetRow = newRow;
			Update();
		}
		if (!mstate.LeftDown()) {
			//  If the mouse cursor is outside the item rect, then dragging should be discarded
			if (dragTargetRow >= 0) {
				r = GetClientRect();
				if (!r.Contains(pt))
					dragTargetRow = -1;
			}
			break;
		}
	}
	if (dragTargetRow >= 0)
		dataSource->DragSelectionToRow(this, dragTargetRow);
	dragTargetRow = -1;
	Update();
}

bool
MyListCtrl::GetItemRectForRowAndColumn(wxRect &rect, int row, int column)
{
	int i, xpos, width, xunit, yunit;
	if (!GetItemRect(row, rect))
		return false;
	GetScrollPixelsPerUnit(&xunit, &yunit);
	xpos = -GetScrollPos(wxHORIZONTAL) * xunit;
	for (i = 0; i < column; i++) {
		width = GetColumnWidth(i);
		xpos += width;
	}
	rect.SetX(xpos);
	rect.SetWidth(GetColumnWidth(column));
	return true;
}

void
MyListCtrl::GetScrollPixelsPerUnit(int *xunit, int *yunit)
{
	int x, y;
	/*  m_mainWin is a protected member in wxGenericListCtrl  */
	((wxScrolledWindow *)m_mainWin)->GetScrollPixelsPerUnit(&x, &y);	
	if (xunit != NULL)
		*xunit = x;
	if (yunit != NULL)
		*yunit = y;
}

bool
MyListCtrl::FindItemAtPosition(const wxPoint &pos, int *row, int *column)
{
	int i, r, ncols, width, xpos, flags, xunit, yunit;
	r = (int)HitTest(pos, flags, NULL);
	if (r == wxNOT_FOUND)
		return false;
	ncols = GetColumnCount();
	GetScrollPixelsPerUnit(&xunit, &yunit);
	xpos = -GetScrollPos(wxHORIZONTAL) * xunit;
	for (i = 0; i < ncols; i++) {
		width = GetColumnWidth(i);
		xpos += width;
		if (pos.x < xpos)
			break;
	}
	if (i >= ncols)
		return false;
	*row = r;
	*column = i;
	return true;
}

void
MyListCtrl::StartEditText(int row, int column)
{
	wxRect rect;
	int x0, x1, dx, size, xpos, ypos, xunit, yunit, yorigin;
	if (editText != NULL && editText->IsShown())
		EndEditText(true);
	if (dataSource == NULL || !dataSource->IsItemEditable(this, row, column))
		return;

	/*  Scroll the list so that the editing item is visible  */
	EnsureVisible(row);
	GetItemRectForRowAndColumn(rect, row, column);
	rect.Inflate(1, 2);
	editRow = row;
	editColumn = column;
	GetScrollPixelsPerUnit(&xunit, &yunit);
	xpos = GetScrollPos(wxHORIZONTAL);
	ypos = GetScrollPos(wxVERTICAL);
	x0 = rect.GetX();
	x1 = x0 + rect.GetWidth();
	if (x0 < 0) {
		/*  Scroll right  */
		dx = x0 / xunit - 1;
		if (xpos + dx < 0)
			dx = -xpos;
	} else if (x1 > (size = GetSize().GetWidth())) {
		/*  Scroll left  */
		dx = ((x1 - size) / xunit) + 1;
	} else dx = 0;
	if (dx != 0) {
		/*  m_mainWin is a protected member in wxGenericListCtrl  */
		((wxScrolledWindow *)m_mainWin)->Scroll(xpos + dx, -1);
		Refresh();
	}

	/*  Reposition the rect relative to the origin of the scrolling area  */
	yorigin = ((wxScrolledWindow *)m_mainWin)->GetPosition().y;
	GetItemRectForRowAndColumn(rect, row, column);
	rect.Inflate(1, 2);
	rect.Offset(0, -yorigin);

	wxString str = dataSource->GetItemText(this, editRow, editColumn);
	if (editText == NULL) {
		editText = new wxTextCtrl(((wxScrolledWindow *)m_mainWin), -1, wxT(""), rect.GetPosition(), rect.GetSize(), wxTE_PROCESS_ENTER | wxTE_PROCESS_TAB);
		editText->Connect(wxID_ANY, wxEVT_KEY_DOWN, wxKeyEventHandler(MyListCtrl::OnKeyDownOnEditText), NULL, this);
		editText->Connect(wxID_ANY, wxEVT_KILL_FOCUS, wxFocusEventHandler(MyListCtrl::OnKillFocusOnEditText), NULL, this);
	} else {
		editText->SetSize(rect);
		editText->Clear();
		editText->Show();
	}
	editText->AppendText(str);
	editText->SetFocus();
//	editText->Connect(wxID_ANY, wxEVT_IDLE, wxIdleEventHandler(MyListCtrl::OnIdle), NULL, this);

	editText->SetSelection(-1, -1);  //  Select all text

	/*  Select only this row  */
	x1 = GetItemCount();
	for (x0 = 0; x0 < x1; x0++)
		SetItemState(x0, (x0 == row ? wxLIST_STATE_SELECTED : 0), wxLIST_STATE_SELECTED);

	//  Call the event handler directly
	//  (Otherwise, the table selection may be updated from the "current" selection in the molecule
	//  before the selection is updated from the table)
	wxCommandEvent dummyEvent;
	OnTableSelectionChanged(dummyEvent);
}

void
MyListCtrl::EndEditTextAndRestart(bool setValueFlag, int newRow, int newColumn)
{
	wxString sval;
	if (editText != NULL && editText->IsShown()) {
		if (setValueFlag && dataSource) {
			sval = editText->GetValue();
		}
		if (wxWindow::FindFocus() == editText) {
			SetFocus();
		}
#if defined(__WXMAC__)
		{
			/*  Erase the focus ring  */
			wxRect rect = editText->GetRect();
			rect = rect.Inflate(5, 5);
			//	Refresh(true, &rect);  /*  This somehow leaves lower side of the focus ring to remain  */
			Refresh();
		}
#endif
		editText->Hide();  /*  Temporarily hide until new editing starts  */

		if (setValueFlag && dataSource)
			dataSource->SetItemText(this, editRow, editColumn, sval);

	}
	
	if (newRow >= 0 && newColumn >= 0) {
		StartEditText(newRow, newColumn);
	} else {
		editRow = editColumn = -1;
#if defined(__WXMAC__)
		if (editText != NULL) {
			editText->Disconnect(wxID_ANY);
			editText->Destroy();
			editText = NULL;
		}
#else
		if (editText != NULL) {
			editText->Move(-1000, -1000);
			editText->Hide();
		}
#endif
	}

}

void
MyListCtrl::EndEditText(bool setValueFlag)
{
	EndEditTextAndRestart(setValueFlag, -1, -1);
}

void
MyListCtrl::OnKillFocusOnEditText(wxFocusEvent &event)
{
	if (editText != NULL && editText->IsShown()) {
		EndEditText(true);
	}
}

void
MyListCtrl::OnIdle(wxIdleEvent &event)
{
	/*
	wxWindow *wp;
	if (editText != NULL && (wp = wxWindow::FindFocus()) != editText) {
		EndEditText(true);
	}
	 */
}

void
MyListCtrl::OnKeyDownOnEditText(wxKeyEvent &event)
{
	int keyCode, ncols, nrows, ecol, erow;
	bool shiftDown;
	if (editText == NULL || !editText->IsShown()) {
		event.Skip();
		return;
	}
	keyCode = event.GetKeyCode();
	ncols = GetColumnCount();
	nrows = GetItemCount();
	shiftDown = (event.GetModifiers() == wxMOD_SHIFT);
	switch (keyCode) {
		case WXK_TAB:
			ecol = editColumn;
			erow = editRow;
			while (1) {
				if (shiftDown) {
					if (ecol == 0) {
						if (erow == 0)
							return;
						ecol = ncols - 1;
						erow--;
					} else {
						ecol--;
					}
				} else {
					if (ecol == ncols - 1) {
						if (erow >= nrows - 1)
							return;
						ecol = 0;
						erow++;
					} else {
						ecol++;
					}
				}
				if (dataSource == NULL || dataSource->IsItemEditable(this, erow, ecol))
					break;
			}
			EndEditTextAndRestart(true, erow, ecol);
			break;
		case WXK_RETURN:
			if (event.GetModifiers() == wxMOD_ALT) {
				printf("alt-return pressed\n"); fflush(stdout);
				EndEditText(true);
				printf("EndEditText completed\n"); fflush(stdout);
				return;
			}
			ecol = editColumn;
			erow = editRow;
			while (1) {
				if (shiftDown) {
					if (erow == 0)
						return;
					erow--;
				} else {
					if (erow == nrows - 1)
						return;
					erow++;
				}
				if (dataSource == NULL || dataSource->IsItemEditable(this, erow, ecol))
					break;
			}
			EndEditTextAndRestart(true, erow, ecol);
			break;
		case WXK_ESCAPE:
			EndEditText(false);
			break;
		default:
			event.Skip();
			break;
	}
}

bool
MyListCtrl::DeleteColumn(int col)
{
	EndEditText(false);
	return wxGenericListCtrl::DeleteColumn(col);
}

bool
MyListCtrl::InsertColumn(long col, const wxString &heading, int format, int width)
{
	EndEditText(false);
	return wxGenericListCtrl::InsertColumn(col, heading, format, width);
}

void
MyListCtrl::OnPopUpMenuSelected(wxCommandEvent &event)
{
	if (dataSource != NULL)
		dataSource->OnPopUpMenuSelected(this, lastPopUpRow, lastPopUpColumn, event.GetId());
}

void
MyListCtrl::OnLeftDClick(wxMouseEvent &event)
{
	int row, col;
	wxPoint pos = event.GetPosition();
	if (!FindItemAtPosition(pos, &row, &col))
		return;
	if (editText != NULL) {
		if (editRow == row && editColumn == col) {
			event.Skip();
			return;
		}
		EndEditTextAndRestart(true, row, col);
	} else {
		StartEditText(row, col);
	}
}

void
MyListCtrl::OnMouseDown(wxMouseEvent &event)
{
	int row, col, i, n;
	char **items;

	if (editText != NULL && editText->IsShown()) {
		//  During the text edit, mouse down outside the textctrl will terminate the editing
		EndEditText();
	}

	wxPoint pos = event.GetPosition();
	if (FindItemAtPosition(pos, &row, &col) && dataSource != NULL && (n = dataSource->HasPopUpMenu(this, row, col, &items)) > 0) {
		wxMenu mnu;
		for (i = 0; i < n; i++) {
			wxString itemStr(items[i], WX_DEFAULT_CONV);
			mnu.Append(i, itemStr);
			free(items[i]);
			items[i] = NULL;
		}
		free(items);
		lastPopUpColumn = col;
		lastPopUpRow = row;
		mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MyListCtrl::OnPopUpMenuSelected), NULL, this);
		PopupMenu(&mnu);
		return;
	}
	
	//  Intercept mouse down event and post selection change notification
	//  (a workaround of wxMSW problem where EVT_LIST_ITEM_SELECTED is not sent in some occasions)
	PostSelectionChangeNotification();
	event.Skip();
}

void
MyListCtrl::OnChar(wxKeyEvent &event)
{
	//  See comments on OnMouseUp()
	PostSelectionChangeNotification();
	event.Skip();
}

void
MyListCtrl::OnItemSelectionChanged(wxListEvent &event)
{
	PostSelectionChangeNotification();
}

void
MyListCtrl::OnTableSelectionChanged(wxCommandEvent &event)
{
	selectionChangeNotificationSent = false;
	if (dataSource == NULL)
		return;
	dataSource->OnSelectionChanged(this);
}

void
MyListCtrl::OnEnableTableSelectionNotification(wxCommandEvent &event)
{
	selectionChangeNotificationEnabled = true;
}
