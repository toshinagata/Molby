/*
 *  GlobalParameterFrame.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 09/11/05.
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

#include "GlobalParameterFrame.h"

#include "wx/menu.h"
#include "wx/regex.h"
#include "wx/colour.h"

#include "MyApp.h"
#include "MyListCtrl.h"
#include "../MolLib/Ruby_bind/Molby_extern.h"
#include "MyMBConv.h"

BEGIN_EVENT_TABLE(GlobalParameterFrame, wxFrame)
	EVT_CLOSE(GlobalParameterFrame::OnCloseWindow)
	EVT_MENU(wxID_CLOSE, GlobalParameterFrame::OnClose)
	EVT_UPDATE_UI(wxID_CLOSE, GlobalParameterFrame::OnUpdateUI)
END_EVENT_TABLE()

GlobalParameterFrame::GlobalParameterFrame(wxWindow *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type):
	wxFrame(parent, wxID_ANY, title, pos, size, type)
{
}

GlobalParameterFrame::~GlobalParameterFrame()
{
	wxGetApp().DocManager()->FileHistoryRemoveMenu(file_history_menu);
}

void
GlobalParameterFrame::OnCreate()
{
	/*  Make a MyListCtrl view  */
	int width, height;
	GetClientSize(&width, &height);
	listCtrl = new MyListCtrl();
	listCtrl->Create(this, wxID_ANY, wxPoint(0, 0), wxSize(width, height));
	listCtrl->SetDataSource(this);
	wxMenuBar *menu_bar = wxGetApp().CreateMenuBar(2, &file_history_menu, &edit_menu);
	
	/*  Associate the menu bar with the frame  */
	SetMenuBar(menu_bar);
}

GlobalParameterFrame *
GlobalParameterFrame::CreateGlobalParameterFrame(wxWindow *parent)
{
#ifdef __WXMSW__
	wxPoint origin(16, 16);
	wxSize size(774, 300);
#else
	wxPoint origin(26, 40);
	wxSize size(774, 300);
#endif
	GlobalParameterFrame *frame = new GlobalParameterFrame(parent, _T("Global Parameters"), FromFrameDIP(parent, origin), FromFrameDIP(parent, size), wxDEFAULT_FRAME_STYLE | wxNO_FULL_REPAINT_ON_RESIZE);
	
	frame->OnCreate();
	return frame;
}

void
GlobalParameterFrame::OnCloseWindow(wxCloseEvent &event)
{
	//  Do not delete this window; it may be reopened later
	this->Hide();
	//  Check if all windows are gone
	wxGetApp().CheckIfAllWindowsAreGone(NULL);
}

void
GlobalParameterFrame::OnClose(wxCommandEvent &event)
{
	//  Why this is not automatically connected?
	this->Close();
}

void
GlobalParameterFrame::OnUpdateUI(wxUpdateUIEvent& event)
{
	//  Why this is not automatically done??
	int uid = event.GetId();
	if (uid == wxID_CLOSE)
		event.Enable(true);
}

#pragma mark ====== MyListCtrl data source ======

int
GlobalParameterFrame::GetItemCount(MyListCtrl *ctrl)
{
	return MainView_numberOfRowsInTable(NULL);
}

wxString
GlobalParameterFrame::GetItemText(MyListCtrl *ctrl, long row, long column) const
{
	char buf[128];
	MainView_valueForTable(NULL, column, row, buf, sizeof buf);
	wxString *str = new wxString(buf, WX_DEFAULT_CONV);
	return *str;
}

int
GlobalParameterFrame::SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value)
{
//	MainView_setValueForTable(NULL, column, row, value.mb_str(WX_DEFAULT_CONV));
	return 0;
}

void
GlobalParameterFrame::DragSelectionToRow(MyListCtrl *ctrl, long row)
{
}

bool
GlobalParameterFrame::IsItemEditable(MyListCtrl *ctrl, long row, long column)
{
	return false;
//	return MainView_isTableItemEditable(NULL, column, row);
}

bool
GlobalParameterFrame::IsDragAndDropEnabled(MyListCtrl *ctrl, long row)
{
	return 0;
}

void
GlobalParameterFrame::OnSelectionChanged(MyListCtrl *ctrl)
{
	MainView_setSelectionFromTable(NULL);
}

int
GlobalParameterFrame::SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg)
{
	if (col == -1) {
		int src = ParameterTableGetItemSource(gBuiltinParameters, row);
		if (src == -2) { /* separator row */
			bg[0] = bg[1] = bg[2] = 0.6;
			return 2;
		}
	}
	return 0;
}
