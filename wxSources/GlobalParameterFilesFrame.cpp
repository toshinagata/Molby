/*
 *  GlobalParameterFilesFrame.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 09/11/14.
 *  Copyright 2009 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "GlobalParameterFilesFrame.h"

#include "wx/menu.h"
#include "wx/regex.h"
#include "wx/colour.h"
#include "wx/button.h"
#include "wx/filedlg.h"
#include "wx/msgdlg.h"

#include "MyApp.h"
#include "MyListCtrl.h"
#include "../MolLib/MolLib.h"
#include "../MolLib/Ruby_bind/Molby_extern.h"
#include "../MolLib/Missing.h"
#include "MyMBConv.h"

// #include "../MolLib/Ruby_bind/Molby_extern.h"

enum {
	myID_addFileButton = 500,
	myID_removeFileButton
};

BEGIN_EVENT_TABLE(GlobalParameterFilesFrame, wxFrame)
EVT_CLOSE(GlobalParameterFilesFrame::OnCloseWindow)
EVT_MENU(wxID_CLOSE, GlobalParameterFilesFrame::OnClose)
EVT_UPDATE_UI(wxID_CLOSE, GlobalParameterFilesFrame::OnUpdateUI)
EVT_BUTTON(myID_addFileButton, GlobalParameterFilesFrame::OnAddGlobalParameterFile)
EVT_BUTTON(myID_removeFileButton, GlobalParameterFilesFrame::OnRemoveGlobalParameterFile)
END_EVENT_TABLE()

#pragma mark ====== Static utility functions ======
									
GlobalParameterFilesFrame::GlobalParameterFilesFrame(wxWindow *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type):
wxFrame(parent, wxID_ANY, title, pos, size, type)
{
}

GlobalParameterFilesFrame::~GlobalParameterFilesFrame()
{
	wxGetApp().DocManager()->FileHistoryRemoveMenu(file_history_menu);
}

void
GlobalParameterFilesFrame::OnCreate()
{
	//  Frame structure:
	//  vertical{ listview, horizontal{ plus_button, minus_button } }
	wxSize size(340, 160);
	wxBoxSizer *sizer0;
	sizer0 = new wxBoxSizer(wxVERTICAL);
	listctrl = new MyListCtrl();
	listctrl->Create(this, wxID_ANY, wxDefaultPosition, size);
	listctrl->InsertColumn(0, _T("name"), wxLIST_FORMAT_LEFT, 80);
	listctrl->InsertColumn(1, _T("directory"), wxLIST_FORMAT_LEFT, 240);
	sizer0->Add(listctrl, 1, wxLEFT | wxRIGHT | wxEXPAND, 0);
	{
		wxBoxSizer *sizer1;
		sizer1 = new wxBoxSizer(wxHORIZONTAL);
		add_button = new wxButton(this, myID_addFileButton, _T("Load..."), wxDefaultPosition, wxDefaultSize);
		remove_button = new wxButton(this, myID_removeFileButton, _T("Unload"), wxDefaultPosition, wxDefaultSize);
		sizer1->Add(add_button, 0, wxALL | wxEXPAND, 3);
		sizer1->Add(remove_button, 0, wxALL | wxEXPAND, 3);
		sizer0->Add(sizer1, 0, wxALL | wxEXPAND, 8);
	}
	this->SetSizer(sizer0);
	listctrl->SetDataSource(this);
	remove_button->Enable(false);

	/*  Create menu bar  */
	wxMenuBar *menu_bar = wxGetApp().CreateMenuBar(2, &file_history_menu, &edit_menu);
	SetMenuBar(menu_bar);
}

GlobalParameterFilesFrame *
GlobalParameterFilesFrame::CreateGlobalParameterFilesFrame(wxWindow *parent)
{
#ifdef __WXMSW__
	wxPoint origin(16, 16);
	wxSize size(700, 240);
#else
	wxPoint origin(26, 40);
	wxSize size(700, 240);
#endif

	GlobalParameterFilesFrame *frame = new GlobalParameterFilesFrame(parent, _T("Load/Unload Global Parameters"), origin, size, wxDEFAULT_FRAME_STYLE | wxNO_FULL_REPAINT_ON_RESIZE);
	
	frame->OnCreate();
	return frame;
}

void
GlobalParameterFilesFrame::OnCloseWindow(wxCloseEvent &event)
{
	//  Do not delete this window; it may be reopened later
	this->Hide();
	//  Check if all windows are gone
	wxGetApp().CheckIfAllWindowsAreGone(NULL);
}

void
GlobalParameterFilesFrame::OnClose(wxCommandEvent &event)
{
	this->Close();
}

void
GlobalParameterFilesFrame::OnUpdateUI(wxUpdateUIEvent& event)
{
	int uid = event.GetId();
	if (uid == wxID_CLOSE)
		event.Enable(true);
}

void
GlobalParameterFilesFrame::OnAddGlobalParameterFile(wxCommandEvent& event)
{
	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Choose Parameter File"), _T(""), _T(""), _T("All files (*.*)|*.*"), wxFD_OPEN | wxFD_CHANGE_DIR | wxFD_FILE_MUST_EXIST);
	if (dialog->ShowModal() == wxID_OK) {
		char *p = strdup((const char *)(dialog->GetPath().mb_str(wxConvFile)));
		char *wbuf;
		int i;
		int src_idx = ParameterCommentIndexForGlobalFileName(p);
		for (i = gGlobalParInfo.count - 1; i >= 0; i--) {
			if (gGlobalParInfo.files[i].src == src_idx)
				break;
		}
		if (i >= 0) {
			char *s;
			asprintf(&s, "File '%s' has already been read. Do you want to replace the parameters?", ParameterGetComment(src_idx));
			wxString mes(s, WX_DEFAULT_CONV);
			if (::wxMessageBox(mes, _T("Reload global parameter file"), wxOK | wxCANCEL | wxICON_EXCLAMATION) != wxOK)
				return;
		}
		ParameterReadFromFile(NULL, p, &wbuf, NULL);
		if (wbuf != NULL) {
			MyAppCallback_setConsoleColor(1);
			MyAppCallback_showScriptMessage("%s", wbuf);
			MyAppCallback_setConsoleColor(0);
			free(wbuf);
		}
		free(p);
		/*  Need to update MD parameters for all molecules  */
		{
			Molecule *mol;
			for (i = 0; (mol = MoleculeCallback_moleculeAtIndex(i)) != NULL; i++) {
				mol->needsMDRebuild = 1;
			}
		}
	}
	dialog->Destroy();
	listctrl->RefreshTable();

	{
		MyListCtrl *parameterListCtrl = wxGetApp().GetGlobalParameterListCtrl();
		if (parameterListCtrl != NULL)
			parameterListCtrl->RefreshTable();
	}
}

void
GlobalParameterFilesFrame::OnRemoveGlobalParameterFile(wxCommandEvent& event)
{
	int i, n, cn;
	wxString files;
	ParFileInfo *ip;
	const char *cp;
	char *s;

	n = listctrl->GetItemCount();
	cn = 0;
	for (i = 0; i < n; i++) {
		if (i >= gGlobalParInfo.builtinCount && listctrl->GetItemState(i, wxLIST_STATE_SELECTED)) {
			ip = &gGlobalParInfo.files[i];
			cp = ParameterGetComment(ip->src);
			if (cn > 0)
				files << _T(", ");
			files << wxString(cp, WX_DEFAULT_CONV);
			cn++;
		}
	}
	if (cn == 0)
		return;
	
	asprintf(&s, "Do you really want to unload the parameter%s %s?", (cn > 1 ? "s" : ""), (const char *)files.mb_str(WX_DEFAULT_CONV));
	wxString mes(s, WX_DEFAULT_CONV);
	if (::wxMessageBox(mes, _T("Unload global parameter file"), wxOK | wxCANCEL | wxICON_EXCLAMATION) != wxOK)
		return;
	for (i = gGlobalParInfo.count - 1; i >= gGlobalParInfo.builtinCount; i--) {
		if (listctrl->GetItemState(i, wxLIST_STATE_SELECTED)) {
			ip = &gGlobalParInfo.files[i];
			ParameterDeleteAllEntriesForSource(gBuiltinParameters, ip->src);
			listctrl->SetItemState(i, 0, wxLIST_STATE_SELECTED);
		}
	}
	listctrl->RefreshTable();

	/*  Request to update MD parameters for all molecules  */
	{
		Molecule *mol;
		for (i = 0; (mol = MoleculeCallback_moleculeAtIndex(i)) != NULL; i++) {
			mol->needsMDRebuild = 1;
		}
	}
	
	{
		MyListCtrl *parameterListCtrl = wxGetApp().GetGlobalParameterListCtrl();
		if (parameterListCtrl != NULL)
			parameterListCtrl->RefreshTable();
	}
	
	remove_button->Enable(false);
}

#pragma mark ====== MyListCtrl data source ======

/*  Get information from gGlobalParInfo (defined in Parameter.c)  */
int
GlobalParameterFilesFrame::GetItemCount(MyListCtrl *ctrl)
{
	return gGlobalParInfo.count;
}

wxString
GlobalParameterFilesFrame::GetItemText(MyListCtrl *ctrl, long row, long column) const
{
	ParFileInfo *ip = &gGlobalParInfo.files[row];
	const char *p = NULL;
	if (column == 0)
		p = ParameterGetComment(ip->src);
	else if (column == 1) {
		if (row < gGlobalParInfo.builtinCount)
			p = "(built-in)";
		else
			p = ip->dir;
	}
	wxString str((p ? p : ""), WX_DEFAULT_CONV);
	return str;
}

int
GlobalParameterFilesFrame::SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value)
{
	return 0;
}

void
GlobalParameterFilesFrame::DragSelectionToRow(MyListCtrl *ctrl, long row)
{
}

bool
GlobalParameterFilesFrame::IsItemEditable(MyListCtrl *ctrl, long row, long column)
{
	return false;
}

bool
GlobalParameterFilesFrame::IsDragAndDropEnabled(MyListCtrl *ctrl, long row)
{
	return 0;
}

void
GlobalParameterFilesFrame::OnSelectionChanged(MyListCtrl *ctrl)
{
	int i, n, cn;
	ctrl->EnableSelectionChangeNotification(false);
	n = ctrl->GetItemCount();
	cn = 0;
	for (i = 0; i < n; i++) {
		if (ctrl->GetItemState(i, wxLIST_STATE_SELECTED)) {
			if (i < gGlobalParInfo.builtinCount)
				ctrl->SetItemState(i, 0, wxLIST_STATE_SELECTED);
			else cn++;
		}
	}
	ctrl->EnableSelectionChangeNotification(true);
	if (cn > 0)
		remove_button->Enable(true);
	else remove_button->Enable(false);
}

int
GlobalParameterFilesFrame::SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg)
{
	if (row >= 0 && row < gGlobalParInfo.builtinCount && col == -1) {
		fg[0] = fg[1] = fg[2] = 0.5;
		return 1;
	} else return 0;
}
