/*
 *  ConsoleFrame.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 08/10/27.
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

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "ConsoleFrame.h"

#include "wx/menu.h"
#include "wx/regex.h"
#include "wx/colour.h"
#include "wx/sizer.h"

#include "MyApp.h"
#include "../MolLib/Ruby_bind/Molby_extern.h"
#include "MyMBConv.h"

BEGIN_EVENT_TABLE(ConsoleFrame, wxMDIChildFrame)
	EVT_UPDATE_UI(-1, ConsoleFrame::OnUpdateUI)
	EVT_CLOSE(ConsoleFrame::OnCloseWindow)
	EVT_MENU(wxID_CLOSE, ConsoleFrame::OnClose)
/*	EVT_MENU(wxID_CUT, ConsoleFrame::OnCut)
	EVT_MENU(wxID_COPY, ConsoleFrame::OnCopy)
	EVT_MENU(wxID_PASTE, ConsoleFrame::OnPaste)
	EVT_MENU(wxID_CLEAR, ConsoleFrame::OnClear) */
	EVT_MENU(wxID_UNDO, ConsoleFrame::OnUndo)
	EVT_MENU(wxID_REDO, ConsoleFrame::OnRedo)
END_EVENT_TABLE()

ConsoleFrame::ConsoleFrame(wxMDIParentFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type):
	wxMDIChildFrame(parent, wxID_ANY, title, pos, size, type)
{
}

ConsoleFrame::~ConsoleFrame()
{
	wxGetApp().DocManager()->FileHistoryRemoveMenu(file_history_menu);
}

void
ConsoleFrame::OnCreate()
{
	//  Make a text view
	int width, height;

	GetClientSize(&width, &height);
	textCtrl = new wxTextCtrl(this, wxID_ANY, _T(""), wxPoint(0, 0), wxSize(100, 100), wxTE_MULTILINE | wxTE_RICH);

	wxBoxSizer *consoleSizer = new wxBoxSizer(wxHORIZONTAL);
	consoleSizer->Add(textCtrl, 1, wxEXPAND);
	this->SetSizer(consoleSizer);
	consoleSizer->SetSizeHints(this);
	
	//  Set the default font (fixed-pitch)
	wxTextAttr attr;
#if defined(__WXMSW__)
	wxFont font = wxSystemSettings::GetFont(wxSYS_ANSI_FIXED_FONT);
#else
	wxFont font(11, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
#endif
	attr.SetFont(font);
	textCtrl->SetDefaultStyle(attr);
	
	//  Connect "OnKeyDown" event handler
	textCtrl->Connect(-1, wxEVT_KEY_DOWN, wxKeyEventHandler(ConsoleFrame::OnKeyDown), NULL, this);
	
	wxMenuBar *menu_bar = wxGetApp().CreateMenuBar(2, &file_history_menu, &edit_menu);
	
	//// Associate the menu bar with the frame
	SetMenuBar(menu_bar);
}

ConsoleFrame *
ConsoleFrame::CreateConsoleFrame(wxMDIParentFrame *parent)
{
#ifdef __WXMSW__
	wxPoint origin(0, 0);
	wxSize size(640, 200);
#else
	wxPoint origin(10, 24);
	wxSize size(640, 200);
#endif
	ConsoleFrame *frame = new ConsoleFrame(parent, _T("Console"), origin, wxDefaultSize, wxDEFAULT_FRAME_STYLE | wxNO_FULL_REPAINT_ON_RESIZE);

	frame->OnCreate();

	if (wxGetApp().IsFilterMode()) {
#if defined(__WXMSW__)
		frame->Maximize();
#else
		frame->SetSize(800, 480);
#endif
	} else {
		frame->SetSize(640, 200);
	}
		
	return frame;
}

bool
GetLineIncludingPosition(wxTextCtrl *ctrl, int pos, int *start, int *end)
{
	int pos1, pos2, posend;
	wxChar posChar;
	
	if (ctrl == NULL)
		return false;
	if (pos == 0)
		pos1 = 0;
	else {
		pos1 = pos;
		while (pos1 > 0) {
			posChar = ctrl->GetRange(pos1 - 1, pos1).GetChar(0);
			if (posChar == '\n')
				break;
			pos1--;
		}
	}
	posend = ctrl->GetLastPosition();
	posChar = ctrl->GetRange(posend - 1, posend).GetChar(0);
	if (pos >= posend)
		pos2 = pos;
	else {
		pos2 = pos;
		while (pos2 < posend) {
			posChar = ctrl->GetRange(pos2, pos2 + 1).GetChar(0);
			if (posChar == '\n') {
				pos2++;
				break;
			}
			pos2++;
		}
	}
	if (start != NULL)
		*start = pos1;
	if (end != NULL)
		*end = pos2;
	return true;
}

void
ConsoleFrame::OnCloseWindow(wxCloseEvent &event)
{
	//  Do not delete this window; it may be reopened later
	this->Hide();
}

void
ConsoleFrame::OnClose(wxCommandEvent &event)
{
	this->Close();
}

//  I do not understand why these functions should be written...
//  Certainly there should be better way to implement these.
void
ConsoleFrame::OnUndo(wxCommandEvent &event)
{
	if (wxWindow::FindFocus() == textCtrl)
		textCtrl->Undo();
	else event.Skip();
}

void
ConsoleFrame::OnRedo(wxCommandEvent &event)
{
	if (wxWindow::FindFocus() == textCtrl)
		textCtrl->Redo();
	else event.Skip();
}

void
ConsoleFrame::OnUpdateUI(wxUpdateUIEvent& event)
{
	int uid = event.GetId();
	if (uid == wxID_CLOSE)
		//  Why this is not automatically done??
		event.Enable(true);
	else if (uid == wxID_UNDO || uid == wxID_REDO) {
		if (wxWindow::FindFocus() == textCtrl) {
			if (uid == wxID_UNDO)
				event.Enable(textCtrl->CanUndo());
			else
				event.Enable(textCtrl->CanRedo());
		} else event.Skip();
	} else event.Skip();
}

void
ConsoleFrame::OnEnterPressed(wxKeyEvent& event)
{
	if (::wxGetKeyState(WXK_ALT)) {
		textCtrl->WriteText(wxT("\n> "));
		return;
	}
	
	int start, pos, end, veryend, lastpos;
	wxChar startChar;

	//  Get the block of script to be executed
	pos = textCtrl->GetInsertionPoint();
	lastpos = textCtrl->GetLastPosition();
	veryend = -1;
	while (pos >= 0) {
		if (!GetLineIncludingPosition(textCtrl, pos, &start, &end) || start == end) {
			start = end = veryend = pos;
			break;
		}
		if (veryend < 0)
			veryend = end;
		startChar = textCtrl->GetRange(start, start + 1).GetChar(0);
		if (startChar == '%') {
			start++;
			break;
		} else if (startChar == '>') {
			pos = start - 1;
			continue;
		} else {
			start = end = veryend = pos;
			break;
		}
	}
	while (start < end && veryend < lastpos) {
		pos = veryend + 1;
		if (!GetLineIncludingPosition(textCtrl, pos, &pos, &end) || pos == end) {
			break;
		}
		startChar = textCtrl->GetRange(pos, pos + 1).GetChar(0);
		if (startChar != '>')
			break;
		veryend = end;
	}

	wxString string = textCtrl->GetRange(start, veryend);
	int len = string.Len();
	
	//  Is there any non-whitespace characters?
	wxChar ch;
	int i;
	for (i = 0; i < len; i++) {
		ch = string[i];
		if (ch != ' ' && ch != '\t' && ch != '\n' && ch != 'r')
			break;
	}
	if (i < len) {
		//  Input is not empty
		if (veryend < lastpos) {
			// Enter is pressed in the block not at the end
			// -> Insert the text at the end
			wxRegEx re1(wxT("[ \t\n]*$"));
			re1.ReplaceFirst(&string, wxT(""));
			wxRegEx re2(wxT("^[ \t\n]*"));
			re2.ReplaceFirst(&string, wxT(""));
			if (textCtrl->GetRange(lastpos - 1, lastpos).GetChar(0) != '\n')
				textCtrl->AppendText(wxT("\n"));
			MyAppCallback_showRubyPrompt();
			MyAppCallback_setConsoleColor(3);
			textCtrl->AppendText(string);
		} else {
			wxTextAttr scriptAttr(*wxBLUE);
			textCtrl->SetStyle(start, veryend, scriptAttr);
		}
		string.Append(wxT("\n"));  //  To avoid choking Ruby interpreter
		wxRegEx re3(wxT("\n>"));
		re3.Replace(&string, wxT("\n"));
		if (textCtrl->GetRange(lastpos - 1, lastpos).GetChar(0) != '\n')
			textCtrl->AppendText(wxT("\n"));
		MyAppCallback_setConsoleColor(0);
		textCtrl->Update();
		
		//  Invoke ruby interpreter
		int status;
		RubyValue val;
		Molecule *mol = MoleculeCallback_currentMolecule();
		MoleculeLock(mol);
		val = Molby_evalRubyScriptOnMolecule(string.mb_str(WX_DEFAULT_CONV), MoleculeCallback_currentMolecule(), NULL, &status);
		MoleculeUnlock(mol);
		if (status != -1) {  /*  Status -1 is already handled  */
			MyAppCallback_setConsoleColor(1);
			if (status != 0) {
				Molby_showError(status);
			} else {
				textCtrl->AppendText(wxT("-->"));
				Molby_showRubyValue(val);
			}
			MyAppCallback_setConsoleColor(0);
			textCtrl->AppendText(wxT("\n"));
			MyAppCallback_showRubyPrompt();
		}
	} else {
		textCtrl->AppendText(wxT("\n"));
		MyAppCallback_showRubyPrompt();
	}
}

void
ConsoleFrame::OnKeyDown(wxKeyEvent &event)
{
	int code = event.GetKeyCode();
	//	printf("OnChar: %d\n", code);
	if (code == WXK_RETURN || code == WXK_NUMPAD_ENTER)
		OnEnterPressed(event);
	else
		event.Skip();
}

void
ConsoleFrame::EmptyBuffer(bool showRubyPrompt)
{
	textCtrl->Clear();
	if (showRubyPrompt)
		MyAppCallback_showRubyPrompt();
}
