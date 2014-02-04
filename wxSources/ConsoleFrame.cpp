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
#include "../MolLib/MolLib.h"
#include "../MolLib/Missing.h"
#include "MyMBConv.h"

#if defined(__WXMSW__)
#include <windows.h>
#include <richedit.h>
#endif

BEGIN_EVENT_TABLE(ConsoleFrame, wxFrame)
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

ConsoleFrame::ConsoleFrame(wxFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type):
	wxFrame(parent, wxID_ANY, title, pos, size, type)
{
	valueHistory = commandHistory = NULL;
	nValueHistory = nCommandHistory = 0;
	valueHistoryIndex = commandHistoryIndex = -1;
	keyInputPos = -1;
	selectionFrom = selectionTo = -1;
}

ConsoleFrame::~ConsoleFrame()
{
	int i;
	for (i = 0; i < nValueHistory; i++)
		free(valueHistory[i]);
	free(valueHistory);
	for (i = 0; i < nCommandHistory; i++)
		free(commandHistory[i]);
	free(commandHistory);	
//	wxGetApp().DocManager()->FileHistoryRemoveMenu(file_history_menu);
}

void
ConsoleFrame::OnCreate()
{
	//  Make a text view
	int width, height;

	GetClientSize(&width, &height);
	textCtrl = new wxTextCtrl(this, wxID_ANY, _T(""), wxPoint(0, 0), wxSize(100, 100), wxTE_MULTILINE | wxTE_RICH | wxTE_PROCESS_ENTER
							  );

#if defined(__WXMSW__)
	{
		HWND hwnd = (HWND)(textCtrl->GetHWND());
		DWORD dwOptions;
		LPARAM result;
		
		/*  Disable dual language mode  */
		dwOptions = SendMessage(hwnd, EM_GETLANGOPTIONS, 0, 0); 
		dwOptions &= ~0x02;   /*  0x02 = IMF_AUTOFONT  */
		result = SendMessage(hwnd, EM_SETLANGOPTIONS, 0, (LPARAM)dwOptions);
		printf("%ld\n", (long)result);		
	}
#endif
		
	wxBoxSizer *consoleSizer = new wxBoxSizer(wxHORIZONTAL);
	consoleSizer->Add(textCtrl, 1, wxEXPAND);
	this->SetSizer(consoleSizer);
	consoleSizer->SetSizeHints(this);
	
	//  Set the default font (fixed-pitch)
#if defined(__WXMSW__)
	default_font = new wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
	if (!default_font->SetFaceName(wxT("Consolas")))
		default_font->SetFaceName(wxT("Courier"));
#elif defined(__WXMAC__)
	default_font = new wxFont(11, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
	default_font->SetFaceName(wxT("Monaco"));
#else
	default_font = new wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
#endif
	current_attr = new wxTextAttr;
	current_attr->SetFont(*default_font);
	textCtrl->SetDefaultStyle(*current_attr);

	//  Connect textCtrl event handler
	textCtrl->Connect(-1, wxEVT_KEY_DOWN, wxKeyEventHandler(ConsoleFrame::OnKeyDown), NULL, this);
	textCtrl->Connect(-1, wxEVT_TEXT_ENTER, wxCommandEventHandler(ConsoleFrame::OnTextEnter), NULL, this);
	textCtrl->Connect(-1, wxEVT_SET_FOCUS, wxFocusEventHandler(ConsoleFrame::OnSetFocus), NULL, this);
	textCtrl->Connect(-1, wxEVT_KILL_FOCUS, wxFocusEventHandler(ConsoleFrame::OnKillFocus), NULL, this);

	wxMenuBar *menu_bar = wxGetApp().CreateMenuBar(2, &file_history_menu, &edit_menu);
	
	//// Associate the menu bar with the frame
	SetMenuBar(menu_bar);
}

ConsoleFrame *
ConsoleFrame::CreateConsoleFrame(wxFrame *parent)
{
#ifdef __WXMSW__
	wxPoint origin(0, 0);
	wxSize size(640, 200);
#else
	wxPoint origin(10, 24);
	wxSize size(640, 200);
#endif
	ConsoleFrame *frame = new ConsoleFrame(parent, _T("Molby Console"), origin, wxDefaultSize, wxDEFAULT_FRAME_STYLE | wxNO_FULL_REPAINT_ON_RESIZE);

	frame->OnCreate();
	frame->SetClientSize(size);
	
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
	//  Check if all windows are gone
	wxGetApp().CheckIfAllWindowsAreGone(NULL);
}

void
ConsoleFrame::OnClose(wxCommandEvent &event)
{
	this->Close();
}

void
ConsoleFrame::OnKillFocus(wxFocusEvent &event)
{
#if defined(__WXMSW__)
	textCtrl->GetSelection(&selectionFrom, &selectionTo);
#endif
	event.Skip();
}

void
ConsoleFrame::OnSetFocus(wxFocusEvent &event)
{
#if defined(__WXMSW__)
	if (selectionFrom >= 0 && selectionTo >= 0) {
		textCtrl->SetSelection(selectionFrom, selectionTo);
	}
#endif
	event.Skip();
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
ConsoleFrame::OnEnterPressed()
{
	int start, pos, end, veryend, lastpos;
	wxChar startChar;
	
	if (::wxGetKeyState(WXK_ALT)) {
		textCtrl->WriteText(wxT("\n> "));
		return;
	}
	
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
			wxTextAttr scriptAttr(*wxBLUE, wxNullColour, *default_font);
			textCtrl->SetStyle(start, veryend, scriptAttr);
		}
		string.Append(wxT("\n"));  //  To avoid choking Ruby interpreter
		wxRegEx re3(wxT("\n>"));
		re3.Replace(&string, wxT("\n"));
		ch = textCtrl->GetRange(lastpos - 1, lastpos).GetChar(0);
		if (ch != '\n')
			textCtrl->AppendText(wxT("\n"));
		textCtrl->Update();
		MyAppCallback_setConsoleColor(0);
		
		//  Invoke ruby interpreter
		int status;
		RubyValue val;
		char *script;
		Molecule *mol = MoleculeCallback_currentMolecule();
		int lock_needed = (mol != NULL && mol->mutex != NULL);
		if (lock_needed)
			MoleculeLock(mol);
		script = strdup(string.mb_str(WX_DEFAULT_CONV));
		val = Molby_evalRubyScriptOnMolecule(script, MoleculeCallback_currentMolecule(), NULL, &status);
		if (lock_needed)
			MoleculeUnlock(mol);
		script[strlen(script) - 1] = 0;  /*  Remove the last newline  */
		AssignArray(&commandHistory, &nCommandHistory, sizeof(char *), nCommandHistory, &script);
		if (nCommandHistory >= MAX_HISTORY_LINES)
			DeleteArray(&commandHistory, &nCommandHistory, sizeof(char *), 0, 1, NULL);
		if (status != -1) {  /*  Status -1 is already handled  */
			char *valueString;
			MyAppCallback_setConsoleColor(1);
			if (status != 0) {
				Molby_showError(status);
			} else {
				textCtrl->AppendText(wxT("-->"));
				Molby_showRubyValue(val, &valueString);
				AssignArray(&valueHistory, &nValueHistory, sizeof(char *), nValueHistory, &valueString);
				if (nValueHistory >= MAX_HISTORY_LINES)
					DeleteArray(&valueHistory, &nValueHistory, sizeof(char *), 0, 1, NULL);
			}
			MyAppCallback_setConsoleColor(0);
			textCtrl->AppendText(wxT("\n"));
			MyAppCallback_showRubyPrompt();
		}
	} else {
		textCtrl->AppendText(wxT("\n"));
		MyAppCallback_showRubyPrompt();
	}
	commandHistoryIndex = valueHistoryIndex = -1;
}

void
ConsoleFrame::ShowHistory(bool up, bool option)
{
	char *p;
	if (commandHistoryIndex == -1 && valueHistoryIndex == -1) {
		if (!up)
			return;
		historyPos = textCtrl->GetLastPosition();
	}
	if (option) {
		if (up) {
			if (valueHistoryIndex == -1) {
				if (nValueHistory == 0)
					return;
				valueHistoryIndex = nValueHistory;
			}
			if (valueHistoryIndex <= 0)
				return; /* Do nothing */
			valueHistoryIndex--;
			p = valueHistory[valueHistoryIndex];
		} else {
			if (valueHistoryIndex == -1)
				return;  /*  Do nothing  */
			if (valueHistoryIndex == nValueHistory - 1) {
				valueHistoryIndex = -1;
				p = "";
			} else {
				valueHistoryIndex++;
				p = valueHistory[valueHistoryIndex];
			}
		}
	} else {
		if (up) {
			if (commandHistoryIndex == -1) {
				if (nCommandHistory == 0)
					return;
				commandHistoryIndex = nCommandHistory;
			}
			if (commandHistoryIndex <= 0)
				return; /* Do nothing */
			commandHistoryIndex--;
			p = commandHistory[commandHistoryIndex];
		} else {
			if (commandHistoryIndex == -1)
				return;  /*  Do nothing  */
			if (commandHistoryIndex == nCommandHistory - 1) {
				commandHistoryIndex = -1;
				p = "";
			} else {
				commandHistoryIndex++;
				p = commandHistory[commandHistoryIndex];
			}
		}
	}
	textCtrl->Replace(historyPos, textCtrl->GetLastPosition(), wxT(""));
	SetConsoleColor(option ? 1 : 3);
	while (isspace(*p))
		p++;
	textCtrl->AppendText(wxString(p, WX_DEFAULT_CONV));
	SetConsoleColor(0);
}

void
ConsoleFrame::OnKeyDown(wxKeyEvent &event)
{
	int code = event.GetKeyCode();

#if __WXOSX_COCOA__
	//  in wxCocoa, the key down event is fired even when the input method still has
	//  marked text. We need to avoid doing our own work when the input method is expecting
	//  key events.
	extern bool textCtrl_hasMarkedText();
	if (textCtrl_hasMarkedText()) {
		event.Skip();
		return;
	}
#endif
	
	if (code == WXK_RETURN || code == WXK_NUMPAD_ENTER)
		OnEnterPressed();
	else if ((code == WXK_UP || code == WXK_DOWN) && (textCtrl->GetInsertionPoint() == textCtrl->GetLastPosition()))
		ShowHistory(code == WXK_UP, event.GetModifiers() == wxMOD_ALT);
	else
		event.Skip();
}

void
ConsoleFrame::OnTextEnter(wxCommandEvent &event)
{
	OnEnterPressed();
}

int
ConsoleFrame::AppendConsoleMessage(const char *mes)
{
	wxString string(mes, WX_DEFAULT_CONV);
	textCtrl->AppendText(string);
	commandHistoryIndex = valueHistoryIndex = -1;
	return string.Len();
}

void
ConsoleFrame::FlushConsoleMessage()
{
	textCtrl->Refresh();
	textCtrl->Update();
}

void
ConsoleFrame::SetConsoleColor(int color)
{
	static const wxColour *col[4];
	if (col[0] == NULL) {
		col[0] = wxBLACK;
		col[1] = wxRED;
		col[2] = wxGREEN;
		col[3] = wxBLUE;
	}
	current_attr->SetTextColour(*col[color % 4]);
	current_attr->SetFont(*default_font);
	textCtrl->SetDefaultStyle(wxTextAttr(*col[color % 4], wxNullColour, *default_font));
}

void
ConsoleFrame::EmptyBuffer(bool showRubyPrompt)
{
	textCtrl->Clear();
	if (showRubyPrompt)
		MyAppCallback_showRubyPrompt();
	commandHistoryIndex = valueHistoryIndex = -1;
}

#pragma mark ====== Plain-C interface ======

int
MyAppCallback_showScriptMessage(const char *fmt, ...)
{
	if (fmt != NULL) {
		char *p;
		va_list ap;
		int retval;
		va_start(ap, fmt);
		if (strchr(fmt, '%') == NULL) {
			/*  No format characters  */
			return wxGetApp().GetConsoleFrame()->AppendConsoleMessage(fmt);
		} else if (strcmp(fmt, "%s") == 0) {
			/*  Direct output of one string  */
			p = va_arg(ap, char *);
			return wxGetApp().GetConsoleFrame()->AppendConsoleMessage(p);
		}
		vasprintf(&p, fmt, ap);
		if (p != NULL) {
			retval = wxGetApp().GetConsoleFrame()->AppendConsoleMessage(p);
			free(p);
			return retval;
		} else return 0;
	} else {
		wxGetApp().GetConsoleFrame()->FlushConsoleMessage();
		return 0;
	}
	return 0;
}

void
MyAppCallback_setConsoleColor(int color)
{
	wxGetApp().GetConsoleFrame()->SetConsoleColor(color);
}

void
MyAppCallback_showRubyPrompt(void)
{
	MyAppCallback_setConsoleColor(0);
	MyAppCallback_showScriptMessage("%% ");
}
