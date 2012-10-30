/*
 *  ConsoleFrame.h
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

#ifndef __ConsoleFrame_h__
#define __ConsoleFrame_h__

#include "wx/mdi.h"
//#include "wx/richtext/richtextctrl.h"
#include "wx/textctrl.h"

class wxMenu;

class ConsoleFrame: public wxMDIChildFrame
{

public:
	wxTextCtrl *textCtrl;
	wxMenu *file_history_menu;
	wxMenu *edit_menu;

	wxFont *default_font;

	ConsoleFrame(wxMDIParentFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type);
	virtual ~ConsoleFrame();

	void OnCreate();
	void OnEnterPressed(wxKeyEvent& event);
	void OnKeyDown(wxKeyEvent &event);

	static ConsoleFrame *CreateConsoleFrame(wxMDIParentFrame *parent);
	void OnCloseWindow(wxCloseEvent &event);
	void OnClose(wxCommandEvent &event);
	void OnUpdateUI(wxUpdateUIEvent& event);

	void OnUndo(wxCommandEvent &event);
	void OnRedo(wxCommandEvent &event);
	void EmptyBuffer(bool showRubyPrompt = true);

private:
	DECLARE_EVENT_TABLE()
};

#endif /* __ConsoleFrame_h__ */
