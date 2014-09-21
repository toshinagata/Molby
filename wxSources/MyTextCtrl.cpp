/*
 *  MyTextCtrl.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 2014/09/21.
 *  Copyright 2014 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "MyTextCtrl.h"

IMPLEMENT_DYNAMIC_CLASS(MyTextCtrl, wxTextCtrl)

BEGIN_EVENT_TABLE(MyTextCtrl, wxTextCtrl)
EVT_KEY_UP(MyTextCtrl::OnKeyUp)
END_EVENT_TABLE()

const wxEventType myTextCtrl_EVT_PROCESS_ESCAPE = wxNewEventType();

MyTextCtrl::MyTextCtrl() : wxTextCtrl()
{}

MyTextCtrl::MyTextCtrl(wxWindow *parent, wxWindowID id, const wxString &value,
					   const wxPoint &pos, const wxSize &size,
					   long style, const wxValidator &validator,
					   const wxString &name)
						: wxTextCtrl(parent, id, value, pos, size, style & ~(MyTextCtrl_Process_Escape), validator)
{
	processEscape = (style & MyTextCtrl_Process_Escape) != 0;
}

MyTextCtrl::~MyTextCtrl()
{}

//  Call escape key handler on _releasing_ the ESC key
//  Note: Overriding OnChar() function does not work!
//  (Escape key is swallowed somewhere between OnKeyDown and OnChar)
void
MyTextCtrl::OnKeyUp(wxKeyEvent &event) {
	if (processEscape && event.GetKeyCode() == WXK_ESCAPE) {
		wxCommandEvent event(myTextCtrl_EVT_PROCESS_ESCAPE, m_windowId);
		InitCommandEvent(event);
		if (HandleWindowEvent(event)) {
			return;
		}
	}
	event.Skip();
	return;
#if defined(__WXMSW__)
	//  OnKeyDown() is private in WXMSW, so we cannot call the superclass handler.
	//  Here we copy the code in wxTextCtrl::OnKeyDown()
	if ( event.GetModifiers() == wxMOD_CONTROL && IsRich() )
    {
        switch ( event.GetKeyCode() )
        {
            case 'C':
                Copy();
                return;
            case 'X':
                Cut();
                return;
            case 'V':
                Paste();
                return;
            default:
                break;
        }
    }
    if ( event.GetKeyCode() == WXK_ESCAPE && IsMultiLine() )
        return;
    event.Skip();	
#else
	//  Pass the event to the superclass handler
	wxTextCtrl::OnKeyDown(event);
#endif
}
