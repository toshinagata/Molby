/*
 *  MyTextCtrl.h
 *
 *  Created by Toshi Nagata on 2014/09/21.
 *  Copyright 2014 Toshi Nagata. All rights reserved.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 */

#ifndef __MyTextCtrl_h__
#define __MyTextCtrl_h__

#include "wx/textctrl.h"

#define MyTextCtrl_Process_Escape    0x10000000

extern const wxEventType myTextCtrl_EVT_PROCESS_ESCAPE;

class MyTextCtrl: public wxTextCtrl {

public:
	MyTextCtrl();
	MyTextCtrl(wxWindow *parent, wxWindowID id, const wxString &value=wxEmptyString,
			   const wxPoint &pos=wxDefaultPosition, const wxSize &size=wxDefaultSize,
			   long style=0, const wxValidator &validator=wxDefaultValidator,
			   const wxString &name=wxTextCtrlNameStr);
	virtual ~MyTextCtrl();
	bool processEscape;
	void OnKeyUp(wxKeyEvent& event);

private:
	DECLARE_DYNAMIC_CLASS(MyTextCtrl)
	DECLARE_EVENT_TABLE()
};

#endif /* __MyListCtrl_h__ */
