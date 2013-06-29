/*
 *  MyProgressIndicator.h
 *  Molby
 *
 *  Created by Toshi Nagata on 13/06/26.
 *  Copyright 2013 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __MyProgressIndicator_h__
#define __MyProgressIndicator_h__

#include "wx/window.h"

class MyProgressIndicator: public wxWindow
{
public:
	bool enabled;
	bool pressed;
	signed char indicatorState;
	
    MyProgressIndicator(wxWindow* parent, wxWindowID id, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = 0);
    void OnPaint(wxPaintEvent &event);
    void OnMouseEvent(wxMouseEvent &event);
	void OnCaptureLost(wxMouseCaptureLostEvent &event);
	void SetEnabled(bool flag);
	bool IsEnabled() { return enabled; }
	void SetPressed(bool flag);
	bool IsPressed() { return pressed; }
	void SetIndicatorState(int state);
	int IndicatorState() { return indicatorState; }
	void ProceedIndicatorState();

private:
    DECLARE_EVENT_TABLE()	
};

#endif  /*  MyProgressIndicator  */

