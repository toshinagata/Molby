/*
 *  MySlider.h
 *  Molby
 *
 *  Created by Toshi Nagata on 08/11/15.
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

#ifndef __MySlider_h__
#define __MySlider_h__

#include "wx/tglbtn.h"
#include "wx/string.h"
#include "wx/slider.h"

extern const wxEventType MySliderEvent;

class MySlider: public wxToggleButton
{
public:
	int m_direction;
	int mouseStatus;  /*  0: mouseUp, 1: mouseDown, 2: dragging  */
	wxPoint mouseDownPoint;
	wxPoint mouseDragPoint;
	
    MySlider(wxWindow* parent, wxWindowID id, int direction, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = 0, const wxValidator& validator = wxDefaultValidator, const wxString& name = wxT("slider"));
    void OnPaint(wxPaintEvent &event);
    void OnMouseEvent(wxMouseEvent &event);
	int GetMouseStatus() { return mouseStatus; }
	float GetFloatValue();
	void OnCaptureLost(wxMouseCaptureLostEvent &event);

private:
    DECLARE_EVENT_TABLE()	
};

#endif /* __MySlider_h__ */
