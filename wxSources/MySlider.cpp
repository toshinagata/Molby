/*
 *  MySlider.cpp
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

#include "MySlider.h"
#include "MyApp.h"

BEGIN_EVENT_TABLE(MySlider, wxToggleButton)
EVT_MOUSE_EVENTS(MySlider::OnMouseEvent)
// EVT_PAINT(MySlider::OnPaint)
EVT_MOUSE_CAPTURE_LOST(MySlider::OnCaptureLost)
END_EVENT_TABLE()

const wxEventType MySliderEvent = wxNewEventType();

MySlider::MySlider(wxWindow* parent, wxWindowID id, int direction, const wxPoint& pos, const wxSize& size, long style, const wxValidator& validator, const wxString& name):
	wxToggleButton(parent, id, wxEmptyString, pos, size, style, validator, name)
{
	m_direction = direction;
	SetCursor(wxCursor(wxCURSOR_HAND));
	mouseStatus = 0;
}

//void
//MySlider::OnPaint(wxPaintEvent &event)
//{
//	wxToggleButton::OnPaint(event);
//}

void
MySlider::OnCaptureLost(wxMouseCaptureLostEvent &event)
{
	mouseStatus = 0;
}

void
MySlider::OnMouseEvent(wxMouseEvent &event)
{
	wxPoint pt(event.GetPosition());
	bool sendAction = false;
	
	if (mouseStatus == 0) {
		if (event.LeftDown()) {
			mouseDownPoint = pt;
			mouseDragPoint = pt;
			mouseStatus = 1;
			sendAction = true;
			CaptureMouse();
		}
	} else {
		if (event.Dragging()) {
			mouseDragPoint = pt;
			mouseStatus = 2;
			sendAction = true;
		} else if (event.LeftUp() /* || (event.Entering() && !event.LeftIsDown()) */) {
			mouseStatus = 0;
			sendAction = true;
			ReleaseMouse();
	//	} else if (event.Leaving()) {
	//		wxGetApp().NotifyMouseUpEvent(this);
		}
	}
	if (sendAction) {
		/*  Send a custom event  */
		wxCommandEvent event(MySliderEvent, GetId());
		event.SetEventObject(this);
		GetEventHandler()->ProcessEvent(event);
	}
}

float
MySlider::GetFloatValue()
{
	int width, height;
	float f;
	GetSize(&width, &height);
	if (m_direction == wxHORIZONTAL) {
		f = (mouseDragPoint.x - mouseDownPoint.x) / (float)width;
	} else {
		/*  Vertical  */
		/*  Note the flipped coordinate! (up is positive internally)  */
		f = -(mouseDragPoint.y - mouseDownPoint.y) / (float)height;
	}
	return f;
}
