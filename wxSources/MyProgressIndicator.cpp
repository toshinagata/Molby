/*
 *  MyProgressIndicator.cpp
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

#include "MyProgressIndicator.h"
#include "MyApp.h"

#include "wx/dcclient.h"
#include "wx/bitmap.h"
#include "wx/image.h"   //  If omitted, then new wxBitmap(char *, int) will fail

BEGIN_EVENT_TABLE(MyProgressIndicator, wxWindow)
EVT_MOUSE_EVENTS(MyProgressIndicator::OnMouseEvent)
EVT_PAINT(MyProgressIndicator::OnPaint)
EVT_MOUSE_CAPTURE_LOST(MyProgressIndicator::OnCaptureLost)
END_EVENT_TABLE()

static wxBitmap *sStopMiniIcons[3];
static wxBitmap *sProgressIndicatorIcons[12];

MyProgressIndicator::MyProgressIndicator(wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style) : wxWindow(parent, id, pos, size, style)
{
	enabled = false;
	pressed = false;
	indicatorState = 0;

	// Initialize icons for progress indicator
	if (sStopMiniIcons[0] == NULL) {
#include "../bitmaps/stop_mini_grey.xpm"
#include "../bitmaps/stop_mini.xpm"
#include "../bitmaps/stop_mini_dark.xpm"
#include "../bitmaps/pi00.xpm"
#include "../bitmaps/pi01.xpm"
#include "../bitmaps/pi02.xpm"
#include "../bitmaps/pi03.xpm"
#include "../bitmaps/pi04.xpm"
#include "../bitmaps/pi05.xpm"
#include "../bitmaps/pi06.xpm"
#include "../bitmaps/pi07.xpm"
#include "../bitmaps/pi08.xpm"
#include "../bitmaps/pi09.xpm"
#include "../bitmaps/pi10.xpm"
#include "../bitmaps/pi11.xpm"
		sStopMiniIcons[0] = new wxBitmap(stop_mini_grey, wxBITMAP_TYPE_XPM);
		sStopMiniIcons[1] = new wxBitmap(stop_mini, wxBITMAP_TYPE_XPM);
		sStopMiniIcons[2] = new wxBitmap(stop_mini_dark, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[0] = new wxBitmap(pi00, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[1] = new wxBitmap(pi01, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[2] = new wxBitmap(pi02, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[3] = new wxBitmap(pi03, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[4] = new wxBitmap(pi04, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[5] = new wxBitmap(pi05, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[6] = new wxBitmap(pi06, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[7] = new wxBitmap(pi07, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[8] = new wxBitmap(pi08, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[9] = new wxBitmap(pi09, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[10] = new wxBitmap(pi10, wxBITMAP_TYPE_XPM);
		sProgressIndicatorIcons[11] = new wxBitmap(pi11, wxBITMAP_TYPE_XPM);
	}
}

void
MyProgressIndicator::OnPaint(wxPaintEvent &event)
{
	int stopState;
	wxPaintDC dc(this);
	
	dc.Clear();
	if (enabled)
		stopState = (pressed ? 2 : 1);
	else stopState = 0;
	dc.DrawBitmap(*sStopMiniIcons[stopState], 0, 0, 1);
	if (enabled)
		dc.DrawBitmap(*sProgressIndicatorIcons[indicatorState % 12], 0, 15, 1);
}

void
MyProgressIndicator::OnMouseEvent(wxMouseEvent &event)
{
	long x, y;
	if (enabled) {
		event.GetPosition(&x, &y);
		if (event.LeftDown()) {
			if (x >= 0 && x <= 12 && y >= 0 && y <= 12) {
				pressed = 1;
				CaptureMouse();
			}
		} else if (event.Dragging()) {
			if (x >= 0 && x <= 12 && y >= 0 && y <= 12) {
				pressed = 1;
			} else pressed = 0;
		} else if (event.LeftUp()) {
			if (pressed) {
				/*  Send button action  */
				wxCommandEvent event(wxEVT_COMMAND_BUTTON_CLICKED, GetId());
				event.SetEventObject(this);
				GetEventHandler()->ProcessEvent(event);
			}
			if (HasCapture())
				ReleaseMouse();
		}
	}
	if (event.LeftDown())
		event.Skip();
}

void
MyProgressIndicator::OnCaptureLost(wxMouseCaptureLostEvent &event)
{
}

void
MyProgressIndicator::SetEnabled(bool flag)
{
	if (flag != enabled) {
		enabled = flag;
		Refresh();
	}
}

void
MyProgressIndicator::SetPressed(bool flag)
{
	if (flag != pressed) {
		pressed = flag;
		Refresh();
	}
}

void
MyProgressIndicator::SetIndicatorState(int state)
{
	if (state != indicatorState) {
		indicatorState = state;
		Refresh();
	}
}

void
MyProgressIndicator::ProceedIndicatorState()
{
	indicatorState = (indicatorState + 1) % 12;
	Refresh();
}


