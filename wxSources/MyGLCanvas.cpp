/*
 *  MyGLCanvas.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 08/10/24.
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

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !wxUSE_DOC_VIEW_ARCHITECTURE
#error "You should have DocView architecture enabled in your wxWidgets installation."
#endif

#include "wx/docview.h"

#include "MyGLCanvas.h"
#include "MoleculeView.h"
#include "MyApp.h"

#include "../MolLib/MolLib.h"

BEGIN_EVENT_TABLE(MyGLCanvas, wxGLCanvas)
    EVT_SIZE(MyGLCanvas::OnSize)
    EVT_ERASE_BACKGROUND(MyGLCanvas::OnEraseBackground)
    EVT_MOUSE_EVENTS(MyGLCanvas::OnMouseEvent)
    EVT_PAINT(MyGLCanvas::OnPaint)
	EVT_MOUSE_CAPTURE_LOST(MyGLCanvas::OnCaptureLost)
END_EVENT_TABLE()

#ifdef __WXMSW__
int *MyGLAttributes = NULL;
#else
int MyGLAttributes[20] = { WX_GL_RGBA, WX_GL_MIN_RED, 1, WX_GL_MIN_GREEN, 1,
        WX_GL_MIN_BLUE, 1, WX_GL_DEPTH_SIZE, 1,
        WX_GL_DOUBLEBUFFER,
#  if defined(__WXMAC__) || defined(__WXCOCOA__)
        GL_NONE
#  else
        None
#  endif
	};
#endif

// Define a constructor for my canvas
MyGLCanvas::MyGLCanvas(MoleculeView *v, wxWindow *frame, const wxPoint& pos, const wxSize& size, long style):
  wxGLCanvas(frame, wxID_ANY, MyGLAttributes, pos, size, style)
{
	view = v;
	context = new wxGLContext(this);
}

MyGLCanvas::~MyGLCanvas()
{
	delete context;
}

// Define the repainting behaviour
void
MyGLCanvas::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
    wxPaintDC dc(this);

#ifndef __WXMOTIF__
//    if (!GetContext()) return;
#endif

    context->SetCurrent(*this);
    if (view)
      view->OnDraw(&dc);
	else {
		glClearColor (0, 0, 0, 0);
		glClear(GL_COLOR_BUFFER_BIT |
				GL_DEPTH_BUFFER_BIT);
	}
    SwapBuffers();
}

void
MyGLCanvas::OnCaptureLost(wxMouseCaptureLostEvent &event)
{
	wxPoint pt;
	int w, h, modifierFlags;
	float p[2];
	MainView *mview;

	if (view == NULL)
		return;

	pt = ScreenToClient(::wxGetMousePosition());
	GetClientSize(&w, &h);
	p[0] = pt.x;
	p[1] = h - pt.y;
	modifierFlags = MainViewCallback_modifierFlags(NULL);
	mview = ((MoleculeView *)view)->mview;
	MoleculeLock(mview->mol);
	MainView_mouseUp(((MoleculeView *)view)->mview, p, modifierFlags, 0);
	MoleculeUnlock(mview->mol);
}

void
MyGLCanvas::OnMouseEvent(wxMouseEvent &event)
{
	if (!view)
		return;

	wxPoint pt(event.GetPosition());
	float p[2];
	int modifierFlags, clickCount;
	MainView *mview = ((MoleculeView *)view)->mview;
	int w, h;
	GetClientSize(&w, &h);

	p[0] = pt.x;
	p[1] = h - pt.y;  /*  The origin of the internal coordinate should be left-bottom  */
	modifierFlags = MainViewCallback_modifierFlags(&event);
	clickCount = MainViewCallback_clickCount(&event);

	if (event.LeftUp() /* || (mview->isDragging && event.Entering() && !event.LeftIsDown()) */) {
		if (HasCapture())
			ReleaseMouse();
		MoleculeLock(mview->mol);
		MainView_mouseUp(mview, p, modifierFlags, clickCount);
		MoleculeUnlock(mview->mol);
	} else if (event.LeftDClick()) {
		if (HasCapture())
			ReleaseMouse();
		MoleculeLock(mview->mol);
		MainView_mouseUp(mview, p, modifierFlags, 2);
		MoleculeUnlock(mview->mol);
	} else if (event.Dragging()) {
		MoleculeLock(mview->mol);
		MainView_mouseDragged(mview, p, modifierFlags);
		MoleculeUnlock(mview->mol);
	} else if (event.LeftDown()) {
		if (wxWindow::FindFocus() != this)
			SetFocus();
		CaptureMouse();
		MoleculeLock(mview->mol);
		MainView_mouseDown(mview, p, modifierFlags);
		MoleculeUnlock(mview->mol);
	} else event.Skip();
}

void
MyGLCanvas::OnSize(wxSizeEvent &event)
{
    // this is also necessary to update the context on some platforms
//    wxGLCanvas::OnSize(event);

    // set GL viewport (not called by wxGLCanvas::OnSize on all platforms...)
    int w, h;
    GetClientSize(&w, &h);
#ifndef __WXMOTIF__
//    if (GetContext())
#endif
    {
        context->SetCurrent(*this);
        glViewport(0, 0, (GLint) w, (GLint) h);
    }
#if defined(__WXMSW__)
	//  On MSW, the window is not repainted upon resize event
	Refresh();
#endif
}

void
MyGLCanvas::OnEraseBackground(wxEraseEvent & WXUNUSED(event))
{
    // Do nothing, to avoid flashing.
}
