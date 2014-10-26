/*
 *  MyGLCanvas.h
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

#ifndef __MyGLCanvas_h__
#define __MyGLCanvas_h__

#include "wx/glcanvas.h"
#include "../MolLib/MolLib.h"

class MoleculeView;

class MyGLCanvas: public wxGLCanvas
{
public:
    MoleculeView *view;
	wxGLContext *context;

    MyGLCanvas(MoleculeView *v, wxWindow *frame, const wxPoint& pos, const wxSize& size, long style = 0);
	~MyGLCanvas();
	void SetCurrent() { context->SetCurrent(*this); }
    void OnPaint(wxPaintEvent &event);
    void OnMouseEvent(wxMouseEvent &event);
    void OnEraseBackground(wxEraseEvent &event);
    void OnSize(wxSizeEvent &event);
	void OnChar(wxKeyEvent &event);
	void OnCaptureLost(wxMouseCaptureLostEvent &event);
	void Update();
	
private:
    DECLARE_EVENT_TABLE()
};


#endif  /* __MyGLCanvas_h__ */

