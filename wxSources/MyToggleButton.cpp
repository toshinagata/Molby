/*
 *  MyToggleButton
 *  Molby
 *
 *  Created by Toshi Nagata on 2022/09/17.
 *  Copyright 2022 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "MyToggleButton.h"
#include "wx/dcclient.h"
#include "wx/toplevel.h"

IMPLEMENT_DYNAMIC_CLASS(MyToggleButton, wxToggleButton)

BEGIN_EVENT_TABLE(MyToggleButton, wxToggleButton)
EVT_PAINT(MyToggleButton::OnPaint)
EVT_LEFT_DOWN(MyToggleButton::OnLeftDown)
EVT_LEFT_UP(MyToggleButton::OnLeftUp)
END_EVENT_TABLE()

void
MyToggleButton::OnPaint(wxPaintEvent &event)
{
    wxPaintDC dc(this);
    wxSize size = GetSize();
    dc.SetPen(*wxGREY_PEN);
    const wxBrush *brush;
    bool isActive;
    int col;
    //  Get the nearest TopLevelWindow and see if it is active or not
    wxWindow *win = this;
    while (win != NULL && !win->IsKindOf(wxCLASSINFO(wxTopLevelWindow)))
        win = win->GetParent();
    if (win != NULL && ((wxTopLevelWindow *)win)->IsActive())
        isActive = true;
    else isActive = false;
    if (IsPressed())
        col = (isActive ? 128 : 180);
    else if (GetValue())
        col = (isActive ? 180 : 220);
    else
        col = 240;
    brush = wxTheBrushList->FindOrCreateBrush(wxColour(col, col, col));
    dc.SetBrush(*brush);
    dc.DrawRectangle(0, 0, size.x, size.y);
    wxString label = GetLabel();
    int w, h, descent;
    int x, y;
    dc.GetTextExtent(label, &w, &h, &descent);
    x = (size.x - w) / 2;
    y = (size.y - h) / 2;
    if (isActive)
        dc.SetTextForeground(*wxBLACK);
    else
        dc.SetTextForeground(wxColour(128, 128, 128));
    dc.DrawText(label, x, y);
}

void
MyToggleButton::OnLeftDown(wxMouseEvent &event)
{
    CaptureMouse();
    SetPressed(true);
    Refresh();
}

void
MyToggleButton::OnLeftUp(wxMouseEvent &event)
{
    if (HasCapture()) {
        ReleaseMouse();
        SetPressed(false);
        int x = event.GetX();
        int y = event.GetY();
        wxSize sz = GetSize();
        if (x > 0 && x < sz.x && y > 0 && y < sz.y) {
            SetValue(!GetValue());
            wxCommandEvent cmdevt(wxEVT_TOGGLEBUTTON, GetId());
            cmdevt.SetInt(GetValue());
            cmdevt.SetEventObject(this);
            ProcessCommand(cmdevt);
            Refresh();
        }
    }
}
