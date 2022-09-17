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
    if (IsPressed())
        brush = wxTheBrushList->FindOrCreateBrush(wxColour(128, 128, 128));
    else if (GetValue())
        brush = wxTheBrushList->FindOrCreateBrush(wxColour(180, 180, 180));
    else
        brush = wxTheBrushList->FindOrCreateBrush(wxColour(240, 240, 240));
    dc.SetBrush(*brush);
    dc.DrawRectangle(0, 0, size.x, size.y);
    wxString label = GetLabel();
    int w, h, descent;
    int x, y;
    dc.GetTextExtent(label, &w, &h, &descent);
    x = (size.x - w) / 2;
    y = (size.y - h) / 2;
    dc.SetPen(*wxBLACK_PEN);
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
