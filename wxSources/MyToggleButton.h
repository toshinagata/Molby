/*
 *  MyToggleButton.h
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

#ifndef __MyToggleButton_h__
#define __MyToggleButton_h__

#include "wx/tglbtn.h"

class MyToggleButton: public wxToggleButton
{
public:
    MyToggleButton() {}
    MyToggleButton (wxWindow *parent, wxWindowID id, const wxString &label, const wxPoint &pos=wxDefaultPosition, const wxSize &size=wxDefaultSize, long style=0, const wxValidator &val=wxDefaultValidator, const wxString &name=wxCheckBoxNameStr)
    : wxToggleButton(parent, id, label, pos, size, style, val, name) {
        m_isPressed = m_value = false;
    }

    void OnPaint(wxPaintEvent &event);
    void OnLeftDown(wxMouseEvent &event);
    void OnLeftUp(wxMouseEvent &event);
    void SetValue(bool val) { m_value = val; Refresh(); }
    bool GetValue() { return m_value; }
    void SetPressed(bool flag) { m_isPressed = flag; Refresh(); }
    bool IsPressed() { return m_isPressed; }
    
private:
    bool m_isPressed;
    bool m_value;
    bool m_isEnabled;
    DECLARE_DYNAMIC_CLASS(MyToggleButton)
    DECLARE_EVENT_TABLE()
};


#endif  /* __MyToggleButton_h__ */

