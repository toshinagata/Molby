/*
 *  GlobalParameterFrame.h
 *  Molby
 *
 *  Created by Toshi Nagata on 09/11/05.
 *  Copyright 2009 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __GlobalParameterFrame_h__
#define __GlobalParameterFrame_h__

#include "wx/mdi.h"
#include "MyListCtrl.h"

class wxMenu;

class GlobalParameterFrame: public wxFrame, public MyListCtrlDataSource
{

public:
	MyListCtrl *listCtrl;
	wxMenu *file_history_menu;
	wxMenu *edit_menu;

	GlobalParameterFrame(wxWindow *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type);
	virtual ~GlobalParameterFrame();

	MyListCtrl *GetListCtrl() { return listCtrl; }

	void OnCreate();

	static GlobalParameterFrame *CreateGlobalParameterFrame(wxWindow *parent);
	
	void OnCloseWindow(wxCloseEvent &event);
	void OnClose(wxCommandEvent &event);
	void OnUpdateUI(wxUpdateUIEvent& event);

	/*  MyListCtrlDataSource functions  */
	virtual int GetItemCount(MyListCtrl *ctrl);
	virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const;
	virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value);
	virtual void DragSelectionToRow(MyListCtrl *ctrl, long row);
	virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column);
	virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl);
	virtual void OnSelectionChanged(MyListCtrl *ctrl);
	virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg);

private:
	DECLARE_EVENT_TABLE()
};

#endif /* __GlobalParameterFrame_h__ */
