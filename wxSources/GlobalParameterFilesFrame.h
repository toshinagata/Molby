/*
 *  GlobalParameterFilesFrame.h
 *  Molby
 *
 *  Created by Toshi Nagata on 09/11/14.
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

#ifndef __GlobalParameterFilesFrame_h__
#define __GlobalParameterFilesFrame_h__

#include "wx/mdi.h"
#include "MyListCtrl.h"
#include "wx/sizer.h"

class wxMenu;
class GlobalParameterFilesFrame: public wxFrame, public MyListCtrlDataSource
{
public:
	MyListCtrl *listctrl;
	wxButton *add_button, *remove_button;
	wxMenu *file_history_menu;
	wxMenu *edit_menu;
	
	GlobalParameterFilesFrame(wxWindow *parent, const wxString& title, const wxPoint& pos, const wxSize& size, long type);
	virtual ~GlobalParameterFilesFrame();
	
	MyListCtrl *GetListCtrl() { return listctrl; }
	
	void OnCreate();
	
	static GlobalParameterFilesFrame *CreateGlobalParameterFilesFrame(wxWindow *parent);
	
	void OnCloseWindow(wxCloseEvent &event);
	void OnClose(wxCommandEvent &event);
	void OnUpdateUI(wxUpdateUIEvent& event);
	void OnAddGlobalParameterFile(wxCommandEvent& event);
	void OnRemoveGlobalParameterFile(wxCommandEvent& event);

	/*  MyListCtrlDataSource functions  */
	virtual int GetItemCount(MyListCtrl *ctrl);
	virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const;
	virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value);
	virtual void DragSelectionToRow(MyListCtrl *ctrl, long row);
	virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column);
	virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl, long row);
	virtual void OnSelectionChanged(MyListCtrl *ctrl);
	virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg);
	
private:
	DECLARE_EVENT_TABLE()
};

#endif /* __GlobalParameterFilesFrame_h__ */
