/*
 *  MoleculeView.h
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

#ifndef __MoleculeView_h__
#define __MoleculeView_h__

#include "wx/docview.h"

#include "../MolLib/MolLib.h"
#include "MyListCtrl.h"

class MyDocument;
class MyGLCanvas;
class wxMenu;
class wxToggleButton;
class wxStaticText;
class wxChoice;

class MoleculeView: public wxView, public MyListCtrlDataSource
{
public:
    wxMDIChildFrame *frame;
    MyGLCanvas *canvas;
    MainView *mview;
	wxChoice *listmenu;
	MyListCtrl *listctrl;
	wxMenu *file_history_menu;
	wxMenu *edit_menu;
	wxToggleButton *tbuttons[6];
	wxStaticText *infotext;
	wxPanel *frameControlPanel;
	wxSlider *frameSlider;
	wxTextCtrl *frameText;

	bool isRebuildingTable;

    MoleculeView() { canvas = (MyGLCanvas *) NULL; frame = (wxMDIChildFrame *) NULL; }
    ~MoleculeView() {}

    MyDocument *MolDocument() { return (MyDocument *)m_viewDocument; }
	MyListCtrl *GetListCtrl() { return listctrl; }

    bool OnCreate(wxDocument *doc, long flags);
    void OnDraw(wxDC *dc);
    void OnUpdate(wxView *sender, wxObject *hint = (wxObject *) NULL);
    bool OnClose(bool deleteWindow = true);

	void OnButtonPressed(wxCommandEvent &event);
	void OnSliderAction(wxCommandEvent &event);

	void OnFrameButtonAction(wxMouseEvent &event);
	void OnFrameSliderAction(wxScrollEvent &event);
	void OnFrameTextAction(wxCommandEvent &event);

	void OnDocumentModified(wxCommandEvent &event);
	void OnChar(wxKeyEvent &event);
	void OnScriptMenuModified(wxCommandEvent& event);
	void OnLeftDClickInListCtrl(wxMouseEvent& event);

	void SelectButtonForMode(int mode);
	void UpdateFrameControlValues();
	void UpdateFrameControls();
	
	void SelectTable(int idx);
	void OnSelectTable(wxCommandEvent &event);

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
    DECLARE_DYNAMIC_CLASS(MoleculeView)
    DECLARE_EVENT_TABLE()
};

#endif
