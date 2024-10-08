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
#include "MyToggleButton.h"

class MyDocument;
class MyGLCanvas;
class wxMenu;
class wxStaticText;
class wxChoice;
class MyProgressIndicator;

class MoleculeView: public wxView, public MyListCtrlDataSource
{
public:
    wxDocChildFrame *frame;
    MyGLCanvas *canvas;
    MainView *mview;
	wxChoice *listmenu;
	MyListCtrl *listctrl;
	wxMenu *file_history_menu;
	wxMenu *edit_menu;
	MyToggleButton *tbuttons[6];
  wxBitmapButton *bbuttons[4];
	wxStaticText *infotext;
	MyProgressIndicator *progress;
	wxPanel *frameControlPanel;
	wxSlider *frameSlider;
	wxTextCtrl *frameText;

	bool isRebuildingTable;

    MoleculeView() { canvas = (MyGLCanvas *) NULL; frame = (wxDocChildFrame *) NULL; }
    virtual ~MoleculeView();

    MyDocument *MolDocument() { return (MyDocument *)m_viewDocument; }
	MyListCtrl *GetListCtrl() { return listctrl; }
	MyToggleButton *GetToggleButtonAtIndex(int i) { return (i >= 0 && i < 6 ? tbuttons[i] : NULL); }

    bool OnCreate(wxDocument *doc, long flags);
    void OnDraw(wxDC *dc);
    void OnUpdate(wxView *sender, wxObject *hint = (wxObject *) NULL);
    bool OnClose(bool deleteWindow = true);
	wxImage *CaptureGLCanvas(float scale = 1.0, int bg_color = -1, int width = 0, int height = 0);
	int	 DoExportGraphic(wxString& fname, float scale, int bg_color, int width, int height);

	virtual void Activate (bool activate);
	virtual wxPrintout *OnCreatePrintout();

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

	void InvalidateProgressIndicator();
	void ProceedProgressIndicator();
	void OnStopProgressPressed(wxCommandEvent& event);

	void SelectTable(int idx);
	void OnSelectTable(wxCommandEvent &event);

	void OnActivate(wxActivateEvent &event);
	
	void OnMoleculeReplaced();  /*  Called when Molecule is replaced within MyDocument  */

	void EnableProgressIndicator(bool flag);

	/*  MyListCtrlDataSource functions  */
	virtual int GetItemCount(MyListCtrl *ctrl);
	virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const;
	virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value);
	virtual void DragSelectionToRow(MyListCtrl *ctrl, long row);
	virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column);
    virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl, long row);
    virtual void OnSelectionChanged(MyListCtrl *ctrl);
	virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg);
    virtual bool IsRowSelectable(MyListCtrl *ctrl, long row);

private:
    DECLARE_DYNAMIC_CLASS(MoleculeView)
    DECLARE_EVENT_TABLE()
};

#endif
