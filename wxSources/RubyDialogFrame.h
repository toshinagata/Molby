/*
 *  RubyDialogFrame.h
 *  Molby
 *
 *  Created by Toshi Nagata on 08/12/05.
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

#ifndef __RubyDialogFrame_h__
#define __RubyDialogFrame_h__

#include "wx/dialog.h"
#include "wx/sizer.h"
#include "wx/panel.h"
#include "wx/timer.h"
#include "../Mollib/Ruby_bind/ruby_dialog.h"
#include "MyListCtrl.h"

/*  MyLayoutPanel: an empty subclass of wxPanel exclusively used in Dialog#layout  */
class MyLayoutPanel: public wxPanel {
public:
	MyLayoutPanel(): wxPanel() {}
	MyLayoutPanel(wxWindow *parent, wxWindowID winid = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize): wxPanel(parent, winid, pos, size) {}
	virtual ~MyLayoutPanel() {}		
private:
	DECLARE_DYNAMIC_CLASS(MyLayoutPanel)
};

/*  MyDrawingPanel: an empty subclass of wxPanel exclusively used as "view" dialog item  */
class MyDrawingPanel: public wxPanel {
public:
	MyDrawingPanel(): wxPanel() {}
	MyDrawingPanel(wxWindow *parent, wxWindowID winid = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize): wxPanel(parent, winid, pos, size) {}
	virtual ~MyDrawingPanel() {}		
private:
	DECLARE_DYNAMIC_CLASS(MyDrawingPanel)
};

class RubyDialogFrame: public wxDialog, public MyListCtrlDataSource {	
public:
	RDItem **ditems;
	int nditems;
	RubyValue dval;  /*  The Ruby value representing this object  */

	wxPanel *contentPanel;
	wxSizer *contentSizer;
	wxSizer *buttonSizer;
	wxBoxSizer *boxSizer;
	wxTimer *myTimer;
	
	wxWindow *currentDrawingItem;
	wxDC *currentContext;

	wxWindow *lastFocusedWindow;

	/*  Auto resizing  */
	RDSize mySize;  /*  Previous size  */
	bool autoResizeEnabled;  /*  true if auto resizing is enabled  */
	
	/*  Message bridge (with Ruby world); obj, event_type, event_id, proc val  */
	void **messageData;
	Int countMessageData;
	
	/*  On key handler (the handler is in the Ruby world)  */
	bool onKeyHandlerEnabled;

	RubyDialogFrame(wxWindow* parent, wxWindowID wid, const wxString& title, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
	virtual ~RubyDialogFrame();

	/*  For internal use (discard data in messageData[])  */
	void DiscardMessageData();
	
	void SetAutoResizeEnabled(int flag) { autoResizeEnabled = (flag != 0); }
	bool IsAutoResizeEnabled() { return autoResizeEnabled; }
	
	int ListenToObject(void *obj, const char *objtype, const char *msg, RubyValue oval, RubyValue pval);
	void HandleDocumentEvent(wxCommandEvent &event);
	void HandlePaintEvent(wxPaintEvent &event);

	int AddDialogItem(RDItem *item);
	RDItem *DialogItemAtIndex(int index);
	int SearchDialogItem(RDItem *item);
	void SetRubyObject(RubyValue val);
	void CreateStandardButtons(const char *oktitle, const char *canceltitle);
	int StartIntervalTimer(int millisec);
	void StopIntervalTimer(void);
	void OnDialogItemAction(wxCommandEvent &event);
	void OnTextUpdated(wxCommandEvent &event);
	void OnEnterProcessedOnText(wxCommandEvent &event);
	void OnKillFocusOnText(wxFocusEvent &event);

	void OnTimerEvent(wxTimerEvent &event);
	void OnDefaultButtonPressed(wxCommandEvent &event);
	void OnSize(wxSizeEvent &event);
	void OnChar(wxKeyEvent &event);
	void OnCloseWindow(wxCloseEvent &event);
	void OnActivate(wxActivateEvent &event);
	void OnChildFocus(wxChildFocusEvent &event);
	void OnCloseFromMenu(wxCommandEvent &event);

	//  MyListCtrlDataSource methods
	virtual int GetItemCount(MyListCtrl *ctrl);
	virtual wxString GetItemText(MyListCtrl *ctrl, long row, long column) const;
	virtual int SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value);
	virtual void DragSelectionToRow(MyListCtrl *ctrl, long row);
	virtual bool IsItemEditable(MyListCtrl *ctrl, long row, long column);
	virtual bool IsDragAndDropEnabled(MyListCtrl *ctrl);
	virtual void OnSelectionChanged(MyListCtrl *ctrl);
	virtual int SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg);
	virtual int HasPopUpMenu(MyListCtrl *ctrl, long row, long column, char ***menu_titles);
	virtual void OnPopUpMenuSelected(MyListCtrl *ctrl, long row, long column, int selected_index);
		
private:
	DECLARE_EVENT_TABLE()
};

#endif /* __RubyDialogFrame_h__ */
