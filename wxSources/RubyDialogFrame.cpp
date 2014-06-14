/*
 *  RubyDialogFrame.cpp
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

#include "wx/gdicmn.h"
#include "wx/stattext.h"
#include "wx/textctrl.h"
#include "wx/button.h"
#include "wx/tglbtn.h"
#include "wx/filedlg.h"
#include "wx/dirdlg.h"
#include "wx/dcclient.h"
#include "wx/choice.h"
#include "wx/checkbox.h"
#include "wx/radiobut.h"
#include "wx/statline.h"
#include "wx/settings.h"
#include "wx/dcmemory.h"

#include "RubyDialogFrame.h"
#include "MyApp.h"
#include "MyMBConv.h"
#include "MyDocument.h"

IMPLEMENT_DYNAMIC_CLASS(MyLayoutPanel, wxPanel)

IMPLEMENT_DYNAMIC_CLASS(MyDrawingPanel, wxPanel)

BEGIN_EVENT_TABLE(RubyDialogFrame, wxDialog)
  EVT_TIMER(-1, RubyDialogFrame::OnTimerEvent)
  EVT_BUTTON(wxID_OK, RubyDialogFrame::OnDefaultButtonPressed)
  EVT_BUTTON(wxID_CANCEL, RubyDialogFrame::OnDefaultButtonPressed)
  EVT_SIZE(RubyDialogFrame::OnSize)
  EVT_CHAR(RubyDialogFrame::OnChar)
  EVT_CLOSE(RubyDialogFrame::OnCloseWindow)
END_EVENT_TABLE()

RubyDialogFrame::RubyDialogFrame(wxWindow* parent, wxWindowID wid, const wxString& title, const wxPoint& pos, const wxSize& size, long style):
	wxDialog(parent, wid, title, pos, size, style)
{
	ditems = NULL;
	nditems = 0;
	dval = NULL;
	mySize = gZeroSize;
	autoResizeEnabled = true;
	messageData = NULL;
	countMessageData = 0;
	onKeyHandlerEnabled = false;
	currentContext = NULL;
	currentDrawingItem = NULL;
	
	//  Create a vertical box sizer that contains a panel containing all controls and a sizer containing
	//  OK/Cancel buttons
	contentSizer = new wxBoxSizer(wxVERTICAL);
	contentPanel = NULL;
	buttonSizer = NULL;  //  Will be created later
	myTimer = NULL;  //  Will be created when necessary
	boxSizer = new wxBoxSizer(wxVERTICAL);
	boxSizer->Add(contentSizer, 1, wxALL | wxEXPAND, 14);
	this->SetSizer(boxSizer);
	boxSizer->Layout();
	this->CentreOnScreen();
}

RubyDialogFrame::~RubyDialogFrame()
{
	if (myTimer != NULL)
		delete myTimer;
	if (ditems != NULL)
		free(ditems);
	DiscardMessageData();
}

void
RubyDialogFrame::DiscardMessageData()
{
	int i;
	for (i = 0; i < countMessageData; i++) {
		if (messageData[i * 5] != NULL) {
			wxEventType eventType = (wxEventType)messageData[i * 5 + 2];
			wxEvtHandler *handler = NULL;
			if (eventType == MyDocumentEvent)
				handler = MyDocumentFromMolecule((Molecule *)messageData[i * 5]);
			if (handler != NULL) {
				handler->Disconnect((int)messageData[i * 5 + 1], eventType, wxCommandEventHandler(RubyDialogFrame::HandleDocumentEvent), NULL, this);
			}
		}
	}
	free(messageData);
	messageData = NULL;
	countMessageData = 0;
}

int
RubyDialogFrame::AddDialogItem(RDItem *item)
{
	if (nditems % 8 == 0) {
		if (nditems == 0)
			ditems = (RDItem **)malloc(sizeof(RDItem *) * 8);
		else
			ditems = (RDItem **)realloc(ditems, sizeof(RDItem *) * (nditems + 8));
	}
	ditems[nditems++] = item;
	if (item != NULL && ((wxWindow *)item)->IsKindOf(CLASSINFO(MyLayoutPanel))) {
		wxSize size = ((MyLayoutPanel *)item)->GetSize();
		if (contentPanel == NULL)
			contentSizer->Add((MyLayoutPanel *)item, 1, wxEXPAND);
		else
			contentSizer->Replace(contentPanel, (MyLayoutPanel *)item);
		contentSizer->SetItemMinSize((MyLayoutPanel *)item, size.GetWidth(), size.GetHeight());
		contentPanel = (MyLayoutPanel *)item;
		boxSizer->Layout();
		Fit();
	}
	return nditems - 1;
}

RDItem *
RubyDialogFrame::DialogItemAtIndex(int index)
{
	if (index >= 0 && index < nditems)
		return ditems[index];
	else return NULL;
}

int
RubyDialogFrame::SearchDialogItem(RDItem *item)
{
	int i;
	for (i = 0; i < nditems; i++) {
		if (item == ditems[i])
			return i;
	}
	return -1;
}

void
RubyDialogFrame::SetRubyObject(RubyValue val)
{
	dval = val;
	if (dval == NULL) {
		/*  Stop message mechanism (because this object is already disconnected from the Ruby world)  */
		DiscardMessageData();
	}
}

/*  Create standard buttons. If oktitle/canceltitle == NULL, then the button is created but set hidden.
    If the title is "", then the default titles are used.  */
void
RubyDialogFrame::CreateStandardButtons(const char *oktitle, const char *canceltitle)
{
	wxSizer *sizer = CreateButtonSizer(wxOK | wxCANCEL);
	if (oktitle != NULL || canceltitle != NULL) {
		if (sizer == NULL)
			return;  /*  Cannot create  */
		if (buttonSizer == NULL) {
			boxSizer->Add(sizer, 0, wxBOTTOM | wxLEFT | wxRIGHT | wxEXPAND, 14);
			buttonSizer = sizer;
		} else {
			boxSizer->Replace(buttonSizer, sizer);
			buttonSizer = sizer;
		}
	} else {
		/*  Buttons are created but sizer is left unregistered  */
		if (buttonSizer != NULL)
			boxSizer->Remove(buttonSizer);
		buttonSizer = NULL;
	}
	while (nditems < 2)
		AddDialogItem(NULL);
	if (ditems[0] != NULL)
		((wxWindow *)ditems[0])->Destroy();
	if (ditems[1] != NULL)
		((wxWindow *)ditems[1])->Destroy();
	ditems[0] = (RDItem *)wxWindow::FindWindowById(wxID_OK, this);
	ditems[1] = (RDItem *)wxWindow::FindWindowById(wxID_CANCEL, this);
	if (oktitle == NULL) {
		((wxWindow *)ditems[0])->Show(false);
		((wxWindow *)ditems[0])->Enable(false);
	} else {
		if (oktitle[0] != 0) {
			wxString label1(oktitle, WX_DEFAULT_CONV);
			((wxButton *)ditems[0])->SetLabel(label1);
		}
		((wxWindow *)ditems[0])->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, this);
	}
	if (canceltitle == NULL) {
		((wxWindow *)ditems[1])->Show(false);
		((wxWindow *)ditems[1])->Enable(false);
	} else {
		if (canceltitle[0] != 0) {
			wxString label2(canceltitle, WX_DEFAULT_CONV);
			((wxButton *)ditems[1])->SetLabel(label2);
		}
		((wxWindow *)ditems[1])->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, this);
	}
}

int
RubyDialogFrame::StartIntervalTimer(int millisec)
{
	if (myTimer == NULL) {
		myTimer = new wxTimer(this);
	}
	return myTimer->Start(millisec);
}

void
RubyDialogFrame::StopIntervalTimer(void)
{
	if (myTimer != NULL)
		myTimer->Stop();
}

void
RubyDialogFrame::OnDialogItemAction(wxCommandEvent &event)
{
	RubyDialog_doItemAction((RubyValue)dval, (RDItem *)(event.GetEventObject()));
}

void
RubyDialogFrame::OnDefaultButtonPressed(wxCommandEvent &event)
{
	/*  Ignore the wxID_OK and wxID_CANCEL requests if the default buttons are hidden  */
	wxWindow *item = wxWindow::FindWindowById(event.GetId(), this);
	if (!item->IsShown())
		return;
	event.Skip();
}

void
RubyDialogFrame::OnTimerEvent(wxTimerEvent &event)
{
	RubyDialog_doTimerAction((RubyValue)dval);
}

static void
sResizeSubWindows(RubyValue dval, wxWindow *win, int dx, int dy)
{
	wxWindowList & children = win->GetChildren();
	wxWindowList::Node *node;
	for (node = children.GetFirst(); node; node = node->GetNext()) {
		int i, d, f, d1, d2, d3, ddx, ddy;
		wxWindow *current = (wxWindow *)node->GetData();
		wxRect frame = current->GetRect();
		int flex = RubyDialog_getFlexFlags(dval, (RDItem *)current);
		if (flex < 0)
			continue;		
		for (i = 0, f = flex; i < 2; i++, f /= 2) {
			if (i == 0)
				d = dx;
			else
				d = dy;
			switch (f & 21) {  /*  left, right, width (or top, bottom, height) */
				case 21:  /*  all flex  */
					d1 = d2 = d / 3;
					d3 = d - d1 - d2;
					break;
				case 5:   /*  left & right  */
					d1 = d / 2;
					d2 = 0;
					d3 = d - d1;
					break;
				case 17:  /*  left & width  */
					d1 = d / 2;
					d2 = d - d1;
					d3 = 0;
					break;
				case 20:  /*  right & width  */
					d1 = 0;
					d2 = d / 2;
					d3 = d - d2;
					break;
				case 1:   /*  left  */
					d1 = d;
					d2 = d3 = 0;
					break;
				case 4:   /*  right */
					d3 = d;
					d1 = d2 = 0;
					break;
				case 16:  /*  width  */
					d2 = d;
					d1 = d3 = 0;
					break;
				default:  /*  no resize  */
					d1 = d2 = d3 = 0;
					break;
			}
			if (i == 0) {
				frame.x += d1;
				frame.width += d2;
				ddx = d2;
			} else {
				frame.y += d1;
				frame.height += d2;
				ddy = d2;
			}
		}
		if (ddx != 0 || ddy != 0)
			sResizeSubWindows(dval, current, ddx, ddy);
		current->SetSize(frame);
#if defined(__WXMSW__)
		if (current->IsKindOf(CLASSINFO(MyDrawingPanel))) {
			//  Cause repaint
			current->Refresh();
		}
#endif
	}
}

void
RubyDialogFrame::OnSize(wxSizeEvent &event)
{
	wxSize size = GetClientSize();
	if (mySize.width != 0 && mySize.height != 0 && /*(mySize.width != size.x || mySize.height != size.y) &&*/ autoResizeEnabled) {
		/*  Resize the subviews  */
		sResizeSubWindows((RubyValue)dval, this, size.x - mySize.width, size.y - mySize.height);
	}
	mySize.width = size.x;
	mySize.height = size.y;
	event.Skip();
}

void
RubyDialogFrame::OnChar(wxKeyEvent &event)
{

	int code = event.GetKeyCode();
	if (onKeyHandlerEnabled) {
		RubyDialog_doKeyAction((RubyValue)dval, code);
	} else
		event.Skip();
}

void
RubyDialogFrame::OnCloseWindow(wxCloseEvent &event)
{
	RubyDialog_doCloseWindow((RubyValue)dval, IsModal());

	//  Check if all windows are gone
	wxGetApp().CheckIfAllWindowsAreGone(NULL);
}

static wxEvtHandler *
sGetEventHandlerFromObjectAndType(void *obj, wxEventType eventType)
{
	if (eventType == MyDocumentEvent)
		return MyDocumentFromMolecule((Molecule *)obj);
	else return NULL;
}

int
RubyDialogFrame::ListenToObject(void *obj, const char *objtype, const char *msg, RubyValue oval, RubyValue pval)
{
	int i, j, eventId = -1;
	wxEventType eventType = -1;
	wxEvtHandler *handler = NULL;
	if (objtype == NULL) {
		if (pval != NULL && pval != RubyNil)
			return -1;  /*  The object type must be specified unless we are removing the registration  */
	} else if (strcmp(objtype, "Molecule") == 0) {
		eventType = MyDocumentEvent;
		handler = sGetEventHandlerFromObjectAndType(obj, eventType);
		if (handler == NULL)
			return -1;  /*  obj is not an event handler  */
		if (msg == NULL)
			eventId = -1;
		else if (strcmp(msg, "documentModified") == 0)
			eventId = MyDocumentEvent_documentModified;
		else if (strcmp(msg, "documentWillClose") == 0)
			eventId = MyDocumentEvent_documentWillClose;
		else return -2; /*  this event type is not supported  */
	} else return -1;
	
	if (pval == NULL || pval == RubyNil) {
		/*  Remove the registration  */
		int ii = -3;
		for (i = 0; i < countMessageData; i++) {
			if ((messageData[i * 5] == obj || obj == NULL) && 
				(messageData[i * 5 + 1] == (void *)eventId || eventId == -1) &&
				(messageData[i * 5 + 2] == (void *)eventType || eventType == -1)) {
					if (obj == NULL)
						handler = sGetEventHandlerFromObjectAndType(messageData[i * 5], (wxEventType)messageData[i * 5 + 2]);
					handler->Disconnect(eventId, eventType, wxCommandEventHandler(RubyDialogFrame::HandleDocumentEvent), NULL, this);
					if ((wxEventType)messageData[i * 5 + 2] == MyDocumentEvent)
						MoleculeRelease((Molecule *)obj);
					messageData[i * 5] = NULL;  /*  Disabled entry  */
					ii = i;
				}
		}
		if (i == countMessageData)
			return -3;  /*  No such message  */
		return ii;
	} else {
		/*  Check the duplicate  */
		j = countMessageData;  /*  The position to store info if it is new  */
		for (i = 0; i < countMessageData; i++) {
			if (messageData[i * 5] == obj && 
				messageData[i * 5 + 1] == (void *)eventId &&
				messageData[i * 5 + 2] == (void *)eventType) {
				/*  Just replace the arguments  */
				messageData[i * 5 + 3] = (void *)oval;
				messageData[i * 5 + 4] = (void *)pval;
				break;
			}
			if (messageData[i * 5] == NULL)
				j = i;
		}
		if (i == countMessageData) {
			/*  Register the data and establish a new connection  */
			if (j == countMessageData) {
				/*  Create a new entry  */
				InsertArray(&messageData, &countMessageData, sizeof(void *) * 5, i, 1, NULL);
			}
			messageData[j * 5] = obj;
			messageData[j * 5 + 1] = (void *)eventId;
			messageData[j * 5 + 2] = (void *)eventType;
			messageData[j * 5 + 3] = (void *)oval;
			messageData[j * 5 + 4] = (void *)pval;
			handler->Connect(eventId, eventType, wxCommandEventHandler(RubyDialogFrame::HandleDocumentEvent), NULL, this);
			if (eventType == MyDocumentEvent)
				MoleculeRetain((Molecule *)obj);
			i = j;
		}
		return i;
	}
}

void
RubyDialogFrame::HandleDocumentEvent(wxCommandEvent &event)
{
	int i;
	int eventId = event.GetId();
	int eventType = event.GetEventType();
	wxObject *eventObject = event.GetEventObject();
	void *obj;

	if (eventType == MyDocumentEvent) {
		if (wxDynamicCast(eventObject, MyDocument) != NULL) {
			obj = ((MyDocument *)eventObject)->GetMolecule();
		} else return;
	} else return;
	
	/*  Look up the message table  */
	for (i = 0; i < countMessageData; i++) {
		if (messageData[i * 5] == obj && 
			messageData[i * 5 + 1] == (void *)eventId &&
			messageData[i * 5 + 2] == (void *)eventType) {
			int status;
			RubyValue oval = (RubyValue)messageData[i * 5 + 3];
			RubyValue pval = (RubyValue)messageData[i * 5 + 4];
			Ruby_funcall2_protect_extern(pval, g_RubyID_call, 1, &oval, &status);
		}
	}
	event.Skip();
}

void
RubyDialogFrame::HandlePaintEvent(wxPaintEvent &event)
{
	wxWindow *win = wxDynamicCast(event.GetEventObject(), wxWindow);
	if (win == NULL)
		return;
	wxPaintDC dc(win);
	currentContext = &dc;
	currentDrawingItem = win;
	RubyDialog_doPaintAction((RubyValue)dval, (RDItem *)win);
	currentContext = NULL;
	currentDrawingItem = NULL;
}

#pragma mark ====== MyListCtrlDataSource methods ======

int
RubyDialogFrame::GetItemCount(MyListCtrl *ctrl)
{
	return RubyDialog_GetTableItemCount((RubyValue)dval, (RDItem *)ctrl);
}

wxString 
RubyDialogFrame::GetItemText(MyListCtrl *ctrl, long row, long column) const
{
	char buf[1024];
	RubyDialog_GetTableItemText((RubyValue)dval, (RDItem *)ctrl, row, column, buf, sizeof buf);
	wxString str(buf, WX_DEFAULT_CONV);
	return str;
}

int 
RubyDialogFrame::SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value)
{
	return RubyDialog_SetTableItemText((RubyValue)dval, (RDItem *)ctrl, row, column, value.mb_str(WX_DEFAULT_CONV));
}

void 
RubyDialogFrame::DragSelectionToRow(MyListCtrl *ctrl, long row)
{
	RubyDialog_DragTableSelectionToRow((RubyValue)dval, (RDItem *)ctrl, row);
}

bool 
RubyDialogFrame::IsItemEditable(MyListCtrl *ctrl, long row, long column)
{
	return RubyDialog_IsTableItemEditable((RubyValue)dval, (RDItem *)ctrl, row, column);
}

bool 
RubyDialogFrame::IsDragAndDropEnabled(MyListCtrl *ctrl)
{
	return RubyDialog_IsTableDragAndDropEnabled((RubyValue)dval, (RDItem *)ctrl);
}

void 
RubyDialogFrame::OnSelectionChanged(MyListCtrl *ctrl)
{
	RubyDialog_OnTableSelectionChanged((RubyValue)dval, (RDItem *)ctrl);
}

int 
RubyDialogFrame::SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg)
{
	return RubyDialog_SetTableItemColor((RubyValue)dval, (RDItem *)ctrl, row, col, fg, bg);
}

int
RubyDialogFrame::HasPopUpMenu(MyListCtrl *ctrl, long row, long column, char ***menu_titles)
{
	return RubyDialog_HasPopUpMenu((RubyValue)dval, (RDItem *)ctrl, row, column, menu_titles);
}

void
RubyDialogFrame::OnPopUpMenuSelected(MyListCtrl *ctrl, long row, long column, int selected_index)
{
	RubyDialog_OnPopUpMenuSelected((RubyValue)dval, (RDItem *)ctrl, row, column, selected_index);
}

#pragma mark ====== Plain C interface ======

RubyDialog *
RubyDialogCallback_new(int style)
{
	RubyDialogFrame *dref;
	int fstyle = wxCAPTION | wxSYSTEM_MENU | wxDIALOG_NO_PARENT;
	if (style & rd_Resizable)
		fstyle |= wxMAXIMIZE_BOX | wxRESIZE_BORDER;
	if (style & rd_HasCloseBox)
		fstyle |= wxCLOSE_BOX;
	dref = new RubyDialogFrame(GetMainFrame(), -1, _T("Ruby Dialog"), wxDefaultPosition, wxDefaultSize, fstyle);
	return (RubyDialog *)dref;
}

void
RubyDialogCallback_release(RubyDialog *dref)
{
	((RubyDialogFrame *)dref)->Destroy();
}

void
RubyDialogCallback_setRubyObject(RubyDialog *dref, RubyValue val)
{
	((RubyDialogFrame *)dref)->SetRubyObject(val);
}

void
RubyDialogCallback_setWindowTitle(RubyDialog *dref, const char *title)
{
	wxString str(title, WX_DEFAULT_CONV);
	((RubyDialogFrame *)dref)->SetLabel(str);
}

int
RubyDialogCallback_runModal(RubyDialog *dref)
{
	int retval = ((RubyDialogFrame *)dref)->ShowModal();
	if (retval == wxID_OK)
		return 0;  /*  OK  */
	else return 1;  /* Cancel */
}

void
RubyDialogCallback_endModal(RubyDialog *dref, int status)
{
	((RubyDialogFrame *)dref)->StopIntervalTimer();
	if (((RubyDialogFrame *)dref)->IsModal())
		((RubyDialogFrame *)dref)->EndModal(status == 0 ? wxID_OK : wxID_CANCEL);
	else {
	/*  This function should not be used with non-modal dialogs, but just in case  */
		((RubyDialogFrame *)dref)->Close();
	}
}

void
RubyDialogCallback_destroy(RubyDialog *dref)
{
	((RubyDialogFrame *)dref)->StopIntervalTimer();  /*  May be unnecessary  */
	((RubyDialogFrame *)dref)->Destroy();
}

void
RubyDialogCallback_close(RubyDialog *dref)
{
	((RubyDialogFrame *)dref)->StopIntervalTimer();
	/*  This function should not be used with modal dialogs, but just in case  */
	if (((RubyDialogFrame *)dref)->IsModal())
		((RubyDialogFrame *)dref)->EndModal(wxID_CANCEL);
	else ((RubyDialogFrame *)dref)->Close();
}

void
RubyDialogCallback_show(RubyDialog *dref)
{
	if (((RubyDialogFrame *)dref)->myTimer != NULL)
		((RubyDialogFrame *)dref)->StartIntervalTimer(-1);
	((RubyDialogFrame *)dref)->Show(true);
	((RubyDialogFrame *)dref)->Raise();
#if defined(__WXMAC__)
	{
		extern void AddWindowsItemWithTitle(const char *title);
		wxString str = ((RubyDialogFrame *)dref)->GetLabel();
		AddWindowsItemWithTitle(str.mb_str(WX_DEFAULT_CONV));
	}
#endif
}

void
RubyDialogCallback_hide(RubyDialog *dref)
{
	((RubyDialogFrame *)dref)->StopIntervalTimer();
	((RubyDialogFrame *)dref)->Show(false);
}

int
RubyDialogCallback_startIntervalTimer(RubyDialog *dref, float interval)
{
	return ((RubyDialogFrame *)dref)->StartIntervalTimer(interval * 1000);
}

void
RubyDialogCallback_stopIntervalTimer(RubyDialog *dref)
{
	((RubyDialogFrame *)dref)->StopIntervalTimer();
}

void
RubyDialogCallback_enableOnKeyHandler(RubyDialog *dref, int flag)
{
	((RubyDialogFrame *)dref)->onKeyHandlerEnabled = (flag != 0);
}


static inline RDRect
RDRectFromwxRect(const wxRect &frame)
{
	RDRect rframe;
	rframe.origin.x = frame.x;
	rframe.origin.y = frame.y;
	rframe.size.width = frame.width;
	rframe.size.height = frame.height;
	return rframe;
}

static inline wxRect
wxRectFromRDRect(RDRect rframe)
{
	wxRect frame((int)rframe.origin.x, (int)rframe.origin.y, (int)rframe.size.width, (int)rframe.size.height);
	return frame;
}

RDSize
RubyDialogCallback_windowMinSize(RubyDialog *dref)
{
	wxSize minSize = ((RubyDialogFrame *)dref)->GetMinSize();
	RDSize rminSize;
	rminSize.width = minSize.GetWidth();
	rminSize.height = minSize.GetHeight();
	return rminSize;
}

void
RubyDialogCallback_setWindowMinSize(RubyDialog *dref, RDSize size)
{
	wxSize minSize;
	minSize.x = (int)size.width;
	minSize.y = (int)size.height;
	((RubyDialogFrame *)dref)->SetMinSize(minSize);
}

RDSize
RubyDialogCallback_windowSize(RubyDialog *dref)
{
	wxSize minSize = ((RubyDialogFrame *)dref)->GetSize();
	RDSize rminSize;
	rminSize.width = minSize.GetWidth();
	rminSize.height = minSize.GetHeight();
	return rminSize;
}
void
RubyDialogCallback_setWindowSize(RubyDialog *dref, RDSize size)
{
	wxSize wsize((int)size.width, (int)size.height);
	((RubyDialogFrame *)dref)->SetSize(wsize);
	((RubyDialogFrame *)dref)->CentreOnScreen();
}

void
RubyDialogCallback_setAutoResizeEnabled(RubyDialog *dref, int flag)
{
	((RubyDialogFrame *)dref)->SetAutoResizeEnabled(flag);
}

int
RubyDialogCallback_isAutoResizeEnabled(RubyDialog *dref)
{
	return ((RubyDialogFrame *)dref)->IsAutoResizeEnabled();
}

int
RubyDialogCallback_Listen(RubyDialog *dref, void *obj, const char *objtype, const char *msg, RubyValue oval, RubyValue pval)
{
	return ((RubyDialogFrame *)dref)->ListenToObject(obj, objtype, msg, oval, pval);
}

void
RubyDialogCallback_createStandardButtons(RubyDialog *dref, const char *oktitle, const char *canceltitle)
{
	((RubyDialogFrame *)dref)->CreateStandardButtons(oktitle, canceltitle);
}

static wxRect
OffsetForItemRect(const char *type)
{
	wxRect offset(0, 0, 0, 0);
	if (strcmp(type, "textfield") == 0) {
#if defined(__WXMAC__)
		offset.height = 4;
#endif
	}
	else if (strcmp(type, "button") == 0) {
#if defined(__WXMAC__)
		offset.width = 24;
		offset.height = 14;
#else
		offset.width = 8;
		offset.height = 0;
#endif
	} else if (strcmp(type, "checkbox") == 0) {
		offset.width = 10;
	}
	return offset;
}

RDItem *
RubyDialogCallback_createItem(RubyDialog *dref, const char *type, const char *title, RDRect frame)
{
	wxWindow *control = NULL;
	wxRect rect, offset;
	RubyDialogFrame *parent = ((RubyDialogFrame *)dref);
	wxString tstr((title ? title : ""), WX_DEFAULT_CONV);
	bool no_action = false;
	
	rect = wxRectFromRDRect(frame);
	offset = OffsetForItemRect(type);
	if (rect.width == 0 && rect.height == 0) {
		rect.SetPosition(wxDefaultPosition);
		rect.SetSize(wxDefaultSize);
	} else {	
		rect.SetX(rect.x + offset.x);
		rect.SetY(rect.y + offset.y);
		rect.SetWidth(rect.width + offset.width);
		rect.SetHeight(rect.height + offset.height);
	}

	if (strcmp(type, "text") == 0) {
		/*  Static text */		
		long style = wxST_NO_AUTORESIZE;
		if (rect.width == wxDefaultSize.x && rect.height == wxDefaultSize.y)
			style = 0;  /*  Allow autoresize  */
		wxStaticText *st = new wxStaticText(parent, -1, tstr, rect.GetPosition(), rect.GetSize(), style);
		control = st;
		no_action = true;
	} else if (strcmp(type, "textfield") == 0) {
		/*  Editable text  */
		wxTextCtrl *tc = new wxTextCtrl(parent, -1, tstr, rect.GetPosition(), rect.GetSize());
		control = tc;
		tc->Connect(-1, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "textview") == 0) {
		/*  Text view  */
		wxTextCtrl *tc = new wxTextCtrl(parent, -1, tstr, rect.GetPosition(), rect.GetSize(), wxTE_MULTILINE | wxTE_RICH);
		control = tc;
		tc->Connect(-1, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "view") == 0) {
		/*  Panel  */
		MyDrawingPanel *pn = new MyDrawingPanel(parent, -1, rect.GetPosition(), rect.GetSize());
		control = pn;
		pn->Connect(-1, wxEVT_PAINT, wxPaintEventHandler(RubyDialogFrame::HandlePaintEvent), NULL, parent);
	} else if (strcmp(type, "layout_view") == 0) {
		/*  Panel (for layout only)  */
		MyLayoutPanel *mpn = new MyLayoutPanel(parent, -1, rect.GetPosition(), rect.GetSize());
		control = mpn;
	} else if (strcmp(type, "line") == 0) {
		/*  Separator line  */
		int direction = (rect.width > rect.height ? wxLI_HORIZONTAL : wxLI_VERTICAL);
		wxStaticLine *ln = new wxStaticLine(parent, -1, rect.GetPosition(), rect.GetSize(), direction);
		control = ln;
	/*	printf("is_vertical = %d\n", (int)ln->IsVertical()); */
	} else if (strcmp(type, "button") == 0) {
		/*  Button  */
		wxButton *bn = new wxButton(parent, -1, tstr, rect.GetPosition(), rect.GetSize());
		control = bn;
		bn->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "togglebutton") == 0) {
		/*  Button  */
		wxToggleButton *bn = new wxToggleButton(parent, -1, tstr, rect.GetPosition(), rect.GetSize(), wxTOGGLEBUTTON_STYLE);
		control = bn;
		bn->Connect(-1, wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "popup") == 0) {
		/*  Popup button (wxChoice)  */
		wxString du[1] = { _T(" ") };
		wxChoice *ch = new wxChoice(parent, -1, rect.GetPosition(), rect.GetSize(), 1, du);
		control = ch;
		ch->Connect(-1, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "checkbox") == 0) {
		/*  Checkbox (wxCheckBox)  */
		wxCheckBox *cb = new wxCheckBox(parent, -1, tstr, rect.GetPosition(), rect.GetSize());
		control = cb;
		cb->Connect(-1, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "radio") == 0) {
		/*  Radio button (not grouped)  */
		wxRadioButton *rb = new wxRadioButton(parent, -1, tstr, rect.GetPosition(), rect.GetSize(), wxRB_SINGLE);
		control = rb;
		rb->Connect(-1, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "table") == 0) {
		/*  Table view = MyListCtrl  */
		MyListCtrl *tb = new MyListCtrl();
		tb->Create(parent, -1, rect.GetPosition(), rect.GetSize());
		control = tb;
		tb->SetDataSource(parent);
#if defined(__WXMSW__)
		//  On Windows, extra resizing seems necessary to properly layout the subwindows
		tb->SetSize(rect.width + 1, rect.height + 1);
		tb->SetSize(rect.width, rect.height);
#endif
	} else return NULL;
	
	if (title[0] != 0 || strcmp(type, "textfield") == 0) {
		/*  Resize the frame rect as necessary  */
		RDSize minSize = RubyDialogCallback_sizeOfString((RDItem *)control, title);
		wxSize size = control->GetSize();
		if (size.GetHeight() < minSize.height)
			size.SetHeight(minSize.height);
		if (size.GetWidth() < minSize.width)
			size.SetWidth(minSize.width);
		size.SetWidth(size.GetWidth() + offset.width);
		size.SetHeight(size.GetHeight() + offset.height);
		control->SetSize(size);
	}
	
	if (wxDynamicCast(control, wxTextCtrl) != NULL) {
		/*  Set default font  */
		wxTextAttr attr = ((wxTextCtrl *)control)->GetDefaultStyle();
		attr.SetFont(wxSystemSettings::GetFont(wxSYS_SYSTEM_FONT));
		((wxTextCtrl *)control)->SetDefaultStyle(attr);
	}

	((RubyDialogFrame *)dref)->AddDialogItem((RDItem *)control);
	
	return (RDItem *)control;
}

RDItem *
RubyDialogCallback_dialogItemAtIndex(RubyDialog *dref, int idx)
{
	if (idx == -1)
		return (RDItem *)dref;
	else return ((RubyDialogFrame *)dref)->DialogItemAtIndex(idx);
}

int
RubyDialogCallback_indexOfItem(RubyDialog *dref, RDItem *item)
{
	return ((RubyDialogFrame *)dref)->SearchDialogItem(item);
}

void
RubyDialogCallback_moveItemUnderView(RDItem *item, RDItem *superView, RDPoint origin)
{
	if (item == NULL || superView == NULL || item == superView)
		return;
	if (((wxWindow *)item)->Reparent((wxWindow *)superView)) {
		((wxWindow *)item)->Move(origin.x, origin.y);
	}
}

RDItem *
RubyDialogCallback_superview(RDItem *item)
{
	return (RDItem *)(((wxWindow *)item)->GetParent());
}

RDRect
RubyDialogCallback_frameOfItem(RDItem *item)
{
	wxRect rect = ((wxWindow *)item)->GetRect();
	if (gRubyDialogIsFlipped) {
		wxWindow *parent = ((wxWindow *)item)->GetParent();
		if (parent != NULL) {
			wxRect superRect = parent->GetRect();
			rect.SetY(superRect.GetHeight() - rect.GetHeight() - rect.GetY());
		}
	}
	return RDRectFromwxRect(rect);
}

void
RubyDialogCallback_setFrameOfItem(RDItem *item, RDRect rect)
{
	wxRect wrect = wxRectFromRDRect(rect);
	if (gRubyDialogIsFlipped) {
		wxWindow *parent = ((wxWindow *)item)->GetParent();
		if (parent != NULL) {
			wxRect srect = parent->GetRect();
			wrect.SetY(srect.GetHeight() - wrect.GetHeight() - wrect.GetY());
		}
	}
	((wxWindow *)item)->SetSize(wrect);
}

void
RubyDialogCallback_setStringToItem(RDItem *item, const char *s)
{
	wxString str(s, WX_DEFAULT_CONV);
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		((wxTextCtrl *)item)->ChangeValue(str);
	//	((wxTextCtrl *)item)->Clear();
	//	((wxTextCtrl *)item)->AppendText(str);
	} else if (wxDynamicCast((wxWindow *)item, wxStaticText) != NULL) {
		((wxStaticText *)item)->SetLabel(str);
	}
}

void
RubyDialogCallback_getStringFromItem(RDItem *item, char *buf, int bufsize)
{
	wxString str;
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		str = ((wxTextCtrl *)item)->GetValue();
	} else if (wxDynamicCast((wxWindow *)item, wxStaticText) != NULL) {
			str = ((wxStaticText *)item)->GetLabel();
	} else {
		buf[0] = 0;
		return;
	}
	strncpy(buf, str.mb_str(WX_DEFAULT_CONV), bufsize - 1);
	buf[bufsize - 1] = 0;
}

char *
RubyDialogCallback_getStringPtrFromItem(RDItem *item)
{
	wxString str;
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		str = ((wxTextCtrl *)item)->GetValue();
	} else {
		return NULL;
	}
	return strdup(str.mb_str(WX_DEFAULT_CONV));
}

char *
RubyDialogCallback_titleOfItem(RDItem *item)
{
	wxString str;
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		str = ((wxTextCtrl *)item)->GetValue();
	} else {
		str = ((wxWindow *)item)->GetLabel();
	}
	return strdup(str.mb_str(WX_DEFAULT_CONV));
}

void
RubyDialogCallback_setTitleToItem(RDItem *item, const char *s)
{
	wxString str(s, WX_DEFAULT_CONV);
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		((wxTextCtrl *)item)->SetValue(str);
	} else {
		((wxWindow *)item)->SetLabel(str);
	}
}

void
RubyDialogCallback_setEnabledForItem(RDItem *item, int flag)
{
	((wxWindow *)item)->Enable(flag);
}

int
RubyDialogCallback_isItemEnabled(RDItem *item)
{
/*	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL)
		return ((wxTextCtrl *)item)->IsEditable();
	else */
	return ((wxWindow *)item)->IsEnabled();
}

void
RubyDialogCallback_setEditableForItem(RDItem *item, int flag)
{
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL)
		return ((wxTextCtrl *)item)->SetEditable(flag);
}

int
RubyDialogCallback_isItemEditable(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL)
		return ((wxTextCtrl *)item)->IsEditable();
	else return 0;
}

void
RubyDialogCallback_setStateForItem(RDItem *item, int state)
{
	if (wxDynamicCast((wxWindow *)item, wxRadioButton) != NULL) {
		((wxRadioButton *)item)->SetValue(state);
	} else if (wxDynamicCast((wxWindow *)item, wxCheckBox) != NULL) {
		((wxCheckBox *)item)->SetValue(state);
	} else if (wxDynamicCast((wxWindow *)item, wxToggleButton) != NULL) {
		((wxToggleButton *)item)->SetValue(state);
	}
}

int
RubyDialogCallback_getStateForItem(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, wxRadioButton) != NULL) {
		return ((wxRadioButton *)item)->GetValue();
	} else if (wxDynamicCast((wxWindow *)item, wxCheckBox) != NULL) {
		return ((wxCheckBox *)item)->GetValue();
	} else if (wxDynamicCast((wxWindow *)item, wxToggleButton) != NULL) {
		return ((wxToggleButton *)item)->GetValue();
	} else return -1;
}

void
RubyDialogCallback_setHiddenForItem(RDItem *item, int flag)
{
	((wxWindow *)item)->Show(flag == 0);
}

int
RubyDialogCallback_isItemHidden(RDItem *item)
{
	return !(((wxWindow *)item)->IsShown());
}

void
RubyDialogCallback_setFontForItem(RDItem *item, int size, int family, int style, int weight)
{
	wxTextCtrl *textctrl;
	wxControl *ctrl;
	wxFont font;
	if ((textctrl = wxDynamicCast((wxWindow *)item, wxTextCtrl)) != NULL) {
		wxTextAttr attr = textctrl->GetDefaultStyle();
		font = attr.GetFont();
	} else if ((ctrl = wxDynamicCast((wxWindow *)item, wxControl)) != NULL) {
		font = ctrl->GetFont();
	} else return;
	if (size == 0)
		size = font.GetPointSize();
	if (family == 0)
		family = font.GetFamily();
	else {
		family = (family == 2 ? wxFONTFAMILY_ROMAN :
				  (family == 3 ? wxFONTFAMILY_SWISS :
				   (family == 4 ? wxFONTFAMILY_MODERN :
					wxFONTFAMILY_DEFAULT)));
	}
	if (style == 0)
		style = font.GetStyle();
	else {
		style = (style == 2 ? wxFONTSTYLE_SLANT :
				 (style == 3 ? wxFONTSTYLE_ITALIC :
				  wxFONTSTYLE_NORMAL));
	}
	if (weight == 0)
		weight = font.GetWeight();
	else {
		weight = (weight == 2 ? wxFONTWEIGHT_BOLD :
				  (weight == 3 ? wxFONTWEIGHT_LIGHT :
				   wxFONTWEIGHT_NORMAL));
	}
	if (textctrl != NULL) {
		wxTextAttr newAttr;
		newAttr.SetFont(wxFont(size, family, style, weight));
		textctrl->SetDefaultStyle(newAttr);
	} else {
		ctrl->SetFont(wxFont(size, family, style, weight));
		wxString label = ctrl->GetLabel();
		ctrl->SetLabel(_(""));
		ctrl->SetLabel(label);  /*  Update the control size  */
	}
}

int
RubyDialogCallback_getFontForItem(RDItem *item, int *size, int *family, int *style, int *weight)
{
	int n;
	wxTextCtrl *ctrl;
	if ((ctrl = wxDynamicCast((wxWindow *)item, wxTextCtrl)) != NULL) {
		wxTextAttr attr = ctrl->GetDefaultStyle();
		wxFont font = attr.GetFont();
		if (size != NULL)
			*size = font.GetPointSize();
		if (family != NULL) {
			n = font.GetFamily();
			*family = (n == wxFONTFAMILY_DEFAULT ? 1 :
					   (n == wxFONTFAMILY_ROMAN ? 2 :
						(n == wxFONTFAMILY_SWISS ? 3 :
						 (n == wxFONTFAMILY_MODERN ? 4 :
						  0))));
		}
		if (style != NULL) {
			n = font.GetStyle();
			*style = (n == wxFONTSTYLE_NORMAL ? 1 :
					  (n == wxFONTSTYLE_SLANT ? 2 :
					   (n == wxFONTSTYLE_ITALIC ? 3 :
						0)));
		}
		if (weight != NULL) {
			n = font.GetWeight();
			*weight = (n == wxFONTWEIGHT_NORMAL ? 1 :
					   (n == wxFONTWEIGHT_BOLD ? 2 :
						(n == wxFONTWEIGHT_LIGHT ? 3 :
						 0)));
		}
		return 1;
	} else return 0;
}

void
RubyDialogCallback_setForegroundColorForItem(RDItem *item, const double *col)
{
	wxColour wcol((int)(col[0] * 255), (int)(col[1] * 255), (int)(col[2] * 255), (int)(col[3] * 255));
	((wxWindow *)item)->SetForegroundColour(wcol);
}

void
RubyDialogCallback_setBackgroundColorForItem(RDItem *item, const double *col)
{
	wxColour wcol((int)(col[0] * 255), (int)(col[1] * 255), (int)(col[2] * 255), (int)(col[3] * 255));
	((wxWindow *)item)->SetBackgroundColour(wcol);
}

void
RubyDialogCallback_getForegroundColorForItem(RDItem *item, double *col)
{
	wxColour wcol = ((wxWindow *)item)->GetForegroundColour();
	col[0] = wcol.Red() / 255.0;
	col[1] = wcol.Green() / 255.0;
	col[2] = wcol.Blue() / 255.0;
	col[3] = wcol.Alpha() / 255.0;
}

void
RubyDialogCallback_getBackgroundColorForItem(RDItem *item, double *col)
{
	wxColour wcol = ((wxWindow *)item)->GetBackgroundColour();
	col[0] = wcol.Red() / 255.0;
	col[1] = wcol.Green() / 255.0;
	col[2] = wcol.Blue() / 255.0;
	col[3] = wcol.Alpha() / 255.0;
}

int
RubyDialogCallback_appendString(RDItem *item, const char *str)
{
	wxTextCtrl *ctrl;
	if ((ctrl = wxDynamicCast((wxWindow *)item, wxTextCtrl)) != NULL) {
		ctrl->AppendText(wxString(str, WX_DEFAULT_CONV));
		return 1;
	} else return 0;
}

void
RubyDialogCallback_setNeedsDisplay(RDItem *item, int flag)
{
	if (flag)
		((wxWindow *)item)->Refresh();
}

void
RubyDialogCallback_setNeedsDisplayInRect(RDItem *item, RDRect rect, int eraseBackground)
{
	wxRect wrect = wxRectFromRDRect(rect);
	((wxWindow *)item)->RefreshRect(wrect, eraseBackground);
}

int
RubyDialogCallback_countSubItems(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL) {
		return ((wxChoice *)item)->GetCount();
	} else return 0;
}

int
RubyDialogCallback_appendSubItem(RDItem *item, const char *s)
{
	wxString str(s, WX_DEFAULT_CONV);
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL) {
		return ((wxChoice *)item)->Append(str);
	} else return -1;
}

int
RubyDialogCallback_insertSubItem(RDItem *item, const char *s, int pos)
{
	wxString str(s, WX_DEFAULT_CONV);
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL && pos >= 0 && pos < ((wxChoice *)item)->GetCount()) {
		return ((wxChoice *)item)->Insert(str, pos);
	} else return -1;
}

int
RubyDialogCallback_deleteSubItem(RDItem *item, int pos)
{	
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL && pos >= 0 && pos < ((wxChoice *)item)->GetCount()) {
		((wxChoice *)item)->Delete(pos);
		return pos;
	} else return -1;
}

char *
RubyDialogCallback_titleOfSubItem(RDItem *item, int pos)
{
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL && pos >= 0 && pos < ((wxChoice *)item)->GetCount()) {
		wxString str = ((wxChoice *)item)->GetString(pos);
		return strdup(str.mb_str(WX_DEFAULT_CONV));
	} else return NULL;
}

void
RubyDialogCallback_setSelectedSubItem(RDItem *item, int pos)
{
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL && pos >= -1 && pos < ((wxChoice *)item)->GetCount()) {
		if (pos == -1)
			pos = wxNOT_FOUND;
		((wxChoice *)item)->SetSelection(pos);
	}
}

int
RubyDialogCallback_selectedSubItem(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL) {
		return ((wxChoice *)item)->GetSelection();
	} else return -1;
}

RDSize
RubyDialogCallback_sizeOfString(RDItem *item, const char *s)
{
	RDSize size;
	wxCoord w, h, descent, leading;
	wxClientDC dc((wxWindow *)item);
	int len;
	const char *s1, *s2, *sfin;
	size.width = size.height = 0;
	s1 = (s == NULL || s[0] == 0 ? " " : s);
	len = strlen(s1);
	sfin = s1 + len;
	while (1) {
		s2 = strchr(s1, '\n');
		if (s2 == NULL)
			s2 = sfin;
		wxString str(s1, WX_DEFAULT_CONV, s2 - s1);
		dc.GetTextExtent(str, &w, &h, &descent, &leading);
		if (size.width < w)
			size.width = w;
		size.height += h + descent + leading;
		if (s2 >= sfin)
			break;
		s1 = s2 + 1;
	}
	return size;
}

RDSize
RubyDialogCallback_resizeToBest(RDItem *item)
{
	wxSize size;
	RDSize rsize;
	size = ((wxWindow *)item)->GetBestSize();
	((wxWindow *)item)->SetSize(size);
	rsize.width = size.GetWidth();
	rsize.height = size.GetHeight();
	return rsize;
}

char
RubyDialogCallback_deleteTableColumn(RDItem *item, int col)
{
	if (wxDynamicCast((wxWindow *)item, MyListCtrl) != NULL) {
		return ((MyListCtrl *)item)->DeleteColumn(col);
	} else return false;
}

char
RubyDialogCallback_insertTableColumn(RDItem *item, int col, const char *heading, int format, int width)
{
	if (wxDynamicCast((wxWindow *)item, MyListCtrl) != NULL) {
		wxString hstr((heading ? heading : ""), WX_DEFAULT_CONV);
		return ((MyListCtrl *)item)->InsertColumn(col, hstr, format, width);
	} else return false;
}

int
RubyDialogCallback_countTableColumn(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, MyListCtrl) != NULL) {
		return ((MyListCtrl *)item)->GetColumnCount();
	} else return -1;
}

char
RubyDialogCallback_isTableRowSelected(RDItem *item, int row)
{
	if (wxDynamicCast((wxWindow *)item, MyListCtrl) != NULL) {
		return ((MyListCtrl *)item)->GetItemState(row, wxLIST_STATE_SELECTED) != 0;
	} else return false;
}

char
RubyDialogCallback_setTableRowSelected(RDItem *item, int row, int flag)
{
	if (wxDynamicCast((wxWindow *)item, MyListCtrl) != NULL) {
		long state = (flag ? wxLIST_STATE_SELECTED : 0);
		return ((MyListCtrl *)item)->SetItemState(row, state, wxLIST_STATE_SELECTED);
	} else return false;
}

void
RubyDialogCallback_refreshTable(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, MyListCtrl) != NULL) {
		((MyListCtrl *)item)->RefreshTable();
	}
}

int
RubyDialogCallback_savePanel(const char *title, const char *dirname, const char *wildcard, char *buf, int bufsize)
{
	int result;
	wxString pstr((dirname ? dirname : ""), WX_DEFAULT_CONV);
	wxString tstr((title ? title : "Choose a file"), WX_DEFAULT_CONV);
	wxString fstr(buf, WX_DEFAULT_CONV);
	wxString wstr((wildcard ? wildcard : "All files (*.*)|*.*"), WX_DEFAULT_CONV);
	wxWindow *old_focus = wxWindow::FindFocus();
	wxFileDialog *dialog = new wxFileDialog(NULL, tstr, pstr, fstr, wstr, wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (dialog->ShowModal() == wxID_OK) {
		strncpy(buf, dialog->GetPath().mb_str(wxConvFile), bufsize - 1);
		buf[bufsize - 1] = 0;
		result = 1;
	} else {
		buf[0] = 0;
		result = 0;
	}
	dialog->Destroy();
	if (old_focus != NULL)
		old_focus->SetFocus();
	return result;
}

int
RubyDialogCallback_openPanel(const char *title, const char *dirname, const char *wildcard, char ***array, int for_directories, int multiple_selection)
{
	int result = 0;
	wxString pstr((dirname ? dirname : ""), WX_DEFAULT_CONV);
	wxString wstr((wildcard ? wildcard : "All files (*.*)|*.*"), WX_DEFAULT_CONV);
	int style = wxFD_OPEN | (multiple_selection ? wxFD_MULTIPLE : 0);
	wxWindow *old_focus = wxWindow::FindFocus();
	if (for_directories) {
		wxString tstr((title ? title : "Choose a directory"), WX_DEFAULT_CONV);
		wxDirDialog *dialog = new wxDirDialog(NULL, tstr, pstr);
		if (dialog->ShowModal() == wxID_OK) {
			*array = (char **)malloc(sizeof(char *));
			(*array)[0] = strdup(dialog->GetPath().mb_str(wxConvFile));
			result = 1;
		}
		dialog->Destroy();
	} else {
		wxString tstr((title ? title : "Choose a file"), WX_DEFAULT_CONV);
		wxFileDialog *dialog = new wxFileDialog(NULL, tstr, pstr, _T(""), wstr, style);
		if (dialog->ShowModal() == wxID_OK) {
			if (multiple_selection) {
				int i, n;
				wxArrayString paths;
				dialog->GetPaths(paths);
				n = paths.GetCount();
				*array = (char **)malloc(sizeof(char *) * n);
				for (i = 0; i < n; i++) {
					(*array)[i] = strdup(paths[i].mb_str(wxConvFile));
				}
				result = n;
			} else {
				*array = (char **)malloc(sizeof(char *));
				(*array)[0] = strdup(dialog->GetPath().mb_str(wxConvFile));
				result = 1;
			}
		}
		dialog->Destroy();
	}
	if (old_focus != NULL)
		old_focus->SetFocus();
	return result;
}

#pragma mark ====== Plain C Interface (Device Context) ======

RDDeviceContext *
RubyDialogCallback_getDeviceContextForRubyDialog(RubyDialog *dref)
{
	return (RDDeviceContext *)(((RubyDialogFrame *)dref)->currentContext);
}

void
RubyDialogCallback_clear(RDDeviceContext *dc)
{
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	dcp->Clear();
}

void
RubyDialogCallback_drawEllipse(RDDeviceContext *dc, float x, float y, float r1, float r2)
{
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	dcp->DrawEllipse(x - r1, y - r2, r1 * 2, r2 * 2);
}

void
RubyDialogCallback_drawLine(RDDeviceContext *dc, int ncoords, float *coords)
{
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	wxPoint *pts = new wxPoint[ncoords];
	int i;
	for (i = 0; i < ncoords; i++) {
		pts[i].x = (int)coords[i * 2];
		pts[i].y = (int)coords[i * 2 + 1];
	}
	dcp->DrawLines(ncoords, pts);
	delete [] pts;
}

void
RubyDialogCallback_drawRectangle(RDDeviceContext *dc, float x, float y, float width, float height, float round)
{
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	if (round > 0.0)
		dcp->DrawRoundedRectangle((int)x, (int)y, (int)width, (int)height, (int)round);
	else
		dcp->DrawRectangle((int)x, (int)y, (int)width, (int)height);
}

void
RubyDialogCallback_drawText(RDDeviceContext *dc, const char *s, float x, float y)
{
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	wxString str(s, WX_DEFAULT_CONV);
	dcp->DrawText(str, (int)x, (int)y);
}

void
RubyDialogCallback_setFont(RDDeviceContext *dc, void **args)
{
	int i, j;
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	wxFont font = dcp->GetFont();
	for (i = 0; args[i] != NULL; i += 2) {
		if (strcmp((const char *)args[i], "size") == 0) {
			float size = *((float *)(args[i + 1]));
			font.SetPointSize((int)size);
		} else if (strcmp((const char *)args[i], "style") == 0) {
			int style = (int)(args[i + 1]);
			switch (style) {
				case 0: style = wxFONTSTYLE_NORMAL; break;
				case 1: style = wxFONTSTYLE_ITALIC; break;
				case 2: style = wxFONTSTYLE_SLANT; break;
				default: style = wxFONTSTYLE_NORMAL; break;
			}
			font.SetStyle(style);
		} else if (strcmp((const char *)args[i], "family") == 0) {
			wxFontFamily family;
			j = (int)(args[i + 1]);
			switch (j) {
				case 0: family = wxFONTFAMILY_DEFAULT; break;
				case 1: family = wxFONTFAMILY_ROMAN; break;
				case 2: family = wxFONTFAMILY_SWISS; break;
				case 3: family = wxFONTFAMILY_MODERN; break;
				default: family = wxFONTFAMILY_DEFAULT; break;
			}
			font.SetFamily(family);
		} else if (strcmp((const char *)args[i], "weight") == 0) {
			wxFontWeight weight;
			j = (int)(args[i + 1]);
			switch (j) {
				case 0: weight = wxFONTWEIGHT_NORMAL; break;
				case 1: weight = wxFONTWEIGHT_LIGHT; break;
				case 2: weight = wxFONTWEIGHT_BOLD; break;
				default: weight = wxFONTWEIGHT_NORMAL; break;
			}
			font.SetWeight(weight);
		}
	}
	dcp->SetFont(font);
}

void
RubyDialogCallback_setPen(RDDeviceContext *dc, void **args)
{
	int i;
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	wxPen pen = wxNullPen;
	if (args != NULL) {
		pen = dcp->GetPen();
		for (i = 0; args[i] != NULL; i += 2) {
			if (strcmp((const char *)args[i], "color") == 0) {
				float *fp = (float *)args[i + 1];
				wxColour col((int)(fp[0] * 255.0), (int)(fp[1] * 255.0), (int)(fp[2] * 255.0), (int)(fp[3] * 255.0));
				pen.SetColour(col);
			} else if (strcmp((const char *)args[i], "width") == 0) {
				float width = *((float *)(args[i + 1]));
				pen.SetWidth((int)width);
			} else if (strcmp((const char *)args[i], "style") == 0) {
				int style = (int)(args[i + 1]);
				switch (style) {
					case 0: style = wxSOLID; break;
					case 1: style = wxTRANSPARENT; break;
					case 2: style = wxDOT; break;
					case 3: style = wxLONG_DASH; break; 
					case 4: style = wxSHORT_DASH; break;
					case 5: style = wxDOT_DASH; break;
					default: style = wxSOLID; break;
				}
				pen.SetStyle(style);
			}
		}
	}
	dcp->SetPen(pen);
}

void
RubyDialogCallback_setBrush(RDDeviceContext *dc, void **args)
{
	int i;
	wxDC *dcp = (wxDC *)dc;
	if (dcp == NULL)
		return;
	wxBrush brush = wxNullBrush;
	if (args != NULL) {
		brush = dcp->GetBrush();
		for (i = 0; args[i] != NULL; i += 2) {
			if (strcmp((const char *)args[i], "color") == 0) {
				float *fp = (float *)args[i + 1];
				wxColour col((int)(fp[0] * 255.0), (int)(fp[1] * 255.0), (int)(fp[2] * 255.0), (int)(fp[3] * 255.0));
				brush.SetColour(col);
			}
		}
	}
	dcp->SetBrush(brush);
}

#pragma mark ====== Bitmap ======

RDBitmap *
RubyDialogCallback_createBitmap(int width, int height, int depth)
{
	wxBitmap *bitmap = new wxBitmap(width, height, depth);
	return (RDBitmap *)bitmap;
}

void
RubyDialogCallback_releaseBitmap(RDBitmap *bitmap)
{
	if (bitmap != NULL)
		delete (wxBitmap *)bitmap;
}


/*  Set focus on a bitmap and execute the given function  */
static RDBitmap *s_temp_dc_pointer = NULL;
int
RubyDialogCallback_executeWithFocusOnBitmap(RDBitmap *bitmap, void (*callback)(void *), void *ptr)
{
	wxMemoryDC temp_dc;
	if (s_temp_dc_pointer != NULL)
		return -1;  /*  Recursive call is not allowed  */
	s_temp_dc_pointer = (RDBitmap *)(&temp_dc);
	temp_dc.SelectObject(*((wxBitmap *)bitmap));
	(*callback)(ptr);
	temp_dc.SelectObject(wxNullBitmap);
	s_temp_dc_pointer = NULL;
	return 0;
}

RDDeviceContext *
RubyDialogCallback_getDeviceContextForBitmap(RDBitmap *bitmap)
{
	return (RDDeviceContext *)s_temp_dc_pointer;
}

int
RubyDialogCallback_saveBitmapToFile(RDBitmap *bitmap, const char *fname)
{
	int len = strlen(fname);
	wxBitmapType type = wxBITMAP_TYPE_PNG;
	MyAppCallback_initImageHandlers();
	if (len >= 4) {
		if (strcasecmp(fname + len - 4, ".png") == 0)
			type = wxBITMAP_TYPE_PNG;
		else if (strcasecmp(fname + len - 4, ".tif") == 0 || (len >= 5 && strcasecmp(fname + len - 5, ".tiff") == 0))
			type = wxBITMAP_TYPE_TIF;
	}
	wxString fn(fname, wxConvFile);
	return ((wxBitmap *)bitmap)->SaveFile(fname, type);
}
