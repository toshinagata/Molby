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
#include "wx/filedlg.h"
#include "wx/dirdlg.h"
#include "wx/dcclient.h"
#include "wx/choice.h"
#include "wx/checkbox.h"
#include "wx/radiobut.h"
#include "wx/statline.h"

#include "RubyDialogFrame.h"
#include "MyApp.h"

BEGIN_EVENT_TABLE(RubyDialogFrame, wxDialog)
//    EVT_TEXT_ENTER(-1, ConsoleFrame::OnEnterPressed)
//	EVT_CHAR(ConsoleFrame::OnChar)
//	EVT_RICHTEXT_RETURN(-1, ConsoleFrame::OnEnterPressed)
END_EVENT_TABLE()

RubyDialogFrame::RubyDialogFrame(wxWindow* parent, wxWindowID wid, const wxString& title, const wxPoint& pos, const wxSize& size, long style):
	wxDialog(parent, wid, title, pos, size, style)
{
	ditems = NULL;
	nditems = 0;
	dval = NULL;
	
	//  Create a vertical box sizer that contains a panel containing all controls and a sizer containing
	//  OK/Cancel buttons
	contentSizer = new wxBoxSizer(wxVERTICAL);
	contentPanel = NULL;
	buttonSizer = NULL;  //  Will be created later
	boxSizer = new wxBoxSizer(wxVERTICAL);
	boxSizer->Add(contentSizer, 1, wxALL | wxEXPAND, 14);
	this->SetSizer(boxSizer);
	boxSizer->Layout();
	this->CentreOnScreen();
}

RubyDialogFrame::~RubyDialogFrame()
{
	if (ditems != NULL)
		free(ditems);
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
	if (item != NULL && ((wxWindow *)item)->IsKindOf(CLASSINFO(wxPanel))) {
		wxSize size = ((wxPanel *)item)->GetSize();
		if (contentPanel == NULL)
			contentSizer->Add((wxPanel *)item, 1, wxEXPAND);
		else
			contentSizer->Replace(contentPanel, (wxPanel *)item);
		contentSizer->SetItemMinSize((wxPanel *)item, size.GetWidth(), size.GetHeight());
		contentPanel = (wxPanel *)item;
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
	} else {
		if (oktitle[0] != 0) {
			wxString label1(oktitle, wxConvUTF8);
			((wxButton *)ditems[0])->SetLabel(label1);
		}
		((wxWindow *)ditems[0])->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, this);
	}
	if (canceltitle == NULL) {
		((wxWindow *)ditems[1])->Show(false);
	} else {
		if (canceltitle[0] != 0) {
			wxString label2(canceltitle, wxConvUTF8);
			((wxButton *)ditems[1])->SetLabel(label2);
		}
		((wxWindow *)ditems[1])->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, this);
	}
}

void
RubyDialogFrame::OnDialogItemAction(wxCommandEvent &event)
{
	RubyDialog_doItemAction((RubyValue)dval, (RDItem *)(event.GetEventObject()));
}

#pragma mark ====== Plain C interface ======

RubyDialog *
RubyDialogCallback_new(void)
{
	/*  RubyDialogFrame should not have a close box  */
	RubyDialogFrame *dref = new RubyDialogFrame(GetMainFrame(), -1, _T("Ruby Dialog"), wxDefaultPosition, wxDefaultSize, wxCAPTION | wxSYSTEM_MENU);
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
	wxString str(title, wxConvUTF8);
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
	((RubyDialogFrame *)dref)->EndModal(status == 0 ? wxID_OK : wxID_CANCEL);
}

void
RubyDialogCallback_close(RubyDialog *dref)
{
	((RubyDialogFrame *)dref)->Close();
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
	wxRect frame(rframe.origin.x, rframe.origin.y, rframe.size.width, rframe.size.height);
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
RubyDialogCallback_setWindowSize(RubyDialog *dref, RDSize size)
{
	wxSize wsize(size.width, size.height);
	((RubyDialogFrame *)dref)->SetSize(wsize);
	((RubyDialogFrame *)dref)->CentreOnScreen();
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
	if (strcmp(type, "textfield") == 0)
		offset.height = 5;
	else if (strcmp(type, "button") == 0) {
		offset.width = 24;
		offset.height = 14;
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
	wxString tstr((title ? title : ""), wxConvUTF8);
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
		wxStaticText *st = new wxStaticText(parent, -1, tstr, rect.GetPosition(), rect.GetSize(), wxST_NO_AUTORESIZE);
		control = st;
		no_action = true;
	} else if (strcmp(type, "textfield") == 0) {
		/*  Editable text  */
		wxTextCtrl *tc = new wxTextCtrl(parent, -1, tstr, rect.GetPosition(), rect.GetSize());
		control = tc;
		tc->Connect(-1, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "textview") == 0) {
		/*  Text view  */
		wxTextCtrl *tc = new wxTextCtrl(parent, -1, tstr, rect.GetPosition(), rect.GetSize(), wxTE_MULTILINE);
		control = tc;
		tc->Connect(-1, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
	} else if (strcmp(type, "view") == 0) {
		/*  Panel  */
		wxPanel *pn = new wxPanel(parent, -1, rect.GetPosition(), rect.GetSize());
		control = pn;
	} else if (strcmp(type, "line") == 0) {
		/*  Separator line  */
		int direction = (frame.size.width > frame.size.height ? wxLI_HORIZONTAL : wxLI_VERTICAL);
		wxStaticLine *ln = new wxStaticLine(parent, -1, rect.GetPosition(), rect.GetSize(), direction);
		control = ln;
	/*	printf("is_vertical = %d\n", (int)ln->IsVertical()); */
	} else if (strcmp(type, "button") == 0) {
		/*  Button  */
		wxButton *bn = new wxButton(parent, -1, tstr, rect.GetPosition(), rect.GetSize());
		control = bn;
		bn->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
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
		wxRadioButton *rb = new wxRadioButton(parent, -2, tstr, rect.GetPosition(), rect.GetSize(), wxRB_SINGLE);
		control = rb;
		rb->Connect(-1, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(RubyDialogFrame::OnDialogItemAction), NULL, parent);
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
	wxString str(s, wxConvUTF8);
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		((wxTextCtrl *)item)->SetValue(str);
	}
}

void
RubyDialogCallback_getStringFromItem(RDItem *item, char *buf, int bufsize)
{
	wxString str;
	if (wxDynamicCast((wxWindow *)item, wxTextCtrl) != NULL) {
		str = ((wxTextCtrl *)item)->GetValue();
	} else {
		buf[0] = 0;
		return;
	}
	strncpy(buf, str.mb_str(wxConvUTF8), bufsize - 1);
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
	return strdup(str.mb_str(wxConvUTF8));
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
	return strdup(str.mb_str(wxConvUTF8));
}

void
RubyDialogCallback_setTitleToItem(RDItem *item, const char *s)
{
	wxString str(s, wxConvUTF8);
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
	}
}

int
RubyDialogCallback_getStateForItem(RDItem *item)
{
	if (wxDynamicCast((wxWindow *)item, wxRadioButton) != NULL) {
		return ((wxRadioButton *)item)->GetValue();
	} else if (wxDynamicCast((wxWindow *)item, wxCheckBox) != NULL) {
		return ((wxCheckBox *)item)->GetValue();
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
RubyDialogCallback_setNeedsDisplay(RDItem *item, int flag)
{
	if (flag)
		((wxWindow *)item)->Refresh();
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
	wxString str(s, wxConvUTF8);
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL) {
		return ((wxChoice *)item)->Append(str);
	} else return -1;
}

int
RubyDialogCallback_insertSubItem(RDItem *item, const char *s, int pos)
{
	wxString str(s, wxConvUTF8);
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
		return strdup(str.mb_str(wxConvUTF8));
	} else return NULL;
}

void
RubyDialogCallback_setSelectedSubItem(RDItem *item, int pos)
{
	if (wxDynamicCast((wxWindow *)item, wxChoice) != NULL && pos >= 0 && pos < ((wxChoice *)item)->GetCount()) {
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
	wxPaintDC dc((wxWindow *)item);
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
		wxString str(s1, wxConvUTF8, s2 - s1);
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

int
RubyDialogCallback_savePanel(const char *title, const char *dirname, const char *wildcard, char *buf, int bufsize)
{
	int result;
	wxString pstr((dirname ? dirname : ""), wxConvUTF8);
	wxString tstr((title ? title : "Choose a file"), wxConvUTF8);
	wxString fstr(buf, wxConvUTF8);
	wxString wstr((wildcard ? wildcard : "All files (*.*)|*.*"), wxConvUTF8);
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
	return result;
}

int
RubyDialogCallback_openPanel(const char *title, const char *dirname, const char *wildcard, char ***array, int for_directories, int multiple_selection)
{
	int result = 0;
	wxString pstr((dirname ? dirname : ""), wxConvUTF8);
	wxString wstr((wildcard ? wildcard : "All files (*.*)|*.*"), wxConvUTF8);
	int style = wxFD_OPEN | (multiple_selection ? wxFD_MULTIPLE : 0);
	if (for_directories) {
		wxString tstr((title ? title : "Choose a directory"), wxConvUTF8);
		wxDirDialog *dialog = new wxDirDialog(NULL, tstr, pstr);
		if (dialog->ShowModal() == wxID_OK) {
			*array = (char **)malloc(sizeof(char *));
			(*array)[0] = strdup(dialog->GetPath().mb_str(wxConvFile));
			result = 1;
		}
		dialog->Destroy();
	} else {
		wxString tstr((title ? title : "Choose a file"), wxConvUTF8);
		wxFileDialog *dialog = new wxFileDialog(NULL, tstr, pstr, _T(""), wstr, style);
		if (dialog->ShowModal() == wxID_OK) {
			*array = (char **)malloc(sizeof(char *));
			(*array)[0] = strdup(dialog->GetPath().mb_str(wxConvFile));
			result = 1;			
		}
		dialog->Destroy();
	}
	return result;
}
