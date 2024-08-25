/*
 *  ProgressFrame.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 09/07/15.
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

#include "ProgressFrame.h"
#include "MyApp.h"

#include "wx/stattext.h"
#include "wx/gauge.h"
#include "wx/sizer.h"
#include "wx/button.h"
#include "wx/evtloop.h"

BEGIN_EVENT_TABLE(ProgressFrame, wxFrame)
EVT_BUTTON(-1, ProgressFrame::OnButtonPressed)
END_EVENT_TABLE()

int ProgressFrame::c_uniqueID = 1000;
wxWindowList ProgressFrame::c_progressFrames;

ProgressFrame::ProgressFrame(const wxString &title, const wxString &mes):
	wxFrame(NULL, c_uniqueID, title, wxDefaultPosition, wxDefaultSize, wxCAPTION)
{
	//  Vertical sizer containing (1) message text, (2) progress gauge, (3) note text
	wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
	
	m_messageText = new wxStaticText(this, -1, wxT("Message"), wxDefaultPosition, wxSize(240, 40), wxST_NO_AUTORESIZE);
	sizer->Add(m_messageText, 0, wxALL | wxEXPAND, 10);   // Can expand horizontally
	
	m_progressGauge = new wxGauge(this, -1, 10000, wxDefaultPosition, wxSize(240, 24), wxGA_HORIZONTAL);
	sizer->Add(m_progressGauge, 0, wxALL | wxEXPAND, 10);
	
  m_cancelButton = new wxButton(this, -1, _T("Cancel"), wxDefaultPosition, wxDefaultSize);
  sizer->Add(m_cancelButton, 0, wxALL | wxEXPAND, 10);
  
	m_value = -1.0;
	m_progressGauge->Pulse();
	m_messageText->SetLabel(mes);
	
	m_interruptValue = 0;

  c_progressFrames.push_back(this);
  c_uniqueID++;
  
	sizer->Layout();
	this->SetSizerAndFit(sizer);
  
  m_disabler = new wxWindowDisabler(this);  //  Make this dialog application modal

  this->Centre();
	this->Show();
  this->Enable();

}

ProgressFrame::~ProgressFrame()
{
  /*  Remove this from the progressFrame list  */
  wxWindowList::iterator itr;
  for (itr = c_progressFrames.begin(); itr != c_progressFrames.end(); itr++) {
    if (*itr == this) {
      c_progressFrames.erase(itr);
      break;
    }
  }
  delete m_disabler;
}

ProgressFrame *
ProgressFrame::FindProgressFrameWithID(int id)
{
  if (id < 0) {
    /*  Returns the last one  */
    if (c_progressFrames.size() > 0) {
      wxWindow *w = c_progressFrames.back();
      return wxDynamicCast(w, ProgressFrame);
    } else return NULL;
  } else {
    wxWindowList::iterator itr;
    for (itr = c_progressFrames.begin(); itr != c_progressFrames.end(); itr++) {
      if ((*itr)->GetId() == id)
        return wxDynamicCast(*itr, ProgressFrame);
    }
    return NULL;  /*  Not found  */
  }
}

void
ProgressFrame::SetProgressMessage(const wxString &mes)
{
	m_messageText->SetLabel(mes);
	if (m_value < 0)
		m_progressGauge->Pulse();
  UpdateAndYield();
}

void
ProgressFrame::SetProgressValue(double value)
{
	m_value = value;
	if (value < 0)
		m_progressGauge->Pulse();
	else
		m_progressGauge->SetValue((int)(value * 10000));
  UpdateAndYield();
}

void
ProgressFrame::SetInterruptValue(int value)
{
	m_interruptValue = value;
}

void
ProgressFrame::UpdateAndYield()
{
  Update();
  wxEventLoopBase * const loop = wxEventLoopBase::GetActive();
  if (loop != NULL)
    loop->YieldFor(wxEVT_CATEGORY_UI);
}

int
ProgressFrame::CheckInterrupt()
{
  UpdateAndYield();
  if (m_interruptValue)
    return m_interruptValue;
  if (::wxGetKeyState(WXK_ESCAPE))
    return 1;
  else return 0;

#if 0
	if (m_interruptValue) {
		int save = m_interruptValue;
		m_interruptValue = 0;
		return save;
	}

#if 1
	wxEventLoopBase * const loop = wxEventLoopBase::GetActive();
	if (loop != NULL)
		loop->YieldFor(wxEVT_CATEGORY_UI);
#else
#if __WXMAC__
	::wxYield();
#else
	{
		wxWindow *activeWin = NULL;
		if (this != NULL)
			activeWin = this;
		// else
		//	activeWin = GetMainFrame()->GetActiveChild();
		::wxSafeYield(activeWin);
	}
#endif
#endif
	if (::wxGetKeyState(WXK_ESCAPE))
		return 1;
	else return 0;
#endif
}

void
ProgressFrame::OnButtonPressed(wxCommandEvent& event)
{
  m_interruptValue = 1;
}
