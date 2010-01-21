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

#if __WXMAC__
#include <Carbon/Carbon.h>
#endif

BEGIN_EVENT_TABLE(ProgressFrame, wxFrame)
//    EVT_TEXT_ENTER(-1, ConsoleFrame::OnEnterPressed)
//	EVT_CHAR(ConsoleFrame::OnChar)
//	EVT_RICHTEXT_RETURN(-1, ConsoleFrame::OnEnterPressed)
END_EVENT_TABLE()

ProgressFrame::ProgressFrame(const wxString &title, const wxString &mes):
	wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxDefaultSize, wxCAPTION)
{
	//  Vertical sizer containing (1) message text, (2) progress gauge, (3) note text
	wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
	
	m_messageText = new wxStaticText(this, -1, wxT("Message"), wxDefaultPosition, wxSize(240, 40), wxST_NO_AUTORESIZE);
	sizer->Add(m_messageText, 0, wxALL | wxEXPAND, 10);   // Can expand horizontally
	
	m_progressGauge = new wxGauge(this, -1, 10000, wxDefaultPosition, wxSize(240, 24), wxGA_HORIZONTAL);
	sizer->Add(m_progressGauge, 0, wxALL | wxEXPAND, 10);
	
	wxStaticText *noteText = new wxStaticText(this, -1, wxT("Press ESC to interrupt"), wxDefaultPosition, wxSize(240, 20), wxALIGN_CENTRE | wxST_NO_AUTORESIZE);
/*	wxFont smallFont(noteText->GetFont());
	int size = smallFont.GetPointSize();
	if (size >= 14)
		size = 12;
	else if (size >= 12)
		size = 10;
	else if (size >= 10)
		size = 9;
	smallFont.SetPointSize(size);
	noteText->SetFont(smallFont); */
	noteText->SetFont(*wxSMALL_FONT);
	sizer->Add(noteText, 0, wxALL | wxEXPAND, 10);

	m_value = -1.0;
	m_progressGauge->Pulse();
	m_messageText->SetLabel(mes);
	
	m_interruptValue = 0;

	sizer->Layout();
	this->SetSizerAndFit(sizer);
	this->Centre();
	this->Show();
	
#if __WXMAC__
	::SetWindowModality(((WindowRef)MacGetWindowRef()), kWindowModalityAppModal, NULL);
#endif
}

ProgressFrame::~ProgressFrame()
{
}

void
ProgressFrame::SetProgressMessage(const wxString &mes)
{
	m_messageText->SetLabel(mes);
	if (m_value < 0)
		m_progressGauge->Pulse();
	Update();
}

void
ProgressFrame::SetProgressValue(double value)
{
	m_value = value;
	if (value < 0)
		m_progressGauge->Pulse();
	else
		m_progressGauge->SetValue((int)(value * 10000));
	Update();
}

void
ProgressFrame::SetInterruptValue(int value)
{
	m_interruptValue = value;
}

int
ProgressFrame::CheckInterrupt()
{
	if (this != NULL && m_interruptValue) {
		int save = m_interruptValue;
		m_interruptValue = 0;
		return save;
	}

#if __WXMAC__
	::wxYield();
#else
	{
		wxWindow *activeWin;
		if (this != NULL)
			activeWin = this;
		else
			activeWin = GetMainFrame()->GetActiveChild();
		::wxSafeYield(activeWin);
	}
#endif
	if (::wxGetKeyState(WXK_ESCAPE))
		return 1;
	else return 0;
}
