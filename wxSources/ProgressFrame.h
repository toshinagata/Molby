/*
 *  ProgressFrame.h
 *  Molby
 *
 *  Created by Toshi Nagata on 09/07/15.
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

#ifndef __ProgressFrame_h__
#define __ProgressFrame_h__

#include "wx/frame.h"
#include "wx/gauge.h"
#include "wx/utils.h"  //  For wxWindowDisabler

class wxStaticText;
//class wxGauge;    //  This forward declaration does not work because in wxMSW wxGauge is #define'ed as wxGauge95

class ProgressFrame: public wxFrame
{

public:
	ProgressFrame(const wxString& title, const wxString &mes);
	virtual ~ProgressFrame();

  static ProgressFrame *FindProgressFrameWithID(int id);

	void SetProgressMessage(const wxString &mes);
	void SetProgressValue(double value);
	void SetInterruptValue(int value);
	int CheckInterrupt();
  void UpdateAndYield();
  void OnButtonPressed(wxCommandEvent& event);

	wxStaticText *m_messageText;
	wxGauge *m_progressGauge;
  wxButton *m_cancelButton;
	double m_value;
	int m_interruptValue;
  wxWindowDisabler *m_disabler;
  
  static int c_uniqueID;
  static wxWindowList c_progressFrames;

private:
	DECLARE_EVENT_TABLE()
};

#endif /* __ProgressFrame_h__ */
