/*
 *  MyDocManager.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 09/11/21.
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

#include "wx/filedlg.h"
#include "MyDocManager.h"
#include "MyMBConv.h"
#include "MyApp.h"
#include "../MolLib/Ruby_bind/Molby_extern.h"

BEGIN_EVENT_TABLE(MyDocManager, wxDocManager)
EVT_MENU(wxID_OPEN, MyDocManager::OnFileOpen)
EVT_MENU(wxID_SAVE, MyDocManager::OnFileSave)
EVT_MENU(wxID_SAVEAS, MyDocManager::OnFileSaveAs)
EVT_MENU(wxID_REVERT, MyDocManager::OnFileRevert)
END_EVENT_TABLE()

static const char *sReadOnlyTypes[] = {
	"out", "fchk", "log", "dat", "ins", "res", NULL
};

void
MyDocManager::SetDocumentTypesEnabled(const char **extensions, bool flag)
{
	wxList &tlist = GetTemplates();
	wxList::iterator iter;
	for (iter = tlist.begin(); iter != tlist.end(); ++iter) {
		int i;
		wxDocTemplate *dt = (wxDocTemplate *)(*iter);
		const char *p = (const char *)(dt->GetDefaultExtension().mb_str(WX_DEFAULT_CONV));
		for (i = 0; extensions[i] != NULL; i++) {
			if (strcmp(extensions[i], p) == 0) {
				dt->SetFlags(flag ? wxTEMPLATE_VISIBLE : wxTEMPLATE_INVISIBLE);
				break;
			}
		}
	}
}

void
MyDocManager::OnFileSave(wxCommandEvent& event)
{
	SetDocumentTypesEnabled(sReadOnlyTypes, false);
	wxDocManager::OnFileSave(event);
	SetDocumentTypesEnabled(sReadOnlyTypes, true);
}

void
MyDocManager::OnFileSaveAs(wxCommandEvent& event)
{
	SetDocumentTypesEnabled(sReadOnlyTypes, false);
	wxDocManager::OnFileSaveAs(event);
	SetDocumentTypesEnabled(sReadOnlyTypes, true);
}

bool
MyDocManager::GetDocumentDescriptionAtIndex(int idx, wxString *outDescription, wxString *outFilter, wxString *outExtension)
{
	wxList &tlist = GetTemplates();
	if (idx >= 0 && idx < tlist.GetCount()) {
		wxList::iterator iter;
		wxDocTemplate *dt = (wxDocTemplate *)(tlist.Item(idx)->GetData());
		if (outDescription != NULL)
			*outDescription = dt->GetDescription();
		if (outFilter != NULL)
			*outFilter = dt->GetFileFilter();
		if (outExtension != NULL)
			*outExtension = dt->GetDefaultExtension();
		return true;
	} else return false;
}

