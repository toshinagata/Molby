/*
 *  MyDocManager.h
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

#ifndef __MyDocManager_h__
#define __MyDocManager_h__

#include "wx/docview.h"

class MyDocManager: public wxDocManager {
public:
	void OnFileSave(wxCommandEvent& event);
	void OnFileSaveAs(wxCommandEvent& event);
	void SetDocumentTypesEnabled(const char **extensions, bool flag);
	bool GetDocumentDescriptionAtIndex(int idx, wxString *outDescription, wxString *outFilter, wxString *outExtension);
private:
	DECLARE_EVENT_TABLE()
};

#endif  /* __MyDocManager_h__ */
