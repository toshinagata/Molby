/*
 *  MyCommand.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 08/10/25.
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

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !wxUSE_DOC_VIEW_ARCHITECTURE
#error "You should have DocView architecture enabled in your wxWidgets installation."
#endif

#include "MoleculeView.h"
#include "MyCommand.h"
#include "MyDocument.h"

MyCommand::MyCommand(Molecule *aMolecule, const wxString& name):
	wxCommand(true, name)
{
	mol = MoleculeRetain(aMolecule);
	undoActions = redoActions = NULL;
	numUndoActions = numRedoActions = 0;
}

MyCommand::~MyCommand()
{
	MoleculeRelease(mol);
	SetUndoActions(NULL, 0);
	SetRedoActions(NULL, 0);
}

void
MyCommand::SetUndoActions(MolAction **actions, int count)
{
	int i;
	if (undoActions != NULL) {
		for (i = 0; i < numUndoActions; i++)
			MolActionRelease(undoActions[i]);
		free(undoActions);
	}
	undoActions = actions;
	numUndoActions = count;
}

void
MyCommand::SetRedoActions(MolAction **actions, int count)
{
	int i;
	if (redoActions != NULL) {
		for (i = 0; i < numRedoActions; i++)
			MolActionRelease(redoActions[i]);
		free(redoActions);
	}
	redoActions = actions;
	numRedoActions = count;
}

bool
MyCommand::Do()
{
	MyDocument *doc = MyDocumentFromMolecule(mol);
	int i;
	if (doc != NULL) {
		bool retval = true;
		if (redoActions != NULL) {
			doc->SetCurrentCommand(this);
			for (i = numRedoActions - 1; i >= 0; i--) {
				if (MolActionPerform(mol, redoActions[i]) != 0) {
					retval = false;
					break;
				}
			}
			doc->CleanUndoStack(retval);
		}
		return retval;
	}
	return false;
}

bool
MyCommand::Undo()
{
	MyDocument *doc = MyDocumentFromMolecule(mol);
	int i;
	if (doc != NULL) {
		bool retval = true;
		if (undoActions != NULL) {
			doc->SetIsUndoing(true);
			doc->SetCurrentCommand(this);
			for (i = numUndoActions - 1; i >= 0; i--) {
				if (MolActionPerform(mol, undoActions[i]) != 0) {
					retval = false;
					break;
				}
			}
			doc->CleanUndoStack(retval);
		}
		return retval;
	}
	return false;
}
