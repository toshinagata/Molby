/*
 *  MyCommand.h
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

#ifndef __MyCommand_h__
#define __MyCommand_h__

#include "wx/cmdproc.h"

#include "../MolLib/MolLib.h"

class MyCommand: public wxCommand
{
public:
	Molecule *mol;
	MolAction **undoActions;
	int numUndoActions;
	MolAction **redoActions;
	int numRedoActions;

    MyCommand(Molecule *aMolecule, const wxString& name = wxEmptyString);
    virtual ~MyCommand();

	void SetUndoActions(MolAction **actions, int count);
	void SetRedoActions(MolAction **actions, int count);
	
    virtual bool Do();
	virtual bool Undo();
};

#endif /* __MyCommand_h__ */
