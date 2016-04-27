/*
 *  MyDocument.h
 *  Molby
 *
 *  Created by Toshi Nagata on 08/10/24.
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

#ifndef __MyDocument_h__
#define __MyDocument_h__

#include "wx/docview.h"
#include "wx/cmdproc.h"
#include "wx/process.h"

#include "../MolLib/MolLib.h"

#include <stdio.h>   /*  For FILE structure  */

class wxProcess;  /*  For running QChem  */

/*  Custom Event for Document handling  */
extern const wxEventType MyDocumentEvent;
enum {
	MyDocumentEvent_willNeedCleanUndoStack = 1,
	MyDocumentEvent_documentModified,
	MyDocumentEvent_scriptMenuModified,
	MyDocumentEvent_updateDisplay,
	MyDocumentEvent_insertFrameFromMD,
	MyDocumentEvent_threadTerminated,
	MyDocumentEvent_openFilesByIPC,
	MyDocumentEvent_documentWillClose,
	MyDocumentEvent_openAuxiliaryDocuments
};

class MyDocument: public wxDocument
{
public:
	Molecule *mol;

	MyDocument();
	virtual ~MyDocument();

	void SetMolecule(Molecule *mol);
	Molecule *GetMolecule() { return mol; }

	//  For "register-undo" based undo/redo mechanism
	bool isUndoEnabled;
	bool isUndoing;
	bool isModifyNotificationSent;

	wxCommand *currentCommand;
	MolAction **undoStack;
	int countUndoStack;

	int undoGroupLevel;
	bool isCleanUndoStackRequested;

	bool hasFile;  /*  wxWidgets does not maintain this info for us  */

	int subThreadKind;  /*  0: none, 1: MM/MD (mutex is in Molecule structure), 2: QChem (also mutex is in Molecule structure)  */
	
	wxProcess *subProcess;  /*  subprocess object for QChem run  */
	long    subProcessPID;
	int		(*endSubProcessCallback)(Molecule *mol, int status);
	int		(*timerSubProcessCallback)(Molecule *mol, int timerCount);
	FILE    *subProcessStdout;  /*  Stdout of the subprocess  */
	FILE    *subProcessStderr;  /*  Stderr of the subprocess  */
	
	MainView	*GetMainView() { return (mol != NULL ? mol->mview : NULL); }
	void    SetIsUndoing(bool flag) { isUndoing = flag; }
	void    SetCurrentCommand(wxCommand *command) { currentCommand = command; }
	void    PushUndoAction(MolAction *action);
	void    CleanUndoStack(bool shouldRegister = true);
	bool	IsUndoEnabled() { return isUndoEnabled; }
	void	SetUndoEnabled(bool flag);	
	void	UpdateModifyFlag();
	void	BeginUndoGrouping();
	void	EndUndoGrouping();

	virtual bool Close();
	
	void    OnNeedCleanUndoStack(wxCommandEvent& event);

	void    OnDocumentModified(wxCommandEvent& event);
	void	OnImport(wxCommandEvent& event);
	void	OnExport(wxCommandEvent& event);
	void	OnExportGraphic(wxCommandEvent &event);

	void	OnCustomClose(wxCommandEvent &event);
	
	void	OnCopy(wxCommandEvent& event);
	void	OnCut(wxCommandEvent& event);
	void	OnPaste(wxCommandEvent& event);
	void	OnDelete(wxCommandEvent& event);

	void	OnCreateNewAtom(wxCommandEvent &event);
	void	OnCreateNewParameter(wxCommandEvent &event);
	void	OnCreatePiAnchor(wxCommandEvent &event);
	
	void	OnSelectAll(wxCommandEvent& event);
	void	OnSelectFragment(wxCommandEvent& event);
	void	OnSelectReverse(wxCommandEvent& event);

	void	OnAddHydrogen(wxCommandEvent& event);

	void	OnFitToScreen(wxCommandEvent& event);
	void	OnCenterSelection(wxCommandEvent& event);
//	void	OnShowMenu(wxCommandEvent& event);
//	void	OnToggleLineMode(wxCommandEvent &event);
//	void	OnShowGraphite(wxCommandEvent &event);

	void	OnShowAllAtoms(wxCommandEvent &event);
	void	OnHideSelected(wxCommandEvent &event);
	void	OnHideUnselected(wxCommandEvent &event);
	void	OnHideReverse(wxCommandEvent &event);
	
	void	DoMDOrMinimize(int minimize);
	void	OnMolecularDynamics(wxCommandEvent &event);
	void	OnMinimize(wxCommandEvent &event);
	void	OnStopMDRun(wxCommandEvent &event);
	
	long	RunSubProcess(const char *cmd, int (*callback)(Molecule *, int), int (*timerCallback)(Molecule *, int), FILE *output, FILE *errout);
	void	OnEndSubProcess(wxProcessEvent &event);
	void	FlushSubProcessOutput();

//	void	OnDefinePeriodicBox(wxCommandEvent &event);
//	void	OnShowPeriodicImage(wxCommandEvent &event);
//	void	OnPressureControl(wxCommandEvent &event);
//	void	OnDefineSymmetry(wxCommandEvent &event);
//	void	OnExpandBySymmetry(wxCommandEvent &event);

//	void	OnGuessUFFParameters(wxCommandEvent &event);
//	void	OnInvokeResp(wxCommandEvent &event);
//	void    OnInvokeAntechamber(wxCommandEvent &event);
//	void	OnCreateSanderInput(wxCommandEvent &event);
//	void	OnImportAmberFrcmod(wxCommandEvent &event);

//	void	OnCreateGamessInput(wxCommandEvent &event);
//	void	OnCreateMOPACInput(wxCommandEvent &event);
//	void	OnCreateMOCube(wxCommandEvent &event);
	
	void	OnInsertFrameFromMD(wxCommandEvent &event);
	void	OnUpdateDisplay(wxCommandEvent &event);
	void	OnSubThreadTerminated(wxCommandEvent &event);
	void	OnOpenAuxiliaryDocuments(wxCommandEvent &event);
	
	void	OnUpdateUI(wxUpdateUIEvent &event);

	void	TimerCallback(int timerCount);

 protected:
	virtual bool DoSaveDocument(const wxString& file);
	virtual bool DoOpenDocument(const wxString& file);
	virtual bool OnCreate(const wxString& path, long flags);
	virtual bool Revert();

 private:
	DECLARE_DYNAMIC_CLASS(MyDocument)
	DECLARE_EVENT_TABLE()
};

MyDocument *MyDocumentFromMolecule(Molecule *mp);

#endif
