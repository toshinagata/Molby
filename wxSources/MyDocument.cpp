/*
 *  MyDocument.cpp
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

#if wxUSE_STD_IOSTREAM
    #include "wx/ioswrap.h"
#else
    #include "wx/txtstrm.h"
#endif

#include "wx/clipbrd.h"
#include "wx/filename.h"
#include "wx/dir.h"

#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include "MyApp.h"
#include "MyDocManager.h"
#include "MyDocument.h"
#include "MoleculeView.h"
#include "MyCommand.h"
#include "MyClipboardData.h"
#include "MyThread.h"
#include "MyMBConv.h"

#include "../MolLib/Ruby_bind/Molby_extern.h"
#include "../MolLib/MD/MDCore.h"
#include "../MolLib/Missing.h"

IMPLEMENT_DYNAMIC_CLASS(MyDocument, wxDocument)

const wxEventType MyDocumentEvent = wxNewEventType();

BEGIN_EVENT_TABLE(MyDocument, wxDocument)
	EVT_COMMAND(MyDocumentEvent_willNeedCleanUndoStack, MyDocumentEvent, MyDocument::OnNeedCleanUndoStack)
	EVT_COMMAND(MyDocumentEvent_documentModified, MyDocumentEvent, MyDocument::OnDocumentModified)
	EVT_COMMAND(MyDocumentEvent_insertFrameFromMD, MyDocumentEvent, MyDocument::OnInsertFrameFromMD)
	EVT_COMMAND(MyDocumentEvent_updateDisplay, MyDocumentEvent, MyDocument::OnUpdateDisplay)
	EVT_COMMAND(MyDocumentEvent_threadTerminated, MyDocumentEvent, MyDocument::OnSubThreadTerminated)
	EVT_COMMAND(MyDocumentEvent_openAuxiliaryDocuments, MyDocumentEvent, MyDocument::OnOpenAuxiliaryDocuments)
	EVT_MENU(myMenuID_Import, MyDocument::OnImport)
	EVT_MENU(myMenuID_Export, MyDocument::OnExport)
	EVT_MENU(myMenuID_ExportGraphic, MyDocument::OnExportGraphic)
	EVT_MENU(wxID_COPY, MyDocument::OnCopy)
	EVT_MENU(wxID_PASTE, MyDocument::OnPaste)
	EVT_MENU(wxID_CUT, MyDocument::OnCut)
	EVT_MENU(wxID_DELETE, MyDocument::OnDelete)
//	EVT_MENU(wxID_CLOSE, MyDocument::OnCustomClose)
	EVT_MENU(myMenuID_CreateNewAtom, MyDocument::OnCreateNewAtom)
	EVT_MENU_RANGE(myMenuID_CreateNewVdwParameter, myMenuID_CreateNewVdwCutoffParameter, MyDocument::OnCreateNewParameter)
	EVT_MENU(myMenuID_CreatePiAnchor, MyDocument::OnCreatePiAnchor)
	EVT_MENU(wxID_SELECTALL, MyDocument::OnSelectAll)
	EVT_MENU(myMenuID_SelectFragment, MyDocument::OnSelectFragment)
	EVT_MENU(myMenuID_SelectReverse, MyDocument::OnSelectReverse)
	EVT_MENU(myMenuID_FitToScreen, MyDocument::OnFitToScreen)
	EVT_MENU(myMenuID_CenterSelection, MyDocument::OnCenterSelection)
//	EVT_MENU(myMenuID_ShowUnitCell, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowPeriodicBox, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowHydrogens, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowDummyAtoms, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowExpandedAtoms, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowEllipsoids, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowRotationCenter, MyDocument::OnShowMenu)
//	EVT_MENU(myMenuID_ShowGraphite, MyDocument::OnShowGraphite)
//	EVT_MENU(myMenuID_LineMode, MyDocument::OnToggleLineMode)
	EVT_MENU_RANGE(myMenuID_AddHydrogenSp3, myMenuID_AddHydrogenBent, MyDocument::OnAddHydrogen)
	EVT_UPDATE_UI_RANGE(myMenuID_MyFirstMenuItem, myMenuID_MyLastMenuItem, MyDocument::OnUpdateUI)
	EVT_MENU(myMenuID_MolecularDynamics, MyDocument::OnMolecularDynamics)
	EVT_MENU(myMenuID_Minimize, MyDocument::OnMinimize)
	EVT_MENU(myMenuID_StopMDRun, MyDocument::OnStopMDRun)
	EVT_MENU(myMenuID_ShowAllAtoms, MyDocument::OnShowAllAtoms)
	EVT_MENU(myMenuID_HideReverse, MyDocument::OnHideReverse)
	EVT_MENU(myMenuID_HideSelected, MyDocument::OnHideSelected)
	EVT_MENU(myMenuID_HideUnselected, MyDocument::OnHideUnselected)
	EVT_END_PROCESS(-1, MyDocument::OnEndSubProcess)
END_EVENT_TABLE()

MyDocument::MyDocument()
{
	mol = MoleculeNew();
	isUndoing = false;
	isUndoEnabled = true;
	isModifyNotificationSent = false;
	currentCommand = NULL;
	undoStack = NULL;
	countUndoStack = 0;
	undoGroupLevel = 0;
	isCleanUndoStackRequested = false;
	hasFile = false;
	subThreadKind = 0;
	subProcess = NULL;
	endSubProcessCallback = NULL;
	timerSubProcessCallback = NULL;
}

MyDocument::~MyDocument()
{
	int i;
	Molecule *mol2 = mol;
	mol = NULL;

	if (subProcess != NULL) {
		subProcess->Detach();
		subProcess->Kill(subProcess->GetPid(), wxSIGTERM, wxKILL_CHILDREN);
	}
	
	/*  May be unnecessary?  */
	MoleculeView *view = (MoleculeView *)GetFirstView();
	if (view != NULL) {
		view->OnMoleculeReplaced();
	}

	if (mol2 != NULL)
		MoleculeRelease(mol2);
	if (undoStack != NULL) {
		for (i = 0; i < countUndoStack; i++)
			MolActionRelease(undoStack[i]);
		free(undoStack);
	}
	
	wxGetApp().DisableTimerForDocument(this);
}

void
MyDocument::SetMolecule(Molecule *aMolecule)
{
	Molecule *mol2 = mol;
	if (mol == aMolecule)
		return;
	mol = aMolecule;
	if (aMolecule != NULL)
		MoleculeRetain(aMolecule);

	MoleculeView *view = (MoleculeView *)GetFirstView();
	if (view != NULL) {
		view->OnMoleculeReplaced();
	}
	if (mol2 != NULL)
		MoleculeRelease(mol2);
}

bool
MyDocument::DoSaveDocument(const wxString& file)
{
	char *buf = NULL;
	char *p = strdup((const char *)file.mb_str(wxConvFile));
	size_t len = strlen(p);
	int retval;
	if (MolActionCreateAndPerform(mol, SCRIPT_ACTION("s"), "molsave", p) != 0) {
		free(p);
		return false;
	}
	retval = 0;
	MoleculeLock(mol);
	if (len > 4 && strcasecmp(p + len - 4, ".psf") == 0) {
		/*  Write as a psf and a pdb file  */
		char *pp = (char *)malloc(len + 2);
		strcpy(pp, p);
		strcpy(pp + len - 4, ".pdb");
		retval = MoleculeWriteToPdbFile(mol, pp, &buf);
		if (retval != 0) {
			free(pp);
			goto exit;
		}
		if (mol->cell != NULL) {
			/*  Write an extended info (bounding box)  */
			strcpy(pp + len - 4, ".info");
			retval = MoleculeWriteExtendedInfo(mol, pp, &buf);
			if (retval != 0) {
				free(pp);
				goto exit;
			}
		}
	}
	GetCommandProcessor()->MarkAsSaved();
	hasFile = true;
	MoleculeSetPath(mol, p);
exit:
	free(p);
	MoleculeUnlock(mol);
	return (retval == 0);
}

bool
MyDocument::DoOpenDocument(const wxString& file)
{
	char *p;
	int len;
	Molecule *newmol;
	p = strdup((const char *)file.mb_str(wxConvFile));
	newmol = MoleculeNew();
	SetMolecule(newmol);
	MoleculeRelease(newmol);
	SetUndoEnabled(false);
	if (MolActionCreateAndPerform(newmol, SCRIPT_ACTION("s"), "molload", p) != 0) {
		free(p);
		SetMolecule(NULL);
		SetUndoEnabled(true);
		return false;
	}
	
	/*  Does this document have multiple representation of molecules?  */
	if (MolActionCreateAndPerform(newmol, SCRIPT_ACTION(";i"), "lambda { @aux_mols ? @aux_mols.count : 0 }", &len) == 0 && len > 0) {
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_openAuxiliaryDocuments);
		wxPostEvent(this, myEvent);
	}
	
	if ((len = strlen(p)) > 4 && strcasecmp(p + len - 4, ".psf") == 0) {
		//  Look for a ".pdb" file with the same basename 
		char *buf;
		strcpy(p + len - 4, ".pdb");
		//  The error will be ignored
		MoleculeReadCoordinatesFromPdbFile(newmol, p, &buf);
		//  Look for an ".info" file with the same basename
		p = (char *)realloc(p, len + 2);
		strcpy(p + len - 4, ".info");
		MoleculeReadExtendedInfo(newmol, p, &buf);
		free(buf);
	}
	free(p);
	Modify(false);
	GetCommandProcessor()->MarkAsSaved();
	hasFile = true;
	if (newmol->natoms > 1000)
		newmol->mview->lineMode = 1;
	if (TrackballGetModifyCount(newmol->mview->track) == 0)
		MainView_resizeToFit(newmol->mview);
	MoleculeCallback_notifyModification(newmol, 0);
	SetUndoEnabled(true);
	return true;
}

void
MyDocument::OnOpenAuxiliaryDocuments(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""),
									"lambda {\n"
									"  fn = self.name\n"
									"  @aux_mols.each_with_index { |am, i| \n"
									"    m = Molecule.open; m.set_molecule(am)\n"
									"    m.set_name(fn + \"[#{i + 2}]\")\n"
									"}; @aux_mols = nil }");
}

bool
MyDocument::Revert()
{
	if (wxDocument::Revert()) {
		MainViewCallback_selectTable(mol->mview, 0);
		return true;
	} else return false;
}

/*  Override to intercept view creation for running script  */
bool
MyDocument::OnCreate(const wxString& path, long flags)
{
	if (path.EndsWith(wxT(".rb")) || path.EndsWith(wxT(".mrb"))) {
		wxGetApp().OnOpenFiles(path);
		return false;  /*  This document will be deleted  */
	} else {
		return wxDocument::OnCreate(path, flags);
	}
}

void
MyDocument::OnImport(wxCommandEvent& event)
{
	wxString wildcard;
	{
		wxString desc, filter, ext;
		int i;
		/*  File filter is built from MyDocManager information  */
		MyDocManager *docm = wxGetApp().DocManager();
		for (i = 0; docm->GetDocumentDescriptionAtIndex(i, &desc, &filter, &ext); i++) {
			if (filter.Contains(_T("*.*"))) {
				i = -1;
				break;
			}
			if (wildcard != _T("")) {
				wildcard += (_T("|"));
			}
			wildcard += (desc + _T(" (") + filter + _T(")|") + filter);
		}
		/*  Insert Import-only file types before "All files"  */
		wildcard += _T("|AMBER mdcrd file (*.crd;*.mdcrd)|*.crd;*.mdcrd");
		wildcard += _T("|DCD file (*.dcd)|*.dcd");
		if (i == -1)
			wildcard += (_T("|") + desc + _T(" (") + filter + _T(")|") + filter);
	}

	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Choose Coordinate File"), _T(""), _T(""), wildcard, wxFD_OPEN | wxFD_CHANGE_DIR | wxFD_FILE_MUST_EXIST);
	if (dialog->ShowModal() == wxID_OK) {
		char *p = strdup((const char *)(dialog->GetPath().mb_str(wxConvFile)));
		MoleculeLock(mol);
		MolActionCreateAndPerform(mol, SCRIPT_ACTION("s"), "molload", p);
		if (gLoadSaveErrorMessage != NULL)
			MyAppCallback_showScriptMessage("On loading %s:\n%s\n", p, gLoadSaveErrorMessage);
		MoleculeUnlock(mol);
		free(p);
	}
	dialog->Destroy();
}

void
MyDocument::OnExport(wxCommandEvent& event)
{
	wxString wildcard;
	wxFileName fname(GetFilename());
    wxString fnstr = GetUserReadableName();
	{
		/*  File filter is built from MyDocManager information  */
		wxString desc, filter, ext;
		int i;
		MyDocManager *docm = wxGetApp().DocManager();
		if ((i = fnstr.Find('.', true)) != wxNOT_FOUND) {
			fnstr = fnstr.Mid(0, i);
		}
		for (i = 0; docm->GetDocumentDescriptionAtIndex(i, &desc, &filter, &ext); i++) {
			if (ext == _T("mbsf") || ext == _T("out") || ext == _T("log") || ext == _T("fchk"))
				continue;
			if (filter.Contains(_T("*.*"))) {
				i = -1;
				break;
			}
			if (wildcard != _T("")) {
				wildcard += (_T("|"));
			}
			wildcard += (desc + _T(" (") + filter + _T(")|") + filter);
		}
		wildcard += _T("|AMBER mdcrd file (*.crd;*.mdcrd)|*.crd;*.mdcrd");
		wildcard += _T("|DCD file (*.dcd)|*.dcd");
		if (i == -1)
			wildcard += (_T("|") + desc + _T(" (") + filter + _T(")|") + filter);
	}
	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Export coordinates"), fname.GetPath(), fnstr + _T(".psf"), wildcard, wxFD_SAVE | wxFD_OVERWRITE_PROMPT | wxFD_CHANGE_DIR);
	if (dialog->ShowModal() == wxID_OK) {
		char *p = strdup((const char *)(dialog->GetPath().mb_str(wxConvFile)));
		MoleculeLock(mol);
		MolActionCallback_setUndoRegistrationEnabled(mol, 0);
		MolActionCreateAndPerform(mol, SCRIPT_ACTION("s"), "molsave", p);
		MolActionCallback_setUndoRegistrationEnabled(mol, 1);
		MoleculeUnlock(mol);
		free(p);
	}
	dialog->Destroy();
}

void
MyDocument::OnExportGraphic(wxCommandEvent& event)
{
	wxString wildcard = _T("PNG File (*.png)|*.png|TIFF File (*.tif)|*.tif|All Files (*.*)|*.*");
	wxFileName fname(GetFilename());
	wxString fnstr;
	Int scale, bg_color, n, i;
	i = MolActionCreateAndPerform(mol, SCRIPT_ACTION(";i"), "ask_graphic_export_scale", &n);
	if (i != 0 || n < 0)
		return;
	i = MyAppCallback_getGlobalSettingsWithType("global.export_graphic_scale", 'i', &scale);
	if (i != 0)
		scale = 4;
	i = MyAppCallback_getGlobalSettingsWithType("global.export_background_color", 'i', &bg_color);
	if (i != 0)
		bg_color = 0;
	fnstr = GetUserReadableName();
	if ((i = fnstr.Find('.', true)) != wxNOT_FOUND) {
		fnstr = fnstr.Mid(0, i);
	}
	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Export Graphic"), fname.GetPath(), fnstr + _T(".png"), wildcard, wxFD_SAVE | wxFD_OVERWRITE_PROMPT | wxFD_CHANGE_DIR);
	if (dialog->ShowModal() == wxID_OK) {
		wxString fnpath = dialog->GetPath();
		MoleculeView *myview = (MoleculeView *)GetFirstView();
		myview->DoExportGraphic(fnpath, scale, bg_color, 0, 0);
	}
	dialog->Destroy();
}

void
MyDocument::SetUndoEnabled(bool flag)
{
	if (flag) {
		isUndoEnabled = true;
	} else {
		//  Remove all registered actions
		wxCommandProcessor *cmdProc = GetCommandProcessor();
		currentCommand = NULL;
		cmdProc->ClearCommands();
		CleanUndoStack(false);
		isUndoEnabled = false;
		//  TODO: mark the document as "edited"
	}
}

void
MyDocument::PushUndoAction(MolAction *action)
{
	if (countUndoStack % 8 == 0) {
		if (undoStack == NULL)
			undoStack = (MolAction **)malloc(sizeof(MolAction *) * 8);
		else
			undoStack = (MolAction **)realloc(undoStack, sizeof(MolAction *) * (countUndoStack + 8));
		if (undoStack == NULL)
			return;
	}
	undoStack[countUndoStack++] = action;
	MolActionRetain(action);
	if (countUndoStack == 1) {
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_willNeedCleanUndoStack);
		wxPostEvent(this, myEvent);
	}
}

/*  Update the modify flag to match the GetCommandProcessor isDirty flag
    (Is this really necessary? It should be handled by wxDocument automatically.  */
void
MyDocument::UpdateModifyFlag()
{
	Modify(GetCommandProcessor()->IsDirty());
}

void
MyDocument::BeginUndoGrouping()
{
	undoGroupLevel++;
}

void
MyDocument::EndUndoGrouping()
{
	if (undoGroupLevel <= 0)
		return;  /* This should not happen  */
	if (--undoGroupLevel == 0) {
		if (isCleanUndoStackRequested) {
			/*  Resend the event so that it can be processed at the next idle time  */
			wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_willNeedCleanUndoStack);
			wxPostEvent(this, myEvent);
			isCleanUndoStackRequested = false;
		}
	}
}

void
MyDocument::CleanUndoStack(bool shouldRegister)
{
	if (undoStack != NULL) {
		if (shouldRegister) {
			MyCommand *cmd = (MyCommand *)currentCommand;
			if (cmd == NULL)
				cmd = new MyCommand(mol, _T(" "));
			if (isUndoing)
				cmd->SetRedoActions(undoStack, countUndoStack);
			else
				cmd->SetUndoActions(undoStack, countUndoStack);
			if (currentCommand == NULL) {
				if (!GetCommandProcessor()->Submit(cmd))
					delete cmd;
				UpdateModifyFlag();
			}
		} else {
			int i;
			for (i = 0; i < countUndoStack; i++)
				MolActionRelease(undoStack[i]);
			free(undoStack);
		}
	}
	isUndoing = false;
	undoStack = NULL;
	countUndoStack = 0;
	currentCommand = NULL;
}

void
MyDocument::OnCustomClose(wxCommandEvent &event)
{
//	RubyValue val;
//	MolActionCreateAndPerform(mol, SCRIPT_ACTION(";r"), "close_all_auxiliary_windows", &val);
//	if (val == NULL || val == RubyNil)
//		event.Skip();
}

bool
MyDocument::Close()
{
	if (mol != NULL && mol->mutex != NULL) {
		const char *msg;
		if (subThreadKind == 1)
			msg = "MM/MD";
		else msg = "Some background process";
		MyAppCallback_errorMessageBox("%s is running: please stop it before closing", msg);
		return false;
	}
	if (wxDocument::Close()) {
		/*  Call close hander in the Ruby world  */
		MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "on_close");
		/*  Send a message that this document will close  */
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_documentWillClose);
		myEvent.SetEventObject(this);
		ProcessEvent(myEvent);
		return true;
	} else return false;		
}

void
MyDocument::OnNeedCleanUndoStack(wxCommandEvent& event)
{
	if (undoGroupLevel == 0)
		CleanUndoStack(true);
	else {
		/*  Do not respond to this event immediately; the same event will be
		    resent when undoGroupLevel becomes 0. See EndUndoGrouping(). */
		isCleanUndoStackRequested = true;
	}
}

void
MyDocument::OnDocumentModified(wxCommandEvent& event)
{
	isModifyNotificationSent = false;
	MoleculeClearModifyCount(GetMainView()->mol);

	/*  Call modified handler in the Ruby world  */
	/*  (Does not if undo is disabled --- e.g. during loading structure)  */
	if (isUndoEnabled)
		MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "on_modified");

	event.Skip();  //  Also pass to other notification handlers
	UpdateModifyFlag();
}

void
MyDocument::OnCopy(wxCommandEvent& event)
{
	wxWindow *focusWindow = wxWindow::FindFocus();
	if (focusWindow->IsKindOf(CLASSINFO(wxTextCtrl))) {
		event.Skip();
		return;
	}
	if (focusWindow == ((MoleculeView *)GetFirstView())->GetListCtrl() && GetMainView()->tableIndex == kMainViewParameterTableIndex) {
		MainView_copyOrCutParameters(GetMainView(), 2);
	} else {
		MoleculeLock(mol);
		MainView_copy(GetMainView());
		MoleculeUnlock(mol);
	}
}

void
MyDocument::OnCut(wxCommandEvent& event)
{
	wxWindow *focusWindow = wxWindow::FindFocus();
	if (focusWindow->IsKindOf(CLASSINFO(wxTextCtrl))) {
		event.Skip();
		return;
	}
	if (focusWindow == ((MoleculeView *)GetFirstView())->GetListCtrl() && GetMainView()->tableIndex == kMainViewParameterTableIndex) {
		MainView_copyOrCutParameters(GetMainView(), 3);
	} else {
		MoleculeLock(mol);
		MainView_cut(GetMainView());
		MoleculeUnlock(mol);
	}
}

void
MyDocument::OnPaste(wxCommandEvent& event)
{
	wxWindow *focusWindow = wxWindow::FindFocus();
	if (focusWindow->IsKindOf(CLASSINFO(wxTextCtrl))) {
		event.Skip();
		return;
	}
	if (focusWindow == ((MoleculeView *)GetFirstView())->GetListCtrl() && GetMainView()->tableIndex == kMainViewParameterTableIndex) {
		MainView_pasteParameters(GetMainView());
	} else {
		MoleculeLock(mol);
		MainView_paste(GetMainView());
		MoleculeUnlock(mol);
	}
}

void
MyDocument::OnDelete(wxCommandEvent& event)
{
	if (wxWindow::FindFocus() == ((MoleculeView *)GetFirstView())->GetListCtrl() && GetMainView()->tableIndex == kMainViewParameterTableIndex) {
		MainView_copyOrCutParameters(GetMainView(), 1);
	} else {
		MoleculeLock(mol);
		MainView_delete(GetMainView());
		MoleculeUnlock(mol);
	}
}

void
MyDocument::OnCreateNewAtom(wxCommandEvent &event)
{
	Int idx, i, j, row;
	char name[6];
	IntGroup *ig = MoleculeGetSelection(mol);
	MainView *mview = GetMainView();
	Atom *ap, arec;

	if (mview == NULL)
		return;

	/*  Make an atom name "Cxxx"  */
	for (i = 0; i < 1000; i++) {
		sprintf(name, "C%03d", i);
		for (j = 0, ap = mol->atoms; j < mol->natoms; j++, ap = ATOM_NEXT(ap)) {
			if (strncmp(ap->aname, name, 4) == 0)
				break;
		}
		if (j >= mol->natoms)
			break;
	}
    memset(&arec, 0, sizeof(arec));
    strncpy(arec.aname, name, 4);
	arec.type = AtomTypeEncodeToUInt("c3");
	arec.element[0] = 'C';
	arec.atomicNumber = 6;
	arec.weight = WeightForAtomicNumber(6);
	arec.occupancy = 1.0;
	
	if (ig != NULL && IntGroupGetCount(ig) > 0) {
		idx = IntGroupGetEndPoint(ig, IntGroupGetIntervalCount(ig) - 1);
	} else {
		idx = mol->natoms;
	}
	
	if (MolActionCreateAndPerform(mol, gMolActionAddAnAtom, &arec, idx, &idx) != 0)
		return;

	/*  Show the atom table and select the newly created atom  */
	MainViewCallback_selectTable(mview, kMainViewAtomTableIndex);
	ig = IntGroupNewWithPoints(idx, 1, -1);
	MoleculeSetSelection(mol, ig);
	IntGroupRelease(ig);
	MainView_refreshTable(mview);
	row = MainView_indexToTableRow(mview, idx);
	MainViewCallback_startEditText(mview, row, 1);
}

void
MyDocument::OnCreatePiAnchor(wxCommandEvent &event)
{
	Int idx, row;
	MainView *mview = GetMainView();
	IntGroup *ig = MoleculeGetSelection(mol), *ig2;
	if (ig == NULL || IntGroupGetCount(ig) < 2)
		return;  /*  Do nothing  */
	if (MolActionCreateAndPerform(mol, SCRIPT_ACTION("G;i"),
			"lambda { |g| create_pi_anchor('AN', atom_group(g) { |ap| ap.atomic_number != 1 }, nil, nil, g.max + 1).index rescue -1 }",
			ig, &idx) != 0)
		return;
	MainViewCallback_selectTable(mview, kMainViewAtomTableIndex);
	ig2 = IntGroupNewWithPoints(idx, 1, -1);
	MoleculeSetSelection(mol, ig2);
	IntGroupRelease(ig2);
	MainView_refreshTable(mview);
	row = MainView_indexToTableRow(mview, idx);
	MainViewCallback_ensureVisible(mview, row);
}

void
MyDocument::OnCreateNewParameter(wxCommandEvent &event)
{
	int uid = event.GetId();
	Int parType, n;
	UnionPar ubuf;
	IntGroup *ig;
	UInt ctype = AtomTypeEncodeToUInt("C");
	Double cweight = WeightForAtomicNumber(6);
	memset(&ubuf, 0, sizeof(ubuf));
	ubuf.bond.src = -1;  /*  Undefined  */
	switch (uid) {
		case myMenuID_CreateNewVdwParameter:
			parType = kVdwParType;
			ubuf.vdw.type1 = ctype;
			ubuf.vdw.atomicNumber = 6;
			ubuf.vdw.weight = cweight;
			break;
		case myMenuID_CreateNewBondParameter:
			parType = kBondParType;
			ubuf.bond.type1 = ubuf.bond.type2 = ctype;
			break;
		case myMenuID_CreateNewAngleParameter:
			parType = kAngleParType;
			ubuf.angle.type1 = ubuf.angle.type2 = ubuf.angle.type3 = ctype;
			break;
		case myMenuID_CreateNewDihedralParameter:
			parType = kDihedralParType;
			ubuf.torsion.type1 = ubuf.torsion.type2 = ubuf.torsion.type3 = ubuf.torsion.type4 = ctype;
			break;
		case myMenuID_CreateNewImproperParameter:
			parType = kImproperParType;
			ubuf.torsion.type1 = ubuf.torsion.type2 = ubuf.torsion.type3 = ubuf.torsion.type4 = ctype;
			break;
		case myMenuID_CreateNewVdwPairParameter:
			parType = kVdwPairParType;
			ubuf.vdwp.type1 = ubuf.vdwp.type2 = ctype;
			break;
		case myMenuID_CreateNewVdwCutoffParameter:
			parType = kVdwCutoffParType;
			ubuf.vdwcutoff.type1 = ubuf.vdwcutoff.type2 = ctype;
			break;			
		default:
			return;
	}
	if (mol->par == NULL) {
		char *errmsg;
		if (MoleculePrepareMDArena(mol, 1, &errmsg) < 0) {
			MyAppCallback_messageBox(errmsg, "MM/MD Setup Error", 1, 3);
			free(errmsg);
			return;
		}
	}
	n = ParameterGetCountForType(mol->par, parType);
	ig = IntGroupNewWithPoints(n, 1, -1);
	MolActionCreateAndPerform(mol, gMolActionAddParameters, parType, ig, 1, &ubuf);
	if (ParameterGetCountForType(mol->par, parType) == n + 1) {
		/*  Successful creation of the parameter  */
		MainView *mview = GetMainView();
		Int row;
		MainViewCallback_selectTable(mview, kMainViewParameterTableIndex);
		MainView_refreshTable(mview);
		row = ParameterTableGetRowFromTypeAndIndex(mol->par, parType, n);
		MainViewCallback_startEditText(mview, row, 1);		
	}
}

void
MyDocument::OnSelectAll(wxCommandEvent& event)
{
	if (wxWindow::FindFocus() == ((MoleculeView *)GetFirstView())->GetListCtrl() && mol->mview->tableIndex == kMainViewParameterTableIndex) {
	} else {
		MoleculeLock(mol);
		MainView_selectAll(GetMainView());
		MoleculeUnlock(mol);
	}
}

void
MyDocument::OnSelectFragment(wxCommandEvent& event)
{
	MoleculeLock(mol);
	MainView_selectFragment(GetMainView());
	MoleculeUnlock(mol);
}

void
MyDocument::OnSelectReverse(wxCommandEvent& event)
{
	MoleculeLock(mol);
	MainView_selectReverse(GetMainView());
	MoleculeUnlock(mol);
}

void
MyDocument::OnAddHydrogen(wxCommandEvent& event)
{
	int uid = event.GetId();
	const char *type;
	IntGroup *ig;
	switch (uid) {
		case myMenuID_AddHydrogenSp3: type = "td"; break;
		case myMenuID_AddHydrogenSp2: type = "tr"; break;
		case myMenuID_AddHydrogenLinear: type = "li"; break;
		case myMenuID_AddHydrogenPyramidal: type = "py"; break;
		case myMenuID_AddHydrogenBent: type = "be"; break;
		default: return;
	}
	MoleculeLock(mol);
	ig = MoleculeGetSelection(mol);
	MolActionCreateAndPerform(mol, SCRIPT_ACTION("Gs"), "add_hydrogen_on_group", ig, type);
	MoleculeUnlock(mol);
}

void
MyDocument::OnFitToScreen(wxCommandEvent& event)
{
	MoleculeLock(mol);
	MainView_resizeToFit(GetMainView());
	MoleculeUnlock(mol);
}

void
MyDocument::OnCenterSelection(wxCommandEvent& event)
{
	MoleculeLock(mol);
	MainView_centerSelection(GetMainView());
	MoleculeUnlock(mol);
}

/*
void
MyDocument::OnShowMenu(wxCommandEvent& event)
{
	int uid = event.GetId();
	if (mol == NULL || mol->mview == NULL)
		return;
	switch (uid) {
		case myMenuID_ShowUnitCell:
			mol->mview->showUnitCell = !mol->mview->showUnitCell;
			break;
		case myMenuID_ShowPeriodicBox:
			mol->mview->showPeriodicBox = !mol->mview->showPeriodicBox;
			break;
		case myMenuID_ShowHydrogens:
			mol->mview->showHydrogens = !mol->mview->showHydrogens;
			break;
		case myMenuID_ShowDummyAtoms:
			mol->mview->showDummyAtoms = !mol->mview->showDummyAtoms;
			break;
		case myMenuID_ShowExpandedAtoms:
			mol->mview->showExpandedAtoms = !mol->mview->showExpandedAtoms;
			break;
		case myMenuID_ShowEllipsoids:
			mol->mview->showEllipsoids = !mol->mview->showEllipsoids;
			break;
		case myMenuID_ShowRotationCenter:
			mol->mview->showRotationCenter = !mol->mview->showRotationCenter;
			break;
	}
	MainViewCallback_setNeedsDisplay(mol->mview, 1);
}
*/

void
MyDocument::OnShowAllAtoms(wxCommandEvent &event)
{
	if (mol == NULL || mol->mview == NULL)
		return;
	MoleculeShowAllAtoms(mol);
}

void
MyDocument::OnHideSelected(wxCommandEvent &event)
{
	IntGroup *ig;
	if (mol == NULL || mol->mview == NULL)
		return;
	ig = MoleculeGetSelection(mol);
	MoleculeHideAtoms(mol, ig);	
}

void
MyDocument::OnHideUnselected(wxCommandEvent &event)
{
	IntGroup *ig;
	if (mol == NULL || mol->mview == NULL)
		return;
	ig = MoleculeGetSelection(mol);
	ig = IntGroupNewFromIntGroup(ig);
	IntGroupReverse(ig, 0, mol->natoms);
	MoleculeHideAtoms(mol, ig);	
	IntGroupRelease(ig);
}

void
MyDocument::OnHideReverse(wxCommandEvent &event)
{
	if (mol == NULL || mol->mview == NULL)
		return;
	MoleculeShowReverse(mol);
}

/*
void
MyDocument::OnShowGraphite(wxCommandEvent &event)
{
	MoleculeLock(mol);
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_show_graphite");
	MoleculeUnlock(mol);
}

void
MyDocument::OnToggleLineMode(wxCommandEvent &event)
{
	mol->mview->lineMode = !mol->mview->lineMode;
	MainViewCallback_setNeedsDisplay(mol->mview, 1);
}
*/

/*  Check whether subthread is running  */
static int
sCheckIsSubThreadRunning(Molecule *mol, int n)
{
	if (mol->mutex != NULL) {
		const char *mes;
		switch (n) {
			case 1: mes = "MM/MD is already running."; break;
			case 2: mes = "Quantum chemistry calculation is already running."; break;
			default: mes = "Some subprocess is already running."; break;
		}
		MyAppCallback_errorMessageBox(mes);
		return 1;
	}
	return 0;
}

/*   Run MD within a subthread  */
static int
sDoMolecularDynamics(void *argptr, int argnum)
{
	MyDocument *doc = (MyDocument *)argptr;
	Molecule *mol = doc->GetMolecule();
	int count, minimize, i, r;
	if (argnum >= 0) {
		count = argnum;
		minimize = 0;
	} else {
		count = -argnum;
		minimize = 1;
	}
	if (count == 0) {
		mol->arena->end_step = mol->arena->start_step;
		md_main(mol->arena, minimize);
	} else if (count > 0) {
		wxCommandEvent insertFrameEvent(MyDocumentEvent, MyDocumentEvent_insertFrameFromMD);
		for (i = 0; i < count; i++) {
			
			mol->arena->end_step = mol->arena->start_step + mol->arena->coord_output_freq;
			r = md_main(mol->arena, minimize);

			if (r == 0) {
				if (mol->requestAbortThread)
					r = -1;
				else {
					/*  Copy the coordinate to the ring buffer  */
					MDRing *ring = mol->arena->ring;
					Vector *rp = ring->buf + ring->size * ring->next;
					Int j;
					Atom *ap;
					MoleculeLock(mol);
					for (j = 0, ap = mol->arena->mol->atoms; j < mol->natoms; j++, ap = ATOM_NEXT(ap)) {
						rp[j] = ap->r;
					}
					if (j < ring->size) {
						XtalCell *cp = mol->arena->mol->cell;
						if (cp != NULL) {
							rp[j++] = cp->axes[0];
							rp[j++] = cp->axes[1];
							rp[j++] = cp->axes[2];
							rp[j++] = cp->origin;
						}
					}
					ring->next = (ring->next + 1) % ring->nframes;
					if (ring->count < ring->nframes)
						ring->count++;
					MoleculeUnlock(mol);
					
					if (minimize && mol->arena->minimize_complete) {
						r = -2;  /*  Minimization complete  */
						break;
					}
					wxPostEvent(doc, insertFrameEvent);
				}
			}
			if (r != 0)
				break;
			if (wxThread::This()->TestDestroy())
				return 0; /* Abnormal termination */
		}
	}
	wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_threadTerminated);
	wxPostEvent(doc, myEvent);
	return 0;
}

void
MyDocument::DoMDOrMinimize(int minimize)
{
	Int n;
	char buf[4096];
	if (sCheckIsSubThreadRunning(mol, subThreadKind))
		return;

	/*  Update the path information of the molecule before MD setup  */
	MoleculeCallback_pathName(mol, buf, sizeof buf);
	MoleculeSetPath(mol, buf);
	
	MolActionCreateAndPerform(mol, SCRIPT_ACTION("b;i"), "cmd_md", minimize, &n);
	if (n < 0)
		return;  /*  Canceled  */
	
	/*  Check whether any bond/angle/torsion are very distant from the equilibrium values  */
	{
	}
	
	mol->mutex = new wxMutex;
	subThreadKind = 1;
	BeginUndoGrouping();
	mol->requestAbortThread = 0;
	MoleculeCallback_disableModificationFromGUI(mol);
	
	if (mol->mview != NULL && mol->mview->ref != NULL) {
		((MoleculeView *)(mol->mview->ref))->EnableProgressIndicator(true);
	}
	wxGetApp().EnableTimerForDocument(this);

	MyThread::DetachNewThread(sDoMolecularDynamics, NULL, (void *)this, (minimize ? -n : n));
}

void
MyDocument::OnMolecularDynamics(wxCommandEvent &event)
{
	DoMDOrMinimize(0);
}

void
MyDocument::OnMinimize(wxCommandEvent &event)
{
	DoMDOrMinimize(1);
}

void
MyDocument::OnStopMDRun(wxCommandEvent &event)
{
	if (mol != NULL && mol->mutex != NULL)
		mol->requestAbortThread = 1;
}

void
MyDocument::OnInsertFrameFromMD(wxCommandEvent &event)
{
	Int i, j, n, old_nframes;
	Atom *ap;
	MDRing *ring;

	/*  Create new frame(s) and copy the new coordinates from the ring buffer  */
	MoleculeLock(mol);
	ring = mol->arena->ring;
	n = ring->count;
	if (n > 0) {
		IntGroup *ig;
		Vector *rp;
		old_nframes = MoleculeGetNumberOfFrames(mol);
		/*  It is more convenient to set cell parameter when inserting frames, whereas 
		    the coordinates can be set afterwards  */
		if (ring->size > mol->natoms) {
			rp = (Vector *)calloc(sizeof(Vector) * 4, n);
			for (i = 0; i < n; i++) {
				j = ((ring->next - n + i + ring->nframes) % ring->nframes) * ring->size + mol->natoms;
				rp[i * 4] = ring->buf[j++];
				rp[i * 4 + 1] = ring->buf[j++];
				rp[i * 4 + 2] = ring->buf[j++];
				rp[i * 4 + 3] = ring->buf[j++];
			}
		} else rp = NULL;
		ig = IntGroupNewWithPoints(old_nframes, n, -1);
		MolActionCreateAndPerform(mol, gMolActionInsertFrames, ig, 0, NULL, (rp != NULL ? n * 4 : 0), rp);
		if (rp != NULL)
			free(rp);
		IntGroupRelease(ig);
		for (i = 0; i < n; i++) {
			MoleculeSelectFrame(mol, old_nframes + i, 1);
			rp = ring->buf + ((ring->next - n + i + ring->nframes) % ring->nframes) * ring->size;
			for (j = 0, ap = mol->atoms; j < mol->natoms; j++, ap = ATOM_NEXT(ap))
				ap->r = rp[j];
		}
		ring->count = 0;
		mol->needsMDCopyCoordinates = 0;  /*  This flag needs to be negated because the coordinates come from the MD run  */
        MainViewCallback_updateCanvas(GetMainView());
	}
	MoleculeUnlock(mol);
}

void
MyDocument::OnUpdateDisplay(wxCommandEvent &event)
{
	MainView *mview = GetMainView();
	MainViewCallback_setNeedsDisplay(mview, 1);
}

void
MyDocument::OnSubThreadTerminated(wxCommandEvent &event)
{
	if (mol != NULL && mol->mutex != NULL) {
		delete (wxMutex *)mol->mutex;
		mol->mutex = NULL;
		mol->requestAbortThread = 0;
		EndUndoGrouping();
		subThreadKind = 0;

		if (mol->mview != NULL && mol->mview->ref != NULL) {
			((MoleculeView *)(mol->mview->ref))->EnableProgressIndicator(false);
		}
		
		wxGetApp().DisableTimerForDocument(this);

		if (mol->arena != NULL && mol->arena->errmsg[0] != 0)
			MyAppCallback_errorMessageBox("MD Error: %s", mol->arena->errmsg);
		
		MoleculeCallback_enableModificationFromGUI(mol);
		if (mol->mview != NULL && mol->mview->ref != NULL) {
			((MoleculeView *)(mol->mview->ref))->InvalidateProgressIndicator();
		}
		
	}
}

/*  Run a subprocess asynchronically  */
long
MyDocument::RunSubProcess(const char *cmd, int (*callback)(Molecule *, int), int (*timerCallback)(Molecule *, int), FILE *output, FILE *errout)
{
	if (sCheckIsSubThreadRunning(mol, subThreadKind))
		return -1;  /*  subProcess (or MM/MD subThread) is already running  */

	mol->mutex = new wxMutex;
	mol->requestAbortThread = 0;
	
	wxString cmdstr(cmd, WX_DEFAULT_CONV);
	subProcess = new wxProcess(this, -1);
	subProcess->Redirect();
	subProcessStdout = output;
	subProcessStderr = errout;

	subProcessPID = ::wxExecute(cmdstr, wxEXEC_ASYNC | wxEXEC_MAKE_GROUP_LEADER, subProcess);
	if (subProcessPID == 0) {
		subProcess->Detach();
		subProcess = NULL;
		delete (wxMutex *)(mol->mutex);
		mol->mutex = NULL;
		subThreadKind = 0;
		return -2;  /*  Cannot start subProcess  */
	}

	if (mol->mview != NULL && mol->mview->ref != NULL) {
		((MoleculeView *)(mol->mview->ref))->EnableProgressIndicator(true);
	}
	subThreadKind = 2;
	mol->requestAbortThread = 0;
	MoleculeCallback_disableModificationFromGUI(mol);
	BeginUndoGrouping();
	wxGetApp().EnableTimerForDocument(this);
	endSubProcessCallback = callback;
	timerSubProcessCallback = timerCallback;

	return subProcessPID;
}

void
MyDocument::FlushSubProcessOutput()
{
	wxInputStream *stream;
	char buf[1024];
	int len;
	if (subProcess == NULL)
		return;  /*  Do nothing  */
	stream = subProcess->GetInputStream();
	if (subProcessStdout != NULL && stream != NULL && stream->CanRead()) {
		stream->Read(buf, sizeof buf - 1);
		len = stream->LastRead();
		if (len > 0) {
			buf[len] = 0;
			if (subProcessStdout == (FILE *)1) {
				MyAppCallback_setConsoleColor(0);
				MyAppCallback_showScriptMessage("%s", buf);
			} else {
				fwrite(buf, 1, len, subProcessStdout);
			}
		}
	}
	stream = subProcess->GetErrorStream();
	if (subProcessStderr != NULL && stream != NULL && stream->CanRead()) {
		stream->Read(buf, sizeof buf - 1);
		len = stream->LastRead();
		if (len > 0) {
			buf[len] = 0;
			if (subProcessStderr == (FILE *)1) {
				MyAppCallback_setConsoleColor(1);
				MyAppCallback_showScriptMessage("%s", buf);
				MyAppCallback_setConsoleColor(0);
			} else {
				fwrite(buf, 1, len, subProcessStderr);
			}
		}
	}	
}

void
MyDocument::OnEndSubProcess(wxProcessEvent &event)
{
	if (mol != NULL && mol->mutex != NULL) {
		
		FlushSubProcessOutput();
		if (subProcessStdout != NULL && subProcessStdout != (FILE *)1)
			fclose(subProcessStdout);
		if (subProcessStderr != NULL && subProcessStderr != (FILE *)1)
			fclose(subProcessStderr);
		subProcessStdout = subProcessStderr = NULL;
	
		delete (wxMutex *)mol->mutex;
		mol->mutex = NULL;
		mol->requestAbortThread = 0;
		EndUndoGrouping();
		subThreadKind = 0;

		delete subProcess;
		subProcess = NULL;
		
		wxGetApp().DisableTimerForDocument(this);
		
		MoleculeCallback_enableModificationFromGUI(mol);
		if (mol->mview != NULL && mol->mview->ref != NULL) {
			((MoleculeView *)(mol->mview->ref))->EnableProgressIndicator(false);
		}
		if (endSubProcessCallback != NULL) {
			(*endSubProcessCallback)(mol, event.GetExitCode());
			endSubProcessCallback = NULL;
		}
		timerSubProcessCallback = NULL;
	}	
}

static wxString
sCreateTemporaryLogDirectoryForAC(const wxString& filename)
{
	char *log_dir;
	int i, status;

	/*  Extract the name  */
	wxFileName fname(filename);
	wxString name = fname.GetName();

	status = MyAppCallback_getGlobalSettingsWithType("antechamber.log_dir", 's', &log_dir);
	if (status) {
		char *hdir = MyAppCallback_getDocumentHomeDir();
		asprintf(&log_dir, "%s/antechamber", (hdir ? hdir : ""));
		if (hdir != NULL)
			free(hdir);
	}
	fix_dosish_path(log_dir);
	
	wxString tdir;

	/*  Prepare the log directory  */
	wxString dirname(log_dir, wxConvFile);
	if (!wxFileName::Mkdir(dirname, 0777, wxPATH_MKDIR_FULL)) {
		MyAppCallback_errorMessageBox("Cannot create log directory '%s'", log_dir);
		free(log_dir);
		return tdir;  /*  empty  */
	}
	free(log_dir);

	for (i = 0; i < 1000; i++) {
		tdir = dirname + wxFileName::GetPathSeparator() + name + wxString::Format(_T("_%04d"), i);
		if (!wxFileName::DirExists(tdir))
			break;
	}
	if (i >= 1000 || !wxFileName::Mkdir(tdir)) {
		MyAppCallback_errorMessageBox("Cannot create temporary files. Please make sure the log directory has enough space for writing.");
		tdir.Empty();
	}
	return tdir;
}

static bool
sRemoveDirectoryRecursively(const wxString &dir)
{
	wxString name, file;
	wxArrayString files;
	int i, n;
	n = 0;
	{
		/*  The GetFirst/GetNext loop should not be mixed with ::wxRemoveFile or ::wxRmdir  */
		wxDir wdir(dir);
		if (wdir.GetFirst(&name)) {
			do {
				file = dir + wxFileName::GetPathSeparator() + name;
				files.Add(file);
				n++;
			} while (wdir.GetNext(&name));
		}
	}
	for (i = 0; i < n; i++) {
		file = files[i];
		if (wxDir::Exists(file)) {
			if (!sRemoveDirectoryRecursively(file))
				return false;
		} else {
			if (!::wxRemoveFile(file))
				return false;
		}
	}
	return ::wxRmdir(dir);
}

static int
sEraseLogFiles(const wxString& tdir, int status)
{
	bool success = true;
	Int log_keep_number, n, i, j;
	char *log_level;
	wxString dir2;

	if (MyAppCallback_getGlobalSettingsWithType("antechamber.log_level", 's', &log_level) != 0)
		log_level = NULL;
	if (MyAppCallback_getGlobalSettingsWithType("antechamber.log_keep_number", 'i', &log_keep_number) != 0)
		log_keep_number = 5;
	if (log_level == NULL || strcmp(log_level, "none") == 0 || (strcmp(log_level, "error_only") == 0 && status == 0)) {
		//  Erase the present log
		if (!sRemoveDirectoryRecursively(tdir)) {
			success = false;
			dir2 = tdir;
		}
	} else if (strcmp(log_level, "latest") == 0) {
		wxString dirname = tdir.BeforeLast(wxFileName::GetPathSeparator());
		wxDir wdir(dirname);
		wxString name;
		wxArrayString files;
		n = 0;
		if (wdir.GetFirst(&name)) {
			do {
				wxString fullname = dirname + wxFileName::GetPathSeparator() + name;
				if (wxDir::Exists(fullname)) {
					files.Add(fullname);
					n++;
				}
			} while (wdir.GetNext(&name));
		}
		if (n > log_keep_number) {
			//  Sort directories by creation date
			struct temp_struct { time_t tm; int idx; } *tp;
			tp = (struct temp_struct *)malloc(sizeof(struct temp_struct) * n);
			for (i = 0; i < n; i++) {
				wxFileName fn(files[i], wxEmptyString);
				wxDateTime dt;
				j = fn.GetTimes(NULL, NULL, &dt);
				tp[i].tm = dt.GetTicks();
				tp[i].idx = i;
			}
			for (i = 0; i < n; i++) {
				struct temp_struct temp;
				int k = i;
				for (j = i + 1; j < n; j++) {
					if (tp[j].tm < tp[k].tm)
						k = j;
				}
				if (k != i) {
					temp = tp[k];
					tp[k] = tp[i];
					tp[i] = temp;
				}
			}
			//  Keep last log_keep_number and delete the rest
			for (i = 0; i < n - log_keep_number; i++) {
				if (!sRemoveDirectoryRecursively(files[tp[i].idx])) {
					success = false;
					dir2 = files[tp[i].idx];
					break;
				}
			}
		}
	}
	
	if (success) {
		return 0;
	} else {
		MyAppCallback_errorMessageBox("Error during deleting log file '%s'", (const char *)dir2.mb_str(wxConvFile));
		return -1;
	}
}

void
MyDocument::OnUpdateUI(wxUpdateUIEvent& event)
{
	int uid = event.GetId();
	IntGroup *ig = MoleculeGetSelection(mol);
	Int nselected = (ig == NULL ? 0 : IntGroupGetCount(ig));
	switch (uid) {
		case wxID_COPY:
		case wxID_CUT:
		case wxID_PASTE:
		case wxID_DELETE:
			event.Enable(true);
			return;			
		case myMenuID_Import:
			event.Enable(true);
			return;
		case wxID_SELECTALL:
			event.Enable(true);
			return;
		case myMenuID_SelectFragment:
			event.Enable(nselected > 0);
			return;
		case myMenuID_SelectReverse:
			event.Enable(true);
			return;
		case myMenuID_CreatePiAnchor:
			event.Enable(nselected > 0);
			return;
		case myMenuID_AddHydrogenSp3:
		case myMenuID_AddHydrogenSp2:
		case myMenuID_AddHydrogenLinear:
		case myMenuID_AddHydrogenPyramidal:
		case myMenuID_AddHydrogenBent:
			event.Enable(nselected > 0);
			return;
		case myMenuID_FitToScreen:
			event.Enable(true);
			return;
	/*	case myMenuID_ShowUnitCell:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showUnitCell != 0);
			return;
		case myMenuID_ShowPeriodicBox:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showPeriodicBox != 0);
			return;
		case myMenuID_ShowHydrogens:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showHydrogens != 0);
			return;
		case myMenuID_ShowDummyAtoms:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showDummyAtoms != 0);
			return;
		case myMenuID_ShowExpandedAtoms:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showExpandedAtoms != 0);
			return;
		case myMenuID_ShowEllipsoids:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showEllipsoids != 0);
			return;
		case myMenuID_ShowRotationCenter:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->showRotationCenter != 0);
			return;
	 */
		case myMenuID_ShowAllAtoms:
		case myMenuID_HideReverse:
			event.Enable(mol->mview != NULL && mol->mview->countHidden > 0);
			return;
		case myMenuID_HideSelected:
		case myMenuID_HideUnselected:
			event.Enable(nselected > 0);
			return;
	/*	case myMenuID_LineMode:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->lineMode != 0);
			return;
	 */
		case myMenuID_MolecularDynamics:
		case myMenuID_Minimize:
			if (mol != NULL && mol->mutex == NULL)
				event.Enable(true);
			else event.Enable(false);
			return;
		case myMenuID_StopMDRun:
			if (mol != NULL && mol->mutex != NULL)
				event.Enable(true);
			else event.Enable(false);
			return;
	}
	event.Skip();
}

void
MyDocument::TimerCallback(int timerCount)
{
	if (mol != NULL && mol->mview != NULL && mol->mview->ref != NULL) {
		((MoleculeView *)(mol->mview->ref))->ProceedProgressIndicator();
		if (subProcess != NULL) {
			FlushSubProcessOutput();
			if (timerSubProcessCallback != NULL) {
				if ((*timerSubProcessCallback)(mol, timerCount) != 0)
					mol->requestAbortThread = 1;
			}
			if (mol->requestAbortThread) {
				/*  Try to terminate the subprocess gently  */
				wxProcess::Kill(subProcessPID, wxSIGTERM, wxKILL_CHILDREN);
			}
		}
	}
}

#pragma mark ====== Plain C Interface ======

MyDocument *
MyDocumentFromMolecule(Molecule *mp)
{
  void *ref;
  if (mp != NULL && mp->mview != NULL && (ref = mp->mview->ref) != NULL)
    return ((MoleculeView *)ref)->MolDocument();
  else return NULL;
}

Molecule *
MoleculeCallback_openNewMolecule(const char *fname)
{
  if (!gUseGUI)
    return NULL;
	wxDocument *doc;
	MyDocManager *manager = wxGetApp().DocManager();
	if (fname == NULL || *fname == 0) {
		doc = manager->CreateDocument(wxT(""), wxDOC_NEW);
	} else {
		wxString fnamestr(fname, wxConvFile);
		doc = manager->CreateDocument(fnamestr, wxDOC_SILENT);
	}
	if (doc == NULL)
		return NULL;
	else return ((MyDocument *)doc)->GetMolecule();
}

void
MoleculeCallback_notifyModification(Molecule *mp, int now_flag)
{
  if (!gUseGUI)
    return;
	MyDocument *doc = MyDocumentFromMolecule(mp);
	if (doc && !doc->isModifyNotificationSent) {
		doc->isModifyNotificationSent = true;
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_documentModified);
		myEvent.SetEventObject(doc);
		if (now_flag)
			doc->ProcessEvent(myEvent);
		else
			wxPostEvent(doc, myEvent);
	}
}

static wxDataObject *
sMoleculePasteboardObjectOfType(const char *type, const void *data, int length)
{
	if (strcmp(type, "TEXT") == 0) {
		wxTextDataObject *tp = new wxTextDataObject();
		if (data != NULL) {
			wxString str((const char *)data, WX_DEFAULT_CONV, length);
			tp->SetText(str);
		}
		return tp;
	} else {
		MyClipboardData *dp = new MyClipboardData(type);
		if (data != NULL)
			dp->SetData(length, data);
		return dp;
	}
}

/*  Write to pasteboard. NOTE: data must be a malloc'ed pointer and its ownership
    will be taken by the pasteboard. */
int
MoleculeCallback_writeToPasteboard(const char *type, const void *data, int length)
{
	int retval = 1;
	if (wxTheClipboard->Open()) {
		wxTheClipboard->SetData(sMoleculePasteboardObjectOfType(type, data, length));
	/*	MyClipboardData *myData = new MyClipboardData();
		if (myData->SetData(length, data)) {
			wxTheClipboard->SetData(myData);
			retval = 0;
		} else
			delete myData; */
		wxTheClipboard->Close();
		retval = 0;
	}
	return retval;
}

int
MoleculeCallback_writeToPasteBoardInMultipleFormats(const char *type1, const void *data1, int length1, const char *type2, ...)
{
	return 0;
}

int
MoleculeCallback_readFromPasteboard(const char *type, void **dptr, int *length)
{
	int retval = 1;
	int len;
	void *p;
	if (wxTheClipboard->Open()) {
		wxDataObject *dp = sMoleculePasteboardObjectOfType(type, NULL, 0);
		if (wxTheClipboard->GetData(*dp)) {
			if (strcmp(type, "TEXT") == 0) {
				wxTextDataObject *tp = (wxTextDataObject *)dp;
				wxString str = tp->GetText();
				const char *cp = str.mb_str(WX_DEFAULT_CONV);
				len = strlen(cp);
				p = malloc(len + 1);
				if (p != NULL) {
					strcpy((char *)p, cp);
					*dptr = p;
					*length = len;
					retval = 0;
				}
				delete tp;
			} else {
				MyClipboardData *mp = (MyClipboardData *)dp;
				len = mp->GetDataSize();
				p = malloc(len); 
				if (p != NULL) {
					mp->GetDataHere(p);
					*dptr = p;
					*length = len;
					retval = 0;
				}
				delete mp;
			}
		}
		wxTheClipboard->Close();
	}
	return retval;
}

int
MoleculeCallback_isDataInPasteboard(const char *type)
{
	if (strcmp(type, "TEXT") == 0)
		return wxTheClipboard->IsSupported(wxDF_TEXT);
	else {
		MyClipboardData myData(type);
		return wxTheClipboard->IsSupported(myData.GetFormat());
	}
}

Molecule *
MoleculeCallback_currentMolecule(void)
{
  if (!gUseGUI)
    return NULL;
  MainView *mview = MainViewCallback_activeView();
  if (mview != NULL)
    return mview->mol;
  else return NULL;
}

Molecule *
MoleculeCallback_moleculeAtIndex(int idx)
{
  if (!gUseGUI)
    return NULL;
  MainView *mview = MainViewCallback_viewWithTag(idx);
  if (mview != NULL)
    return mview->mol;
  else return NULL;
}

Molecule *
MoleculeCallback_moleculeAtOrderedIndex(int idx)
{
  return MoleculeCallback_moleculeAtIndex(idx);
}

void
MoleculeCallback_displayName(Molecule *mol, char *buf, int bufsize)
{
  if (!gUseGUI) {
    buf[0] = 0;
    return;
  }
  MyDocument *doc = MyDocumentFromMolecule(mol);
  if (doc != NULL) {
    wxString fname;
    fname = doc->GetUserReadableName();
    strncpy(buf, (const char*)fname.mb_str(wxConvFile), bufsize - 1);
    buf[bufsize - 1] = 0;
  } else {
    buf[0] = 0;
  }
}

void
MoleculeCallback_pathName(Molecule *mol, char *buf, int bufsize)
{
  if (!gUseGUI) {
    if (mol != NULL && mol->path != NULL) {
      strncpy(buf, mol->path, bufsize - 1);
      buf[bufsize - 1] = 0;
    } else buf[0] = 0;
    return;
  }
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL && doc->hasFile)
		MainViewCallback_getFilename(mol->mview, buf, bufsize);
	else buf[0] = 0;
}

int
MoleculeCallback_setDisplayName(Molecule *mol, const char *name)
{
  if (!gUseGUI) {
    return 0;
  }
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc == NULL || doc->hasFile)
		return 1; /*  Cannot change file-associated window title  */
	wxString fname(name, wxConvFile);
	doc->SetTitle(fname);
	doc->GetFirstView()->OnChangeFilename();
	return 0;
}

void
MoleculeCallback_lockMutex(void *mutex)
{
  if (gUseGUI)
    ((wxMutex *)mutex)->Lock();
}

void
MoleculeCallback_unlockMutex(void *mutex)
{
  if (gUseGUI)
    ((wxMutex *)mutex)->Unlock();
}

void
MoleculeCallback_disableModificationFromGUI(Molecule *mol)
{
	mol->dontModifyFromGUI = 1;
	if (mol->mview != NULL) {
		if (mol->mview->mode == kTrackballCreateMode || mol->mview->mode == kTrackballEraseMode) {
			MainView_setMode(mol->mview, kTrackballSelectionMode);
			MainViewCallback_selectMatrixCellForMode(mol->mview, kTrackballSelectionMode);
		}
		MainViewCallback_enableToggleButton(mol->mview, kTrackballCreateMode, false);
		MainViewCallback_enableToggleButton(mol->mview, kTrackballEraseMode, false);
	}
}

void
MoleculeCallback_enableModificationFromGUI(Molecule *mol)
{
	mol->dontModifyFromGUI = 0;
	if (mol->mview != NULL) {
		MainViewCallback_enableToggleButton(mol->mview, kTrackballCreateMode, true);
		MainViewCallback_enableToggleButton(mol->mview, kTrackballEraseMode, true);
	}
}

void
MoleculeCallback_cannotModifyMoleculeDuringMDError(Molecule *mol)
{
	MyAppCallback_errorMessageBox("Cannot modify molecule during MD");
}

int
MoleculeCallback_callSubProcessAsync(Molecule *mol, const char *cmd, int (*callback)(Molecule *, int), int (*timerCallback)(Molecule *, int), FILE *output, FILE *errout)
{
  if (!gUseGUI)
    return -1;
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL)
		return doc->RunSubProcess(cmd, callback, timerCallback, output, errout);
	else return -1;
}

void
MolActionCallback_registerUndo(Molecule *mol, MolAction *action)
{
  if (!gUseGUI)
    return;
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL && doc->IsUndoEnabled())
		doc->PushUndoAction(action);
}

int
MolActionCallback_setUndoRegistrationEnabled(Molecule *mol, int flag)
{
  if (!gUseGUI)
    return 0;
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL) {
		doc->SetUndoEnabled(flag);
		return (doc->IsUndoEnabled() ? 1 : 0);
	} else return 0;
}

int
MolActionCallback_isUndoRegistrationEnabled(Molecule *mol)
{
  if (!gUseGUI)
    return 0;
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL && doc->IsUndoEnabled())
		return 1;
	else return 0;
}

