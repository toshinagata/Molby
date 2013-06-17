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
	EVT_MENU(myMenuID_Import, MyDocument::OnImport)
	EVT_MENU(myMenuID_Export, MyDocument::OnExport)
	EVT_MENU(wxID_COPY, MyDocument::OnCopy)
	EVT_MENU(wxID_PASTE, MyDocument::OnPaste)
	EVT_MENU(wxID_CUT, MyDocument::OnCut)
	EVT_MENU(wxID_DELETE, MyDocument::OnDelete)
	EVT_MENU(myMenuID_CreateNewAtom, MyDocument::OnCreateNewAtom)
	EVT_MENU_RANGE(myMenuID_CreateNewVdwParameter, myMenuID_CreateNewVdwCutoffParameter, MyDocument::OnCreateNewParameter)
	EVT_MENU(myMenuID_CreatePiAnchor, MyDocument::OnCreatePiAnchor)
	EVT_MENU(wxID_SELECTALL, MyDocument::OnSelectAll)
	EVT_MENU(myMenuID_SelectFragment, MyDocument::OnSelectFragment)
	EVT_MENU(myMenuID_SelectReverse, MyDocument::OnSelectReverse)
	EVT_MENU(myMenuID_FitToScreen, MyDocument::OnFitToScreen)
	EVT_MENU(myMenuID_CenterSelection, MyDocument::OnCenterSelection)
	EVT_MENU(myMenuID_ShowUnitCell, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowPeriodicBox, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowHydrogens, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowDummyAtoms, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowExpandedAtoms, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowEllipsoids, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowRotationCenter, MyDocument::OnShowMenu)
	EVT_MENU(myMenuID_ShowGraphite, MyDocument::OnShowGraphite)
	EVT_MENU(myMenuID_LineMode, MyDocument::OnToggleLineMode)
	EVT_MENU_RANGE(myMenuID_AddHydrogenSp3, myMenuID_AddHydrogenBent, MyDocument::OnAddHydrogen)
	EVT_UPDATE_UI_RANGE(myMenuID_MyFirstMenuItem, myMenuID_MyLastMenuItem, MyDocument::OnUpdateUI)
	EVT_MENU(myMenuID_MolecularDynamics, MyDocument::OnMolecularDynamics)
	EVT_MENU(myMenuID_Minimize, MyDocument::OnMinimize)
	EVT_MENU(myMenuID_StopMDRun, MyDocument::OnStopMDRun)
	EVT_MENU(myMenuID_DefinePeriodicBox, MyDocument::OnDefinePeriodicBox)
	EVT_MENU(myMenuID_ShowPeriodicImage, MyDocument::OnShowPeriodicImage)
	EVT_MENU(myMenuID_PressureControl, MyDocument::OnPressureControl)
	EVT_MENU(myMenuID_DefineSymmetry, MyDocument::OnDefineSymmetry)
	EVT_MENU(myMenuID_ExpandBySymmetry, MyDocument::OnExpandBySymmetry)
	EVT_MENU(myMenuID_RunAntechamber, MyDocument::OnInvokeAntechamber)
	EVT_MENU(myMenuID_RunResp, MyDocument::OnInvokeResp)
	EVT_MENU(myMenuID_CreateSanderInput, MyDocument::OnCreateSanderInput)
	EVT_MENU(myMenuID_ImportAmberFrcmod, MyDocument::OnImportAmberFrcmod)
	EVT_MENU(myMenuID_CreateGamessInput, MyDocument::OnCreateGamessInput)
	EVT_MENU(myMenuID_CreateMOCube, MyDocument::OnCreateMOCube)
	EVT_MENU(myMenuID_ShowAllAtoms, MyDocument::OnShowAllAtoms)
	EVT_MENU(myMenuID_HideReverse, MyDocument::OnHideReverse)
	EVT_MENU(myMenuID_HideSelected, MyDocument::OnHideSelected)
	EVT_MENU(myMenuID_HideUnselected, MyDocument::OnHideUnselected)
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
}

MyDocument::~MyDocument()
{
	int i;
	Molecule *mol2 = mol;
	mol = NULL;

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
}

/*
MainView *
MyDocument::GetMainView()
{
	MoleculeView *view = (MoleculeView *)GetFirstView();
	if (view != NULL)
		return view->mview;
	else return NULL;
}
*/

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
	wxString fnstr;
	GetPrintableName(fnstr);
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
//	printf("isDirty = %d\n", (GetCommandProcessor()->IsDirty()));
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
				cmd = new MyCommand(mol);
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
	return wxDocument::Close();
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
//	printf("MyDocument::OnDocumentModified invoked\n");
	isModifyNotificationSent = false;
	MoleculeClearModifyCount(GetMainView()->mol);
	
	event.Skip();  //  Also pass to other notification handlers
	UpdateModifyFlag();
}

void
MyDocument::OnCopy(wxCommandEvent& event)
{
	wxWindow *focusWindow = wxWindow::FindFocus();
/*	printf("focus window class = %ls\n", focusWindow->GetClassInfo()->GetClassName());  */
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
/*	MainViewCallback_ensureVisible(mview, row); */ /* Invoked from startEditText */
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
			"proc { |g| create_pi_anchor('AN', g).index rescue -1 }",
			ig, &idx) != 0)
		return;
	MainViewCallback_selectTable(mview, kMainViewAtomTableIndex);
	ig2 = IntGroupNewFromIntGroup(ig);
	IntGroupAdd(ig2, idx, 1);
	MoleculeSetSelection(mol, ig2);
	IntGroupRelease(ig2);
	MainView_refreshTable(mview);
	row = MainView_indexToTableRow(mview, idx);
	MainViewCallback_startEditText(mview, row, 1);
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

/*  Check whether subthread is running  */
static int
sCheckIsSubThreadRunning(Molecule *mol, int n)
{
	if (mol->mutex != NULL) {
		const char *mes;
		switch (n) {
			case 1: mes = "MM/MD is already running."; break;
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
					
					/*  Create a new frame and copy the new coordinates  */
				/*	MoleculeLock(mol);
					ig = IntGroupNewWithPoints(MoleculeGetNumberOfFrames(mol), 1, -1);
					MolActionCreateAndPerform(mol, gMolActionInsertFrames, ig, 0, NULL);
					IntGroupRelease(ig);					
					md_copy_coordinates_from_internal(mol->arena);
					MoleculeUnlock(mol); */

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
		
		if (mol->arena != NULL && mol->arena->errmsg[0] != 0)
			MyAppCallback_errorMessageBox("MD Error: %s", mol->arena->errmsg);
	}
}

void
MyDocument::OnDefinePeriodicBox(wxCommandEvent &event)
{
	MoleculeLock(mol);
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_define_unit_cell");
	MoleculeUnlock(mol);
}

void
MyDocument::OnShowPeriodicImage(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_show_periodic_image");
}

void
MyDocument::OnPressureControl(wxCommandEvent &event)
{
	MoleculeLock(mol);
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_pressure_control");
	MoleculeUnlock(mol);
}

void
MyDocument::OnDefineSymmetry(wxCommandEvent &event)
{
}

void
MyDocument::OnExpandBySymmetry(wxCommandEvent &event)
{
}

static wxString
sCreateTemporaryLogDirectoryForAC(const wxString& filename)
{
//	char *ante_dir;
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
MyDocument::OnInvokeAntechamber(wxCommandEvent &event)
{
	char *ante_dir, buf[256];
	Int net_charge, i, n, calc_charge, use_residue;
	int status;

	/*  Find the ambertool directory  */
	wxString ante = MyApp::FindResourcePath() + wxFileName::GetPathSeparator() + _T("amber11") + wxFileName::GetPathSeparator() + _T("bin");
	ante_dir = strdup(ante.mb_str(wxConvFile));
	fix_dosish_path(ante_dir);

	/*  Ask for antechamber options and log directory  */
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(";i"), "cmd_antechamber", &n);
	if (!n)
		return;

	if ((status = MyAppCallback_getGlobalSettingsWithType("antechamber.nc", 'i', &net_charge))
		|| (status = MyAppCallback_getGlobalSettingsWithType("antechamber.calc_charge", 'i', &calc_charge))
		|| (status = MyAppCallback_getGlobalSettingsWithType("antechamber.use_residue", 'i', &use_residue))) {
		Molby_showError(status);
		return;
	}

	/*  Prepare the log directory  */
	wxString tdir = sCreateTemporaryLogDirectoryForAC(GetFilename());
	if (tdir.IsEmpty())
		return;

	/*  Move to the temporary directory and export the molecule as a pdb  */
	wxString cwd = wxFileName::GetCwd();
	if (!wxFileName::SetCwd(tdir)) {
		MyAppCallback_errorMessageBox("Cannot move to the temporary directory '%s'", (const char *)tdir.mb_str(WX_DEFAULT_CONV));
		return;
	}
	{
		Int *resno;
		char *resnames;
		Atom *ap;
		char *errbuf;
		resno = (Int *)calloc(sizeof(Int), mol->natoms);
		resnames = (char *)calloc(sizeof(char), mol->natoms * 4);
		if (resno == NULL || resnames == NULL) {
			MyAppCallback_errorMessageBox("Cannot save current residue informations (out of memory)");
			return;
		}
		if (!use_residue) {
			for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
				resno[i] = ap->resSeq;
				memmove(resnames + i * 4, ap->resName, 4);
				ap->resSeq = 1;
				memmove(resnames + i * 4, "RES", 4);
			}
		}
		n = MoleculeWriteToPdbFile(mol, "mol.pdb", &errbuf);
		if (!use_residue) {
			for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
				ap->resSeq = resno[i];
				memmove(ap->resName, resnames + i * 4, 4);
			}
		}
		free(resno);
		free(resnames);
		if (n != 0) {
			MyAppCallback_errorMessageBox("PDB export error: %s", errbuf);
			free(errbuf);
			wxFileName::SetCwd(cwd);
			return;
		}
	}
	
	{
		/*  Run antechamber and parmck  */
		char *p;

		/*  Set AMBERHOME environment variable if necessary  */
		n = strlen(ante_dir);
		if (n >= 4) {
			p = ante_dir + n - 3;
			if ((p[-1] == '\\' || p[-1] == '/') && (strcmp(p, "bin") == 0 || strcmp(p, "exe") == 0)) {
				n -= 4;
			}
		}
		snprintf(buf, sizeof buf, "%.*s", n, ante_dir);
		p = getenv("AMBERHOME");
		if (p == NULL || strcmp(p, buf) != 0) {
			asprintf(&p, "AMBERHOME=%s", buf);
			putenv(p);
		}
		
		if (calc_charge) {
			snprintf(buf, sizeof buf, "-nc %d -c bcc", net_charge);
		} else buf[0] = 0;

		asprintf(&p, "\"%s/antechamber\" -i mol.pdb -fi pdb -o mol.ac -fo ac %s", ante_dir, buf);

		status = MyAppCallback_callSubProcess(p, "antechamber", NULL, NULL);
		if (status != 0) {
			MyAppCallback_errorMessageBox("Antechamber failed: status = %d.", status);
		} else {
			asprintf(&p, "\"%s/parmchk\" -i mol.ac -f ac -o frcmod", ante_dir);
			status = MyAppCallback_callSubProcess(p, "parmchk", NULL, NULL);
			if (status != 0)
				MyAppCallback_errorMessageBox("Parmchk failed: status = %d.", status);
		}
	}

	if (status == 0) {
		wxString acfile = tdir + wxFileName::GetPathSeparator() + _T("mol.ac");
		status = MolActionCreateAndPerform(mol, SCRIPT_ACTION("s"), "import_ac", (const char *)acfile.mb_str(wxConvFile));
		if (status != 0) {
			MyAppCallback_errorMessageBox("Cannot import antechamber output.");
		}
	}
	
	if (calc_charge && status == 0) {
		wxString sqmfile = tdir + wxFileName::GetPathSeparator() + _T("sqm.out");
		if (wxFileName::FileExists(sqmfile)) {
			status = MolActionCreateAndPerform(mol, SCRIPT_ACTION("s"), "import_sqmout", (const char *)sqmfile.mb_str(wxConvFile));
			if (status != 0) {
				MyAppCallback_errorMessageBox("Cannot import sqm output.");
			}
		}
	}

	if (status == 0) {
		wxString frcmodfile = tdir + wxFileName::GetPathSeparator() + _T("frcmod");
		status = MolActionCreateAndPerform(mol, SCRIPT_ACTION("s"), "import_frcmod", (const char *)frcmodfile.mb_str(wxConvFile));
		if (status != 0) {
			MyAppCallback_errorMessageBox("Cannot import parmchk output.");
		}
	}
	
	if (status == 0) {
		/*  Remove improper torsions (they should be rebuilt)  */
		if (mol->nimpropers > 0) {
			IntGroup *ig;
			ig = IntGroupNewWithPoints(0, mol->nimpropers, -1);
			MolActionCreateAndPerform(mol, gMolActionDeleteImpropers, ig);
			IntGroupRelease(ig);
		}
	}
	
	wxFileName::SetCwd(cwd);

	/*  Erase log files  */
	sEraseLogFiles(tdir, status);

	if (status == 0) {
		((MoleculeView *)GetFirstView())->GetListCtrl()->Update();
		MyAppCallback_messageBox("Antechamber succeeded.", "Success", 0, 0);
	}
	MyAppCallback_showRubyPrompt();
}

void
MyDocument::OnInvokeResp(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_gamess_resp");
}

void
MyDocument::OnCreateSanderInput(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "export_prmtop");
}

void
MyDocument::OnImportAmberFrcmod(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_import_frcmod");
}

void
MyDocument::OnCreateGamessInput(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_create_gamess_input");
}

void
MyDocument::OnCreateMOCube(wxCommandEvent &event)
{
	MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), "cmd_create_cube");	
}

void
MyDocument::OnUpdateUI(wxUpdateUIEvent& event)
{
	int uid = event.GetId();
	IntGroup *ig = MoleculeGetSelection(mol);
	Int nselected = (ig == NULL ? 0 : IntGroupGetCount(ig));
//	wxMenuItem *item = (wxMenuItem *)event.GetEventObject();
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
		case myMenuID_ShowUnitCell:
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
		case myMenuID_ShowAllAtoms:
		case myMenuID_HideReverse:
			event.Enable(mol->mview != NULL && mol->mview->countHidden > 0);
			return;
		case myMenuID_HideSelected:
		case myMenuID_HideUnselected:
			event.Enable(nselected > 0);
			return;
		case myMenuID_LineMode:
			event.Enable(true);
			event.Check(mol != NULL && mol->mview != NULL && mol->mview->lineMode != 0);
			return;
		case myMenuID_MolecularDynamics:
		case myMenuID_Minimize:
		case myMenuID_DefinePeriodicBox:
		case myMenuID_PressureControl:
			if (mol != NULL && mol->mutex == NULL)
				event.Enable(true);
			else event.Enable(false);
			return;
		case myMenuID_StopMDRun:
			if (mol != NULL && mol->mutex != NULL)
				event.Enable(true);
			else event.Enable(false);
			return;
		case myMenuID_ShowPeriodicImage:
			event.Enable(true);
			return;
		case myMenuID_RunAntechamber:
		case myMenuID_RunResp:
		case myMenuID_CreateSanderInput:
			if (mol != NULL && mol->natoms > 0)
				event.Enable(true);
			else event.Enable(false);
			return;			
		case myMenuID_CreateGamessInput:
			if (mol != NULL && mol->natoms > 0)
				event.Enable(true);
			else event.Enable(false);
			return;
		case myMenuID_CreateMOCube:
			if (mol == NULL || mol->bset == NULL || mol->bset->natoms == 0)
				event.Enable(false);
			else
				event.Enable(true);
			return;
	}
	event.Skip();
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
	MyDocument *doc = MyDocumentFromMolecule(mp);
	if (doc && !doc->isModifyNotificationSent) {
		doc->isModifyNotificationSent = true;
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_documentModified);
		if (now_flag)
			doc->ProcessEvent(myEvent);
		else
			wxPostEvent(doc, myEvent);
	}
}

/*
static NSString *
sMoleculePasteboardType(const char *type)
{
	static NSMutableArray *array;
	NSString *str;
	unsigned int idx;
	if (array == nil)
		array = [[NSMutableArray array] retain];
	str = [NSString stringWithUTF8String: type];
	idx = [array indexOfObject: str];
	if (idx == NSNotFound) {
		[array addObject: str];
		return str;
	} else {
		return [array objectAtIndex: idx];
	}
}
*/

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
  MainView *mview = MainViewCallback_activeView();
  if (mview != NULL)
    return mview->mol;
  else return NULL;
}

Molecule *
MoleculeCallback_moleculeAtIndex(int idx)
{
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
  MyDocument *doc = MyDocumentFromMolecule(mol);
  if (doc != NULL) {
    wxString fname;
    doc->GetPrintableName(fname);
    strncpy(buf, (const char*)fname.mb_str(wxConvFile), bufsize - 1);
    buf[bufsize - 1] = 0;
  } else {
    buf[0] = 0;
  }
}

void
MoleculeCallback_pathName(Molecule *mol, char *buf, int bufsize)
{
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL && doc->hasFile)
		MainViewCallback_getFilename(mol->mview, buf, bufsize);
	else buf[0] = 0;
}

int
MoleculeCallback_setDisplayName(Molecule *mol, const char *name)
{
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
	((wxMutex *)mutex)->Lock();
}

void
MoleculeCallback_unlockMutex(void *mutex)
{
	((wxMutex *)mutex)->Unlock();
}

void
MoleculeCallback_cannotModifyMoleculeDuringMDError(Molecule *mol)
{
	MyAppCallback_errorMessageBox("Cannot modify molecule during MD");
}

void
MolActionCallback_registerUndo(Molecule *mol, MolAction *action)
{
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL && doc->IsUndoEnabled())
		doc->PushUndoAction(action);
}

int
MolActionCallback_setUndoRegistrationEnabled(Molecule *mol, int flag)
{
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL) {
		doc->SetUndoEnabled(flag);
		return (doc->IsUndoEnabled() ? 1 : 0);
	} else return 0;
}

int
MolActionCallback_isUndoRegistrationEnabled(Molecule *mol)
{
	MyDocument *doc = MyDocumentFromMolecule(mol);
	if (doc != NULL && doc->IsUndoEnabled())
		return 1;
	else return 0;
}

