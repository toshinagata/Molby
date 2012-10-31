/*
 *  MyApp.h
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

#ifndef __MyApp_H__
#define __MyApp_H__

#include "wx/app.h"
#include "wx/mdi.h"
#include "wx/docview.h"
#include "wx/docmdi.h"
#include "wx/hashmap.h"
#include "wx/process.h"

#if defined(__WXMSW__)
#include "wx/snglinst.h"
#endif

#include "MyDocManager.h"

class MyDocManager;

#if defined(__WXMSW__)
class MyServer;
class MyClient;
#endif

class wxMenuBar;
class wxMenu;
class wxProgressDialog;

class ConsoleFrame;
class ProgressFrame;
class GlobalParameterFrame;
class GlobalParameterFilesFrame;
class MyListCtrl;

enum {
	myMenuID_MyFirstMenuItem = 100,
	myMenuID_Import = 101,
	myMenuID_Export = 102,
	myMenuID_SelectFragment = 103,
	myMenuID_SelectReverse = 104,
	myMenuID_CreateNewAtom = 105,
	myMenuID_CreateNewParameter = 106,
	myMenuID_CreateNewVdwParameter = 107,
	myMenuID_CreateNewBondParameter = 108,
	myMenuID_CreateNewAngleParameter = 109,
	myMenuID_CreateNewDihedralParameter = 110,
	myMenuID_CreateNewImproperParameter = 111,
	myMenuID_CreateNewVdwPairParameter = 112,
	myMenuID_CreateNewVdwCutoffParameter = 113,
	myMenuID_CreatePiAnchor = 114,
	myMenuID_AddHydrogen = 120,
	myMenuID_AddHydrogenSp3 = 121,
	myMenuID_AddHydrogenSp2 = 122,
	myMenuID_AddHydrogenLinear = 123,
	myMenuID_AddHydrogenPyramidal = 124,
	myMenuID_AddHydrogenBent = 125,
	myMenuID_FitToScreen = 150,
	myMenuID_CenterSelection = 151,
	myMenuID_ShowUnitCell = 152,
	myMenuID_ShowPeriodicBox = 153,
	myMenuID_ShowHydrogens = 154,
	myMenuID_ShowDummyAtoms = 155,
	myMenuID_ShowExpandedAtoms = 156,
	myMenuID_ShowEllipsoids = 157,
	myMenuID_ShowRotationCenter = 158,
	myMenuID_ShowGraphite = 159,
	myMenuID_LineMode = 160,
	myMenuID_ShowAllAtoms = 161,
	myMenuID_HideSelected = 162,
	myMenuID_HideUnselected = 163,
	myMenuID_HideReverse = 164,
	myMenuID_MolecularDynamics = 200,
	myMenuID_Minimize = 201,
	myMenuID_StopMDRun = 202,
//	myMenuID_ReadParameters = 143,
	myMenuID_ViewGlobalParameters = 204,
	myMenuID_ViewParameterFilesList = 205,
	myMenuID_DefinePeriodicBox = 206,
	myMenuID_ShowPeriodicImage = 207,
	myMenuID_PressureControl = 208,
	myMenuID_DefineSymmetry = 209,
	myMenuID_ExpandBySymmetry = 210,
	myMenuID_MDTools = 211,
	myMenuID_RunAntechamber = 212,
	myMenuID_RunResp = 213,
	myMenuID_CreateSanderInput = 214,
	myMenuID_ImportAmberLib = 215,
	myMenuID_ImportAmberFrcmod = 216,
	myMenuID_CreateGamessInput = 250,
	myMenuID_CreateMOCube = 251,
	myMenuID_ExecuteScript = 300,
	myMenuID_OpenConsoleWindow = 301,
	myMenuID_EmptyConsoleWindow = 302,
	myMenuID_CustomScript = 303,
	myMenuID_PredefinedFragment = 350,
	myMenuID_MyLastMenuItem = 499
};

enum {
	myMenuIndex_File = 0,
	myMenuIndex_Edit = 1,
	myMenuIndex_Show = 2,
	myMenuIndex_MMMD = 3,
	myMenuIndex_QChem = 4,
	myMenuIndex_Script = 5
};

WX_DECLARE_STRING_HASH_MAP( wxString, MyStringHash );

// Define a new application
class MyApp: public wxApp
{
  public:
    MyApp(void);
    bool OnInit(void);
    int OnExit(void);

	ConsoleFrame *GetConsoleFrame() { return consoleFrame; }
	
	void ShowProgressPanel(const char *mes);
	void HideProgressPanel();
	void SetProgressValue(double dval);
	void SetProgressMessage(const char *mes);
	int IsInterrupted();
	ProgressFrame *GetProgressFrame() { return m_progressFrame; }

    MyDocManager *DocManager() { return m_docManager; }

	wxMenuBar *CreateMenuBar(int kind, wxMenu **out_file_history_menu = NULL, wxMenu **out_edit_menu = NULL);

//	void OnReadParameters(wxCommandEvent& event);

	static wxString FindResourcePath();
	static wxString DefaultSettingsPath();

	void LoadDefaultSettings();
	void SaveDefaultSettings();
	void SetDefaultSetting(const wxString& key, const wxString& value);
	wxString& GetDefaultSetting(const wxString& key);

	int LookupScriptMenu(const char *title);
	void RegisterScriptMenu(const char *cmd, const char *title);
	void UpdateScriptMenu(wxMenuBar *mbar);
	void OnScriptMenuModified(wxCommandEvent& event);
	void OnScriptMenuSelected(wxCommandEvent& event);
	void UpdatePredefinedFragmentMenu(wxMenuBar *mbar);
	void OnFragmentMenuSelected(wxCommandEvent& event);
	
	void OnUpdateUI(wxUpdateUIEvent &event);
	void OnExecuteScript(wxCommandEvent &event);
	void OnOpenConsoleWindow(wxCommandEvent &event);
	void OnEmptyConsoleWindow(wxCommandEvent &event);
	void OnViewGlobalParameters(wxCommandEvent &event);
	void OnViewParameterFilesList(wxCommandEvent &event);

	void OnImportAmberLib(wxCommandEvent &event);

	void OnEndProcess(wxProcessEvent &event);
	int CallSubProcess(const char *cmdline, const char *procname);

	void OnActivate(wxActivateEvent &event);

	void RequestOpenFilesByIPC(wxString& files);
	void OnOpenFilesByIPC(wxCommandEvent& event);
	
	bool OnOpenFiles(const wxString &files);

	MyListCtrl *GetGlobalParameterListCtrl();
	
#if defined(__WXMAC__)
	virtual void MacNewFile();
	virtual void MacOpenFile(const wxString &fileName);
	virtual short MacHandleAEODoc(const WXEVENTREF event, WXEVENTREF WXUNUSED(reply));
#endif
	
protected:
    MyDocManager* m_docManager;
	ProgressFrame *m_progressFrame;
	bool m_progressCanceled;
	int m_progressValue;
	MyStringHash m_defaultSettings;

	bool m_processTerminated;
	int m_processExitCode;

	ConsoleFrame *consoleFrame;
	GlobalParameterFrame *parameterFrame;
	GlobalParameterFilesFrame *parameterFilesFrame;
	
	int countNonCustomScriptMenu;
	int countScriptMenu;
	char **scriptMenuCommands;
	char **scriptMenuTitles;
	bool scriptMenuModifiedEventPosted;

	int m_CountNamedFragments;
	char **m_NamedFragments;

	wxString *m_pendingFilesToOpen;  /*  Files to be processed by OnOpenFilesByIPC()  */

#if defined(__WXMSW__)
public:
	wxSingleInstanceChecker *m_checker;
	wxString *m_ipcServiceName;
	MyServer *m_server;
	MyClient *m_client;
#endif
	
private:
	DECLARE_EVENT_TABLE()
};

DECLARE_APP(MyApp)

// Define a new frame
class MyFrame: public wxDocMDIParentFrame
{
	DECLARE_CLASS(MyFrame)
public:
	wxMenu *editMenu;
  
	MyFrame(wxDocManager *manager, wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size, long type);

	void OnAbout(wxCommandEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

extern MyFrame *GetMainFrame(void);
extern bool singleWindowMode;

#endif  /* __MyApp_H__ */
