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
#include "wx/timer.h"

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
class MyDocument;

#if defined(__WXOSX_COCOA__)
#define wxTOGGLEBUTTON_STYLE wxBORDER_SIMPLE
#else
#define wxTOGGLEBUTTON_STYLE 0
#endif

enum {
	myMenuID_MyFirstMenuItem = 100,
	myMenuID_Import = 101,
	myMenuID_Export = 102,
	myMenuID_ExportGraphic = 103,
	myMenuID_SelectFragment = 104,
	myMenuID_SelectReverse = 105,
	myMenuID_CreateNewAtom = 106,
	myMenuID_CreateNewParameter = 107,
	myMenuID_CreateNewVdwParameter = 108,
	myMenuID_CreateNewBondParameter = 109,
	myMenuID_CreateNewAngleParameter = 110,
	myMenuID_CreateNewDihedralParameter = 111,
	myMenuID_CreateNewImproperParameter = 112,
	myMenuID_CreateNewVdwPairParameter = 113,
	myMenuID_CreateNewVdwCutoffParameter = 114,
	myMenuID_CreatePiAnchor = 115,
	myMenuID_AddHydrogen = 120,
	myMenuID_AddHydrogenSp3 = 121,
	myMenuID_AddHydrogenSp2 = 122,
	myMenuID_AddHydrogenLinear = 123,
	myMenuID_AddHydrogenPyramidal = 124,
	myMenuID_AddHydrogenBent = 125,
	myMenuID_FitToScreen = 150,
	myMenuID_CenterSelection = 151,
/*	myMenuID_ShowUnitCell = 152,
	myMenuID_ShowPeriodicBox = 153,
	myMenuID_ShowHydrogens = 154,
	myMenuID_ShowDummyAtoms = 155,
	myMenuID_ShowExpandedAtoms = 156,
	myMenuID_ShowEllipsoids = 157,
	myMenuID_ShowRotationCenter = 158, */
//	myMenuID_ShowGraphite = 159,
//	myMenuID_LineMode = 160,
	myMenuID_ShowAllAtoms = 161,
	myMenuID_HideSelected = 162,
	myMenuID_HideUnselected = 163,
	myMenuID_HideReverse = 164,
	myMenuID_MolecularDynamics = 200,
	myMenuID_Minimize = 201,
	myMenuID_StopMDRun = 202,
	myMenuID_ViewGlobalParameters = 205,
	myMenuID_ViewParameterFilesList = 206,
	myMenuID_ExecuteScript = 300,
	myMenuID_OpenConsoleWindow = 301,
	myMenuID_EmptyConsoleWindow = 302,
	myMenuID_CustomScript = 303,
	myMenuID_MyLastMenuItem = 499,
	myMenuID_PredefinedFragment = 500,
	myMenuID_MyLastFragment = 699,
	myMenuID_Internal_CheckIfAllWindowsAreGone = 900
};

enum {
	myMenuIndex_File = 0,
	myMenuIndex_Edit = 1,
	myMenuIndex_Show = 2,
	myMenuIndex_MMMD = 3,
	myMenuIndex_Script = 4
};

//  Global Setting Keys
extern const char *gSettingQuitOnCloseLastWindow;

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

	static wxString FindResourcePath();
	static wxString DefaultSettingsPath();

	void LoadDefaultSettings();
	void SaveDefaultSettings();
	void SetDefaultSetting(const wxString& key, const wxString& value);
	wxString& GetDefaultSetting(const wxString& key);

	int LookupScriptMenu(const char *title);
	int RegisterScriptMenu(const char *title);
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

	void OnEndProcess(wxProcessEvent &event);
	int CallSubProcess(const char *cmdline, const char *procname, int (*callback)(void *) = NULL, void *callback_data = NULL, FILE *fpout = NULL, FILE *fperr = NULL);

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
	
	void EnableTimerForDocument(MyDocument *doc);
	void DisableTimerForDocument(MyDocument *doc);
	void TimerInvoked(wxTimerEvent &event);

	void CheckIfAllWindowsAreGoneHandler(wxCommandEvent &event);
	void CheckIfAllWindowsAreGone(wxTopLevelWindow *frame);

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
	char **scriptMenuTitles;
	int *scriptMenuPositions;
	bool scriptMenuModifiedEventPosted;

	int m_CountNamedFragments;
	char **m_NamedFragments;

	int m_CountTimerDocs;
	MyDocument **m_TimerDocs;
	wxTimer *m_Timer;

	wxString *m_pendingFilesToOpen;  /*  Files to be processed by OnOpenFilesByIPC()  */

	wxTopLevelWindow *m_frameToBeDestroyed;   /*  Used in CheckIfAllWindowsAreGone()  */
	
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
class MyFrame: public wxDocParentFrame
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
