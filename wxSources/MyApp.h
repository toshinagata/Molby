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

#if __WXMSW__
#include "wx/snglinst.h"
#endif

#include "MyDocManager.h"

class MyDocManager;
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
	myMenuID_AddHydrogen = 105,
	myMenuID_AddHydrogenSp3 = 106,
	myMenuID_AddHydrogenSp2 = 107,
	myMenuID_AddHydrogenLinear = 108,
	myMenuID_AddHydrogenPyramidal = 109,
	myMenuID_AddHydrogenBent = 110,
	myMenuID_FitToScreen = 120,
	myMenuID_CenterSelection = 121,
	myMenuID_ShowUnitCell = 122,
	myMenuID_ShowPeriodicBox = 123,
	myMenuID_ShowHydrogens = 124,
	myMenuID_ShowDummyAtoms = 125,
	myMenuID_ShowExpandedAtoms = 126,
	myMenuID_ShowEllipsoids = 127,
	myMenuID_ShowRotationCenter = 128,
	myMenuID_ShowGraphite = 129,
	myMenuID_LineMode = 130,
	myMenuID_ShowAllAtoms = 131,
	myMenuID_HideSelected = 132,
	myMenuID_HideUnselected = 133,
	myMenuID_HideReverse = 134,
	myMenuID_MolecularDynamics = 140,
	myMenuID_Minimize = 141,
	myMenuID_StopMDRun = 142,
//	myMenuID_ReadParameters = 143,
	myMenuID_ViewGlobalParameters = 144,
	myMenuID_ViewParameterFilesList = 145,
	myMenuID_DefinePeriodicBox = 146,
	myMenuID_ShowPeriodicImage = 147,
	myMenuID_PressureControl = 148,
	myMenuID_DefineSymmetry = 149,
	myMenuID_ExpandBySymmetry = 150,
	myMenuID_MDTools = 151,
	myMenuID_RunAntechamber = 152,
	myMenuID_RunResp = 153,
	myMenuID_CreateSanderInput = 154,
	myMenuID_CreateGamessInput = 160,
	myMenuID_CreateMOCube = 161,
	myMenuID_ExecuteScript = 200,
	myMenuID_OpenConsoleWindow = 201,
	myMenuID_CustomScript = 202,
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

	int AppendConsoleMessage(const char *mes);
	void FlushConsoleMessage();
	void SetConsoleColor(int color);

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

	void RegisterScriptMenu(const char *cmd, const char *title);
	void UpdateScriptMenu(wxMenuBar *mbar);
	void OnScriptMenuModified(wxCommandEvent& event);
	void OnScriptMenuSelected(wxCommandEvent& event);
	void OnUpdateUI(wxUpdateUIEvent &event);
	void OnExecuteScript(wxCommandEvent &event);
	void OnOpenConsoleWindow(wxCommandEvent &event);
	void OnViewGlobalParameters(wxCommandEvent &event);
	void OnViewParameterFilesList(wxCommandEvent &event);

	void OnEndProcess(wxProcessEvent &event);
	int CallSubProcess(const char *cmdline, const char *procname);

	void OnActivate(wxActivateEvent &event);

	MyListCtrl *GetGlobalParameterListCtrl();
#if __WXMAC__
	virtual void MacNewFile();
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
	
#if __WXMSW__
	wxSingleInstanceChecker *m_checker;
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
