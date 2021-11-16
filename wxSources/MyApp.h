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
	myMenuID_MyFirstMenuItem = 12000,
	myMenuID_Import = 12001,
	myMenuID_Export = 12002,
	myMenuID_ExportGraphic = 12003,
	myMenuID_SelectFragment = 12004,
	myMenuID_SelectReverse = 12005,

	myMenuID_CreateNewAtom = 12010,
	myMenuID_CreateNewParameter = 12011,
	myMenuID_CreateNewVdwParameter = 12012,
	myMenuID_CreateNewBondParameter = 12013,
	myMenuID_CreateNewAngleParameter = 12014,
	myMenuID_CreateNewDihedralParameter = 12015,
	myMenuID_CreateNewImproperParameter = 12016,
	myMenuID_CreateNewVdwPairParameter = 12017,
	myMenuID_CreateNewVdwCutoffParameter = 12018,
	myMenuID_CreatePiAnchor = 12019,

	myMenuID_AddHydrogen = 12030,
	myMenuID_AddHydrogenSp3 = 12031,
	myMenuID_AddHydrogenSp2 = 12032,
	myMenuID_AddHydrogenLinear = 12033,
	myMenuID_AddHydrogenPyramidal = 12034,
	myMenuID_AddHydrogenBent = 12035,

	myMenuID_FitToScreen = 12050,
	myMenuID_CenterSelection = 12051,

	myMenuID_ShowAllAtoms = 12060,
	myMenuID_HideSelected = 12061,
	myMenuID_HideUnselected = 12062,
	myMenuID_HideReverse = 12063,

	myMenuID_MolecularDynamics = 12070,
	myMenuID_Minimize = 12071,
	myMenuID_StopMDRun = 12072,
	
	myMenuID_ViewGlobalParameters = 12080,
	myMenuID_ViewParameterFilesList = 12081,
	
	myMenuID_BringAllWindowsToFront = 12090,
	
	myMenuID_ExecuteScript = 12100,
	myMenuID_OpenConsoleWindow = 12101,
	myMenuID_EmptyConsoleWindow = 12102,

	myMenuID_PredefinedFragment = 12200,
	myMenuID_MyLastFragment = 12999,
	
	//  The ID of script menu has "1000 * depth" offset
	myMenuID_CustomScript = 13000,
	
	myMenuID_MyLastMenuItem = 29999,

	myMenuID_Internal_CheckIfAllWindowsAreGone = 30000
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

//  A support class for wxProcess
//  When the process terminates, the exit status is kept inside the object
class wxBetterProcess : public wxProcess
{
public:
    wxBetterProcess(wxEvtHandler *parent, int id) : wxProcess(parent, id)
    {
        m_status = 0;
        m_terminated = false;
        m_killSignal = wxSIGNONE;
    }
    wxBetterProcess(int flags) : wxProcess(flags)
    {
        m_status = 0;
        m_terminated = false;
        m_killSignal = wxSIGNONE;
    }
    virtual ~wxBetterProcess() {}
    virtual void OnTerminate(int pid, int status);
    wxKillError KillProcess(wxSignal sig = wxSIGTERM, int flags = wxKILL_NOCHILDREN);
    int PutLine(wxString str);
    int GetLine(wxString &outStr);
    int GetLineSub(wxString &outStr, wxInputStream *stream, wxMemoryBuffer &mbuf);
    int GetErrorLine(wxString &outStr);
    void CloseOutput();
    bool IsTerminated() { return m_terminated; }
    int GetStatus() { return m_status; }
    int GetKillSignal() { return m_killSignal; }
protected:
    bool m_terminated;
    int m_status;
    int m_killSignal;
    wxMemoryBuffer m_stdout;
    wxMemoryBuffer m_stderr;
    wxMemoryBuffer m_stdin;
};


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
	// ProgressFrame *GetProgressFrame() { return m_progressFrame; }
    
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
	void OnBringAllWindowsToFront(wxCommandEvent &event);
	
//	void OnEndProcess(wxProcessEvent &event);
	int CallSubProcess(const char *cmdline, const char *procname, int (*callback)(void *) = NULL, void *callback_data = NULL, FILE *fpout = NULL, FILE *fperr = NULL, int *exitstatus_p = NULL, int *pid_p = NULL);

	void OnActivate(wxActivateEvent &event);

	void RequestOpenFilesByEvent(wxString& files);
	void OnOpenFilesByEvent(wxCommandEvent& event);
	
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

    void OnHelp(wxCommandEvent& event);

protected:
    MyDocManager* m_docManager;
	wxProgressDialog *m_progressDialog;
//	bool m_progressCanceled;
//	int m_progressValue;
	MyStringHash m_defaultSettings;

	//  For CallSubProcess()
	wxBetterProcess *m_process;
//	bool m_processTerminated;
//	int m_processExitCode;

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

	wxString *m_pendingFilesToOpen;  /*  Files to be processed by OnOpenFilesByEvent()  */

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

// About dialog
class AboutDialog: public wxDialog
{
    DECLARE_CLASS(AboutDialog)
public:
    AboutDialog();
    void OnOKPressed(wxCommandEvent &event);
private:
    DECLARE_EVENT_TABLE()
};

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
