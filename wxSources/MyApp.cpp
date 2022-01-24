/*
 *  MyApp.cpp
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

#if !wxUSE_MDI_ARCHITECTURE
#error "You should have MDI architecture enabled in your wxWidgets installation."
#endif

#include "wx/filename.h"
#include "wx/progdlg.h"
#include "wx/sysopt.h"
#include "wx/regex.h"
#include "wx/stdpaths.h"
#include "wx/textfile.h"
#include "wx/process.h"
#include "wx/utils.h"
#include "wx/sound.h"
#include "wx/time.h"

#include "MyApp.h"
#include "MyDocument.h"
#include "MoleculeView.h"
#include "ConsoleFrame.h"
#include "ProgressFrame.h"
#include "GlobalParameterFrame.h"
#include "GlobalParameterFilesFrame.h"
#include "RubyDialogFrame.h"
#include "MyMBConv.h"

#if defined(__WXMSW__)
#include "MyIPCSupport.h"
#endif

#include "../MolLib/MolLib.h"
#include "../MolLib/Ruby_bind/Molby_extern.h"
#include "../MolLib/Missing.h"

#include <wchar.h>
#include <stdio.h>

#if defined(__WXMAC__)
#include <CoreFoundation/CoreFoundation.h>
#undef T_DATA
#include <Carbon/Carbon.h>
// #include "wx/mac/carbon/private.h"
#include <sys/wait.h>  /* for waitpid()  */
#endif

#pragma mark ====== MyApp ======

const char *gSettingQuitOnCloseLastWindow = "quit_on_close_last_window";

MyFrame *frame = (MyFrame *) NULL;

bool gInitCompleted = false;

int gSuppressConsole = 0;  //  Non-zero if console output should be suppressed in non-GUI mode
int gUseGUI = 1;

IMPLEMENT_APP(MyApp)

//IMPLEMENT_CLASS(MyApp, wxApp)

BEGIN_EVENT_TABLE(MyApp, wxApp)
	EVT_COMMAND(MyDocumentEvent_scriptMenuModified, MyDocumentEvent, MyApp::OnScriptMenuModified)
	EVT_COMMAND(MyDocumentEvent_openFilesByEvent, MyDocumentEvent, MyApp::OnOpenFilesByEvent)
	EVT_UPDATE_UI_RANGE(myMenuID_MyFirstMenuItem, myMenuID_MyLastMenuItem, MyApp::OnUpdateUI)
	EVT_MENU(myMenuID_ExecuteScript, MyApp::OnExecuteScript)
	EVT_MENU(myMenuID_OpenConsoleWindow, MyApp::OnOpenConsoleWindow)
	EVT_MENU(myMenuID_EmptyConsoleWindow, MyApp::OnEmptyConsoleWindow)
	EVT_MENU(myMenuID_ViewGlobalParameters, MyApp::OnViewGlobalParameters)
	EVT_MENU(myMenuID_ViewParameterFilesList, MyApp::OnViewParameterFilesList)
	EVT_MENU(myMenuID_BringAllWindowsToFront, MyApp::OnBringAllWindowsToFront)
    EVT_MENU(wxID_HELP, MyApp::OnHelp)
#if defined(__WXMAC__)
	EVT_ACTIVATE(MyApp::OnActivate)
#endif
//	EVT_END_PROCESS(-1, MyApp::OnEndProcess)
	EVT_TIMER(-1, MyApp::TimerInvoked)
	EVT_COMMAND(myMenuID_Internal_CheckIfAllWindowsAreGone, MyDocumentEvent, MyApp::CheckIfAllWindowsAreGoneHandler)
END_EVENT_TABLE()

//  Find the path of the directory where the relevant resources are to be found.
//  Mac: the "Resources" directory in the application bundle.
//  Windows: the directory in which the application executable is located.
//  UNIX: ?
wxString
MyApp::FindResourcePath()
{
#if defined(__WXMAC__)
	CFBundleRef mainBundle = CFBundleGetMainBundle();
	CFURLRef ref = CFBundleCopyResourcesDirectoryURL(mainBundle);
	if (ref != NULL) {
		UInt8 buffer[256];
		if (CFURLGetFileSystemRepresentation(ref, true, buffer, sizeof buffer)) {
			wxString dirname((const char *)buffer, WX_DEFAULT_CONV);
			CFRelease(ref);
			return dirname;
		}
		CFRelease(ref);
	}
	return wxEmptyString;
#elif defined(__WXMSW__)
    wxString str;
	wxString argv0 = wxTheApp->argv[0];
	//  Fix dosish path (when invoked from MSYS console, the path may be unix-like)
	//  Note: absolute paths like /c/Molby/... (== c:\Molby\...) is not supported
	{
		char *p = strdup(argv0.mb_str(wxConvFile));
		fix_dosish_path(p);
		wxString argv0_fixed(p, wxConvFile);
		argv0 = argv0_fixed;
	}
	//  Is it an absolute path?
    if (wxIsAbsolutePath(argv0)) {
        return wxPathOnly(argv0);
    } else {
        //  Is it a relative path?
        wxString currentDir = wxGetCwd();
        if (currentDir.Last() != wxFILE_SEP_PATH)
            currentDir += wxFILE_SEP_PATH;		
        str = currentDir + argv0;
        if (wxFileExists(str))
            return wxPathOnly(str);
    }
	//  Search PATH
    wxPathList pathList;
    pathList.AddEnvList(wxT("PATH"));
    str = pathList.FindAbsoluteValidPath(argv0);
    if (!str.IsEmpty())
        return wxPathOnly(str);
    return wxEmptyString;
#else
#error "FindResourcePath is not defined for UNIXes."
#endif
}

MyApp::MyApp(void)
{
    m_docManager = NULL;
	m_progressDialog = NULL;
	m_process = NULL;
//	m_processTerminated = false;
//	m_processExitCode = 0;
	countScriptMenu = 0;
	scriptMenuTitles = NULL;
	scriptMenuPositions = NULL;
	scriptMenuModifiedEventPosted = false;
	parameterFrame = NULL;
	parameterFilesFrame = NULL;
	consoleFrame = NULL;
	m_CountNamedFragments = 0;
	m_NamedFragments = (char **)(-1);  /*  Will be set to NULL after Ruby interpreter is initialized  */
	m_pendingFilesToOpen = NULL;
	m_CountTimerDocs = 0;
	m_TimerDocs = NULL;
	m_Timer = NULL;

#if defined(__WXMSW__)
	m_checker = NULL;
	m_ipcServiceName = NULL;
	m_server = NULL;
	m_client = NULL;
#endif
}

//  We override Initialize() instead of OnInit, because wxAppBase::Initialize() calls OnInitGui(), which
//  makes the program run as a GUI application.
//  So, we intercept here, and go directly to the execution in the batch mode.
//  Otherwise, we call the inherited version of Initialize() and the program will run as a normal application.
bool MyApp::Initialize(int& argc, wxChar **argv)
{
    //  Called with a batch mode?
    if (argc > 1 && wcscmp(argv[1], L"-b") == 0) {

        //  Disable any wxLog functionality (otherwise ::exit() may crash)
        wxLog::EnableLogging(false);

        gUseGUI = 0;
        gSuppressConsole = 1;
        
        if (argc > 2 && wcscmp(argv[2], L"-v") == 0)
            gSuppressConsole = 0;

        //  We need these parameters in FindResourcePath(), so we assign them here
        this->argc = argc;
        this->argv = argv;

        static const char fname[] = "startup.rb";
        wxString dirname = FindResourcePath();
        
        dirname += wxFILE_SEP_PATH;
        dirname += wxT("Scripts");
        wxString cwd = wxGetCwd();
        wxSetWorkingDirectory(dirname);
        
        wxString fnamestr(fname, wxConvFile);
        Molby_startup(wxFileExists(fnamestr) ? fname : NULL, (const char *)dirname.mb_str(wxConvFile));
        
        wxSetWorkingDirectory(cwd);
        
        //  Build ARGV
        int c = (gSuppressConsole ? 2 : 3);
        int i = 1;
        int status;
        if (c >= argc) {
            if (gSuppressConsole) {
                fprintf(stderr, "The script is not given\n");
                exit(1);
            } else exit(0);  //  Show startup message and exit
        }
        wxString argv_script;
        while (i + c < argc) {
            wxString arg(argv[i + c]);
            arg.Replace(wxT("\'"), wxT("\\\'"));
            argv_script += wxString::Format(wxT("ARGV[%d] = \'"), i - 1);
            argv_script += arg;
            argv_script += wxT("\'\n");
            i++;
        }
        gSuppressConsole = 0;  //  Console output is no longer suppressed (startup is done)
        status = Molby_loadScript(argv_script.mb_str(wxConvFile), 0);
        if (status == 0) {
            wxString arg2(argv[c]);
            status = Molby_loadScript(arg2.mb_str(wxConvFile), 1);
        }
        if (status != 0) {
            Ruby_showError(status);
        }
        //  Force exit
        ::exit(status);
    } else {
        //  Call the inherited version
        return wxApp::Initialize(argc, argv);
    }
}

bool MyApp::OnInit(void)
{
	//  Set defaults
#ifdef __WXMAC__
	wxSystemOptions::SetOption(wxT("mac.listctrl.always_use_generic"), 1);
	wxSystemOptions::SetOption(wxT("osx.openfiledialog.always-show-types"), 1);
#endif

#if __WXMSW__
	{
		//  Check if the same application is already running
		char *buf, *p;
		asprintf(&buf, "Molby-%s", (const char *)wxGetUserId().mb_str(WX_DEFAULT_CONV));
		wxString name(buf, WX_DEFAULT_CONV);
		m_ipcServiceName = new wxString(name);
		m_ipcServiceName->Prepend(wxT("IPC-"));
		free(buf);
		m_checker = new wxSingleInstanceChecker(name);
		if (m_checker->IsAnotherRunning()) {
			//  Make a connection with the other instance and ask for opening the file(s)
			if (argc > 1) {
				wxString files;
				wxConnectionBase *connection;
				int i;
				for (i = 1; i < argc; i++) {
					files.append(argv[i]);
					files.append(wxT("\n"));
				}
				m_client = new MyClient;
				connection = m_client->MakeConnection(wxT("localhost"), *m_ipcServiceName, MOLBY_IPC_TOPIC);
				if (connection == NULL) {
					wxLogError(wxT("Molby is already running; please shut it down and retry"));
					delete m_client;
					return false;
				}
				connection->Execute(files);
				delete m_client;
			}
			return false;
		} else {
			m_server = new MyServer;
			if (m_server->Create(*m_ipcServiceName) == false) {
				delete m_server;
				m_server = NULL;
			}
		}
	}
#endif
	
	// Create a document manager
	m_docManager = new MyDocManager;

	// Create templates relating drawing documents to their views
	new wxDocTemplate(m_docManager, _T("Molby Structure File"), _T("*.mbsf"), _T(""), _T("mbsf"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Protein Structure File"), _T("*.psf"), _T(""), _T("psf"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Protein Data Bank File"), _T("*.pdb"), _T(""), _T("pdb"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Gaussian Input File"), _T("*.com;*.gjf"), _T(""), _T("com"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Gaussian/GAMESS Log File"), _T("*.out;*.log"), _T(""), _T("out"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Gaussian Checkpoint File"), _T("*.fchk;*.fch"), _T(""), _T("fchk"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("GAMESS Input File"), _T("*.inp"), _T(""), _T("inp"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("GAMESS DAT File"), _T("*.dat"), _T(""), _T("dat"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("ORTEP Input File"), _T("*.tep"), _T(""), _T("tep"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("SHELX Input File"), _T("*.ins;*.res"), _T(""), _T("ins"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Crystallographic Information File"), _T("*.cif"), _T(""), _T("cif"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Cartesian"), _T("*.xyz"), _T(""), _T("xyz"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Z Matrix"), _T("*.zmat"), _T(""), _T("zmat"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Any Molecule"), _T("*.*"), _T(""), _T(""), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));

    // Init image handlers
    MyAppCallback_initImageHandlers();

	// Create the main frame window
	frame = new MyFrame((wxDocManager *) m_docManager, (wxFrame *) NULL,
                      _T("Molby"), wxPoint(0, 0), wxSize(800, 600),
                      wxDEFAULT_FRAME_STYLE | wxNO_FULL_REPAINT_ON_RESIZE);

	// Give it an icon (this is ignored in MDI mode: uses resources)
#ifdef __WXMSW__
	frame->SetIcon(wxIcon(_T("doc")));
#endif
#ifdef __X__
	frame->SetIcon(wxIcon(_T("doc.xbm")));
#endif
	
	wxMenuBar *menu_bar = CreateMenuBar(0);
	
#ifdef __WXMAC__
	wxMenuBar::MacSetCommonMenuBar(menu_bar);
#endif

	// Associate the menu bar with the frame
	frame->SetMenuBar(menu_bar);

	frame->Centre(wxBOTH);

#if defined(__WXMAC__) || defined(__WXMSW__)
	frame->Move(-10000, -10000);  //  Set invisible
	frame->Show(false);
#else
	frame->Show(true);
#endif
	
	SetTopWindow(frame);
	
	//  Load default settings from the preference file
	LoadDefaultSettings();
	
	//  Create a console window
	consoleFrame = ConsoleFrame::CreateConsoleFrame(frame);
	consoleFrame->Show(true);

	/*  Initialize Ruby interpreter with the startup script  */
	/*  (Also read startup information)  */
	{
		extern int gRevisionNumber;
		static const char fname[] = "startup.rb";
		wxString dirname = FindResourcePath();
	
		dirname += wxFILE_SEP_PATH;
		dirname += wxT("Scripts");
	/*	wxString cwd = wxGetCwd(); */
		wxSetWorkingDirectory(dirname);

		wxString fnamestr(fname, wxConvFile);
		Molby_startup(wxFileExists(fnamestr) ? fname : NULL, (const char *)dirname.mb_str(wxConvFile));
		
		{
			wxString docHome(MyAppCallback_getDocumentHomeDir(), wxConvFile);
			wxSetWorkingDirectory(docHome);
			
		}

		/*  Pasteboard type strings (includes the revision number)  */
		asprintf(&gMoleculePasteboardType, "Molecule_%d", gRevisionNumber);
		asprintf(&gParameterPasteboardType, "Parameter_%d", gRevisionNumber);

		MyAppCallback_showScriptMessage("%% ");
		
		/*  Build the predefined fragments menu  */
		m_NamedFragments = NULL;
	//	UpdatePredefinedFragmentMenu(GetMainFrame()->GetMenuBar());
		UpdatePredefinedFragmentMenu(consoleFrame->GetMenuBar());

	}
	
	/*  Open given files: Ruby script is executed, other files are opened as a document  */
	if (argc == 1) {
#if defined(__WXMSW__)
		m_docManager->CreateDocument(wxEmptyString, wxDOC_NEW);
#endif
	} else {
		int i;
		wxString files;
		for (i = 1; i < argc; i++) {
			files.append(argv[i]);
			files.append(wxT("\n"));
		}
		RequestOpenFilesByEvent(files);
	}
	
	gInitCompleted = true;
	
	return true;
}

//  Create Menu Bars
//  kind == 0: main menu
//  kind == 1: molecule window
//  kind == 2: console window
//  kind == 3: Ruby dialog (non-modal)
wxMenuBar *
MyApp::CreateMenuBar(int kind, wxMenu **out_file_history_menu, wxMenu **out_edit_menu)
{

#if __WXMSW__
	if (kind == 3) {

		//  Simplified menu
		wxMenu *file_menu = new wxMenu;
		file_menu->Append(wxID_CLOSE, _T("&Close\tCtrl-W"));

		wxMenu *edit_menu = new wxMenu;
		edit_menu->Append(wxID_UNDO, _T("&Undo\tCtrl-Z"));
		edit_menu->Append(wxID_REDO, _T("&Redo"));
		edit_menu->AppendSeparator();
		edit_menu->Append(wxID_CUT, _T("Cut\tCtrl-X"));
		edit_menu->Append(wxID_COPY, _T("Copy\tCtrl-C"));
		edit_menu->Append(wxID_PASTE, _T("Paste\tCtrl-V"));
		edit_menu->Append(wxID_CLEAR, _T("Clear"));
		edit_menu->AppendSeparator();
		edit_menu->Append(wxID_SELECTALL, _T("Select All\tCtrl-A"));
		
		wxMenu *help_menu = new wxMenu;
		help_menu->Append(wxID_ABOUT, _T("&About...\tF1"));
        help_menu->Append(wxID_HELP, _T("&Molby Help"));

		wxMenuBar *menu_bar = new wxMenuBar;
		
		menu_bar->Append(file_menu, _T("&File"));
		menu_bar->Append(edit_menu, _T("&Edit"));
		menu_bar->Append(help_menu, _T("&Help"));
		
		return menu_bar;
	}
#endif
	
	wxMenu *file_menu = new wxMenu;
	wxMenu *file_history_menu = NULL;

	file_menu->Append(wxID_NEW, _T("&New...\tCtrl-N"));
	file_menu->Append(wxID_OPEN, _T("&Open...\tCtrl-O"));
	if (out_file_history_menu != NULL) {
		file_history_menu = new wxMenu;
		*out_file_history_menu = file_history_menu;
		file_menu->AppendSubMenu(*out_file_history_menu, _T("Open Recent"));
		m_docManager->FileHistoryAddFilesToMenu(*out_file_history_menu);
		m_docManager->FileHistoryUseMenu(*out_file_history_menu);  //  Should be removed when menu is discarded
	}
	/*  Build "Open Predefined"  */
	wxMenu *predefined_menu = new wxMenu;
	file_menu->Append(myMenuID_PredefinedFragment, _T("Open Predefined"), predefined_menu);

	file_menu->AppendSeparator();
	file_menu->Append(wxID_CLOSE, _T("&Close\tCtrl-W"));
	file_menu->Append(wxID_SAVE, _T("&Save\tCtrl-S"));
	file_menu->Append(wxID_SAVEAS, _T("Save &As..."));	
	file_menu->Append(wxID_REVERT, _T("Revert to Saved"));	
	
	file_menu->AppendSeparator();
	file_menu->Append(myMenuID_Import, _T("Import..."));	
	file_menu->Append(myMenuID_Export, _T("Export..."));	
	file_menu->Append(myMenuID_ExportGraphic, _T("Export Graphic..."));	
	
	file_menu->AppendSeparator();
	file_menu->Append(wxID_PRINT, _T("&Print...\tCtrl-P"));
	file_menu->Append(wxID_PRINT_SETUP, _T("Print &Setup..."));
	file_menu->Append(wxID_PREVIEW, _T("Print Pre&view"));
	
	file_menu->AppendSeparator();
#if defined(__WXMAC__)
	file_menu->Append(wxID_EXIT, _T("E&xit\tCtrl-Q"));
#else
	file_menu->Append(wxID_EXIT, _T("E&xit\tAlt-X"));
#endif

	wxMenu *edit_menu = new wxMenu;
	edit_menu->Append(wxID_UNDO, _T("&Undo\tCtrl-Z"));
	edit_menu->Append(wxID_REDO, _T("&Redo"));
	edit_menu->AppendSeparator();
	edit_menu->Append(wxID_CUT, _T("Cut\tCtrl-X"));
	edit_menu->Append(wxID_COPY, _T("Copy\tCtrl-C"));
	edit_menu->Append(wxID_PASTE, _T("Paste\tCtrl-V"));
	edit_menu->Append(wxID_CLEAR, _T("Clear"));
	edit_menu->AppendSeparator();
	edit_menu->Append(wxID_SELECTALL, _T("Select All\tCtrl-A"));
	edit_menu->Append(myMenuID_SelectFragment, _T("Select Fragment\tCtrl-F"));
	edit_menu->Append(myMenuID_SelectReverse, _T("Select Reverse"));
	edit_menu->AppendSeparator();
	wxMenu *create_parameter_menu = new wxMenu;
	create_parameter_menu->Append(myMenuID_CreateNewVdwParameter, _T("Vdw"));
	create_parameter_menu->Append(myMenuID_CreateNewBondParameter, _T("Bond"));
	create_parameter_menu->Append(myMenuID_CreateNewAngleParameter, _T("Angle"));
	create_parameter_menu->Append(myMenuID_CreateNewDihedralParameter, _T("Dihedral"));
	create_parameter_menu->Append(myMenuID_CreateNewImproperParameter, _T("Improper"));
	create_parameter_menu->Append(myMenuID_CreateNewVdwPairParameter, _T("Vdw Pair"));
	create_parameter_menu->Append(myMenuID_CreateNewVdwCutoffParameter, _T("Vdw Cutoff"));
	edit_menu->Append(myMenuID_CreateNewAtom, _T("Create New Atom\tCtrl-I"));
	edit_menu->Append(myMenuID_CreateNewParameter, _T("Create New Parameter"), create_parameter_menu);
	edit_menu->Append(myMenuID_CreatePiAnchor, _T("Create Pi Anchor"));
	edit_menu->AppendSeparator();
	wxMenu *add_hydrogen_menu = new wxMenu;
	add_hydrogen_menu->Append(myMenuID_AddHydrogenSp3, _T("Tetrahedral sp3"));
	add_hydrogen_menu->Append(myMenuID_AddHydrogenSp2, _T("Trigonal sp2"));
	add_hydrogen_menu->Append(myMenuID_AddHydrogenLinear, _T("Linear sp"));
	add_hydrogen_menu->Append(myMenuID_AddHydrogenPyramidal, _T("Pyramidal (like NH2)"));
	add_hydrogen_menu->Append(myMenuID_AddHydrogenBent, _T("Bent (like OH)"));
	edit_menu->Append(myMenuID_AddHydrogen, _T("Add Hydrogen"), add_hydrogen_menu);
	
	if (out_edit_menu != NULL)
		*out_edit_menu = edit_menu;	// Should be associated with the command processor if available
	
	wxMenu *view_menu = new wxMenu;
	view_menu->Append(myMenuID_FitToScreen, _T("Fit To Screen\tCtrl-T"));
	view_menu->Append(myMenuID_CenterSelection, _T("Center Selection"));
/*	view_menu->AppendSeparator();
	view_menu->Append(myMenuID_ShowUnitCell, _T("Show Unit Cell"), _T(""), wxITEM_CHECK);
	view_menu->Append(myMenuID_ShowHydrogens, _T("Show Hydrogen Atoms"), _T(""), wxITEM_CHECK);
	view_menu->Append(myMenuID_ShowDummyAtoms, _T("Show Dummy Atoms"), _T(""), wxITEM_CHECK);
	view_menu->Append(myMenuID_ShowExpandedAtoms, _T("Show Expanded Atoms"), _T(""), wxITEM_CHECK);
	view_menu->Append(myMenuID_ShowEllipsoids, _T("Show Ellipsoids"), _T(""), wxITEM_CHECK);
	view_menu->Append(myMenuID_ShowRotationCenter, _T("Show Rotation Center"), _T(""), wxITEM_CHECK); */
	view_menu->AppendSeparator();
	view_menu->Append(myMenuID_HideSelected, _T("Hide Selected"), _T(""));
	view_menu->Append(myMenuID_HideUnselected, _T("Hide Unselected"), _T(""));
	view_menu->Append(myMenuID_HideReverse, _T("Hide Reverse"), _T(""));
	view_menu->Append(myMenuID_ShowAllAtoms, _T("Show All Atoms"), _T(""));
//	view_menu->AppendSeparator();
//	view_menu->Append(myMenuID_ShowGraphite, _T("Show Graphite..."));
//	view_menu->AppendSeparator();
//	view_menu->Append(myMenuID_LineMode, _T("Line Mode"), _T(""), wxITEM_CHECK);

	wxMenu *md_menu = new wxMenu;
	md_menu->Append(myMenuID_MolecularDynamics, _T("Molecular Dynamics..."));
	md_menu->Append(myMenuID_Minimize, _T("Minimize..."));
	md_menu->Append(myMenuID_StopMDRun, _T("Stop\tCtrl-."));
	md_menu->AppendSeparator();
	md_menu->Append(myMenuID_ViewGlobalParameters, _T("View Global Parameters..."));
	md_menu->Append(myMenuID_ViewParameterFilesList, _T("Load/Unload Global Parameters..."));
	
	wxMenu *script_menu = new wxMenu;
	script_menu->Append(myMenuID_ExecuteScript, _T("Execute Script..."));
	script_menu->Append(myMenuID_OpenConsoleWindow, _T("Open Console Window"));
	script_menu->Append(myMenuID_EmptyConsoleWindow, _T("Empty Console Window"));
	script_menu->AppendSeparator();
	countNonCustomScriptMenu = script_menu->GetMenuItemCount();

	wxMenu *help_menu = new wxMenu;
	help_menu->Append(wxID_ABOUT, _T("&About...\tF1"));
    help_menu->Append(wxID_HELP, _T("&Molby Help"));

	wxMenuBar *menu_bar = new wxMenuBar;
	
	menu_bar->Append(file_menu, _T("&File"));
	menu_bar->Append(edit_menu, _T("&Edit"));
	menu_bar->Append(view_menu, _T("View"));
	menu_bar->Append(md_menu, _T("MM/MD"));
	menu_bar->Append(script_menu, _T("&Script"));
	
#if defined(__WXMAC__)
	wxMenu *window_menu = new wxMenu;
	window_menu->Append(myMenuID_BringAllWindowsToFront, _T("Bring All to Front"));
	window_menu->AppendSeparator();
	menu_bar->Append(window_menu, _T("Window"));
#endif
	
	menu_bar->Append(help_menu, _T("&Help"));
	
	UpdateScriptMenu(menu_bar);
	if (m_NamedFragments != (char **)(-1))
		UpdatePredefinedFragmentMenu(menu_bar);
	
	return menu_bar;
}

#if __WXMAC__
/*  When the application is launched without any documents, an empty document is opened.
    This should be implemented by overriding this special method; parsing argc/argv does
    not work, because the list of files is passed through an Apple Event.  */
void
MyApp::MacNewFile()
{
	if (m_docManager == NULL)
		return;  //  Initialization is not yet complete
	m_docManager->CreateDocument(_T(""), wxDOC_NEW);
}

void
MyApp::MacOpenFile(const wxString &fileName)
{
    wxString file(fileName);
	RequestOpenFilesByEvent(file);
}

/*  Open given files: instead of calling MacOpenFile() for each entry, build a file list
    and call MyApp::OnOpenFiles()  */
short
MyApp::MacHandleAEODoc(const WXEVENTREF event, WXEVENTREF WXUNUSED(reply))
{
    AEDescList docList;
    AEKeyword keywd;
    DescType returnedType;
    Size actualSize;
    long itemsInList;
    OSErr err;
    short i;
	
	return noErr;  /*  TODO: handle open Apple event  */
	
    err = AEGetParamDesc((AppleEvent *)event, keyDirectObject, typeAEList, &docList);
    if (err != noErr)
        return err;
	
    err = AECountItems(&docList, &itemsInList);
    if (err != noErr)
        return err;
	
    ProcessSerialNumber PSN ;
    PSN.highLongOfPSN = 0 ;
    PSN.lowLongOfPSN = kCurrentProcess ;
    SetFrontProcess( &PSN ) ;
	
    wxString fName, fNameList;
    FSRef theRef ;
	
    for (i = 1; i <= itemsInList; i++)
    {
        AEGetNthPtr(
					&docList, i, typeFSRef, &keywd, &returnedType,
					(Ptr)&theRef, sizeof(theRef), &actualSize);
    //    fName = wxMacFSRefToPath( &theRef ) ;
		fNameList.append(fName);
		fNameList.append(wxT("\n"));
    }
	
	OnOpenFiles(fNameList);
	
    return noErr;
}

#endif

int
MyApp::OnExit(void)
{
	SaveDefaultSettings();
    delete m_docManager;
#if __WXMSW__
	delete m_checker;
	delete m_server;
#endif
    return 0;
}

static void
sModifyMenuForFilterMode(wxMenuBar *mbar)
{
	int idx, i, n, id;
	wxMenu *menu;
	wxMenuItem *item;
	idx = mbar->FindMenu(wxT("Show"));
	if (idx != wxNOT_FOUND)
		delete mbar->Remove(idx);
	idx = mbar->FindMenu(wxT("MM/MD"));
	if (idx != wxNOT_FOUND)
		delete mbar->Remove(idx);
	idx = mbar->FindMenu(wxT("QChem"));
	if (idx != wxNOT_FOUND)
		delete mbar->Remove(idx);
	idx = mbar->FindMenu(wxT("Script"));
	if (idx != wxNOT_FOUND) {
		menu = mbar->GetMenu(idx);
		n = menu->GetMenuItemCount();
		for (i = n - 1; i >= 0; i--) {
			item = menu->FindItemByPosition(i);
			id = item->GetId();
			if (id != myMenuID_OpenConsoleWindow && id != myMenuID_EmptyConsoleWindow && id != myMenuID_ExecuteScript) {
				menu->Remove(item);
				delete item;
			}
		}
	}
	
	idx = mbar->FindMenu(wxT("Edit"));
	if (idx != wxNOT_FOUND) {
		menu = mbar->GetMenu(idx);
		n = menu->GetMenuItemCount();
		for (i = n - 1; i >= 0; i--) {
			item = menu->FindItemByPosition(i);
			id = item->GetId();
			if (id == wxID_SELECTALL)
				break;
			menu->Remove(item);
			delete item;
		}
	}
	
	idx = mbar->FindMenu(wxT("File"));
	if (idx != wxNOT_FOUND) {
		menu = mbar->GetMenu(idx);
		n = menu->GetMenuItemCount();
		for (i = n - 1; i >= 0; i--) {
			item = menu->FindItemByPosition(i);
			id = item->GetId();
			if (id != wxID_OPEN && id != wxID_EXIT) {
				menu->Remove(item);
				delete item;
			}
		}
	}
	
}

void
MyApp::ShowProgressPanel(const char *mes)
{
    wxString string((mes ? mes : ""), WX_DEFAULT_CONV);
    if (m_progressDialog == NULL) {
        m_progressDialog = new wxProgressDialog(wxT("Progress"), mes, 100, NULL, wxPD_APP_MODAL | wxPD_CAN_ABORT);
    }
/*
#if __WXMAC__
		{
			wxMenuBar *mbar = ((wxFrame *)GetTopWindow())->GetMenuBar();
			wxMenuItem *quitMenuItem = mbar->FindItem(wxID_EXIT);
			if (quitMenuItem != NULL)
				quitMenuItem->Enable(false);
			mbar->Enable(false);
		}
#endif
		m_progressFrame = new ProgressFrame(_T("Progress"), string);
		m_progressCanceled = false;
		m_progressValue = -1;
*/
}

void
MyApp::HideProgressPanel()
{
    if (m_progressDialog != NULL) {
        m_progressDialog->Hide();
        m_progressDialog->Destroy();
        m_progressDialog = NULL;
    }
/*
    if (m_progressFrame != NULL) {
		m_progressFrame->Hide();
		m_progressFrame->Destroy();
		m_progressFrame = NULL;
#if __WXMAC__
		{
			wxMenuBar *mbar = ((wxFrame *)GetTopWindow())->GetMenuBar();
			mbar->Enable(true);
			wxMenuItem *quitMenuItem = mbar->FindItem(wxID_EXIT);
			if (quitMenuItem != NULL)
				quitMenuItem->Enable(true);
		}
#endif
	}
*/
}

void
MyApp::SetProgressValue(double dval)
{
    if (m_progressDialog != NULL) {
        if (dval >= 0)
            m_progressDialog->Update((int)(dval * 100));
        else
            m_progressDialog->Pulse();
    }
}

void
MyApp::SetProgressMessage(const char *mes)
{
	if (m_progressDialog != NULL) {
		wxString string((mes ? mes : ""), WX_DEFAULT_CONV);
		m_progressDialog->Update(0, string);
	}
}

int
MyApp::IsInterrupted()
{
    if (m_progressDialog != NULL)
        return m_progressDialog->WasCancelled();
    else {
        if (::wxGetKeyState(WXK_ESCAPE))
            return 1;
        else return 0;
    }
}

void
MyApp::OnOpenConsoleWindow(wxCommandEvent& event)
{
	consoleFrame->Show(true);
	consoleFrame->Raise();
}

void
MyApp::OnEmptyConsoleWindow(wxCommandEvent& event)
{
	consoleFrame->EmptyBuffer();
}

void
MyApp::OnViewGlobalParameters(wxCommandEvent& event)
{
	if (parameterFrame == NULL) {
		parameterFrame = GlobalParameterFrame::CreateGlobalParameterFrame(GetMainFrame());
		MainView_createColumnsForTableAtIndex(NULL, kMainViewParameterTableIndex);
	}
	MainView_refreshTable(NULL);
	parameterFrame->Show(true);
	parameterFrame->Raise();
}

void
MyApp::OnViewParameterFilesList(wxCommandEvent &event)
{
	if (parameterFilesFrame == NULL) {
		parameterFilesFrame = GlobalParameterFilesFrame::CreateGlobalParameterFilesFrame(GetMainFrame());
	}
	parameterFilesFrame->Show(true);
	parameterFilesFrame->Raise();
}

void
MyApp::OnBringAllWindowsToFront(wxCommandEvent &event)
{
	int size = 0, n;
	wxWindowList::iterator iter;
	wxTopLevelWindow **wins;
	size = wxTopLevelWindows.size();
	if (size > 0) {
		wins = (wxTopLevelWindow **)calloc(sizeof(wxTopLevelWindow *), size);
		for (iter = wxTopLevelWindows.begin(), n = 0; iter != wxTopLevelWindows.end(); ++iter, ++n) {
			wins[n] = (wxTopLevelWindow *)(*iter);
		}
		for (n = 0; n < size; n++) {
			if (wins[n]->IsShown())
				wins[n]->Raise();
		}
	}
}

int
MyApp::LookupScriptMenu(const char *title)
{
	int i;
	if (title == NULL)
		return -1;
	for (i = 0; i < countScriptMenu; i++) {
		if (strcmp(title, scriptMenuTitles[i]) == 0)
			return i;
	}
	return -1;
}

int
MyApp::RegisterScriptMenu(const char *title)
{
	int i;

	//  Already registered? (If it is not a separator item)
	const char *p;
	p = strrchr(title, '\t');
	if (p == NULL)
		p = title;
	else p++;
	if (p[0] != 0 && p[0] != '-') {
		for (i = 0; i < countScriptMenu; i++) {
			if (strcmp(title, scriptMenuTitles[i]) == 0) {
				return i;
			}
		}
	}
	
	//  Not yet
	if (countScriptMenu == 0) {
		scriptMenuTitles = (char **)malloc(sizeof(char *));
		scriptMenuPositions = (int *)malloc(sizeof(int));
	} else {
		scriptMenuTitles = (char **)realloc(scriptMenuTitles, sizeof(char *) * (countScriptMenu + 1));
		scriptMenuPositions = (int *)realloc(scriptMenuPositions, sizeof(int) * (countScriptMenu + 1));
	}
	scriptMenuTitles[countScriptMenu] = strdup(title);
	scriptMenuPositions[countScriptMenu] = -1;
	countScriptMenu++;

	if (!scriptMenuModifiedEventPosted) {
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_scriptMenuModified);
		wxPostEvent(this, myEvent);	
		scriptMenuModifiedEventPosted = true;
	}
	
	return countScriptMenu - 1;
}

void
MyApp::UpdateScriptMenu(wxMenuBar *mbar)
{
	int i, mid;
	bool checkable;
	for (i = 0; i < countScriptMenu; i++) {
		//  Find Menu to be inserted
		const char *s = scriptMenuTitles[i], *p;
		wxMenu *menu = NULL;
		int depth = 1;
		while ((p = strchr(s, '\t')) != NULL) {
			//  Find existing menu
			wxString menuTitle(s, WX_DEFAULT_CONV, p - s);
			if (menu == NULL) {
				mid = mbar->FindMenu(menuTitle);
				if (mid == wxNOT_FOUND) {
					//  Create a new menu before "Script" menu
					wxMenu *newMenu = new wxMenu;
					int sid = mbar->FindMenu(_T("Script"));
					if (sid == wxNOT_FOUND) {
						mbar->Append(newMenu, menuTitle);
						sid = mbar->GetMenuCount() - 1;
					} else {
						mbar->Insert(sid, newMenu, menuTitle);
					}
					menu = mbar->GetMenu(sid);
				} else {
					menu = mbar->GetMenu(mid);
				}
			} else {
				mid = menu->FindItem(menuTitle);
				if (mid == wxNOT_FOUND) {
					//  Create a new menu at the end
					wxMenu *newMenu1 = new wxMenu;
					menu->Append(myMenuID_CustomScript + i + depth * 1000, menuTitle, newMenu1);
					menu = newMenu1;
				} else {
					menu = menu->FindItem(mid)->GetSubMenu();
				}
			}
			s = p + 1;
			depth++;
		}
		if (menu == NULL) {
			//  The new item should be under "Script" menu
			mid = mbar->FindMenu(_T("Script"));
			menu = mbar->GetMenu(mid);
		}
		if (*s == 0 || *s == '-') {
			//  Separator item
			wxMenuItem *sitem;
			int count = menu->GetMenuItemCount();
			if (scriptMenuPositions[i] >= 0 && scriptMenuPositions[i] < count) {
				sitem = menu->FindItemByPosition(scriptMenuPositions[i]);
				if (sitem != NULL && sitem->IsSeparator())
					continue;  //  Already present
			}
			if (count != 0 && !menu->FindItemByPosition(count - 1)->IsSeparator()) {
				menu->AppendSeparator();
				scriptMenuPositions[i] = count;
			}
			continue;
		}
		if (*s == '^') {
			//  Checkable item
			checkable = true;
			s++;
		} else checkable = false;
		//  Check if the item is already existing
		wxString itemTitle(s, WX_DEFAULT_CONV);
		wxMenu *omenu;
		wxMenuItem *item;
		item = mbar->FindItem(myMenuID_CustomScript + i, &omenu);
		if (item != NULL) {
			if (omenu == menu && item->GetItemLabel() == itemTitle) {
				//  The menu is already existing with correct position and title
				continue;
			}
			//  The menu title does not match; remove this menu item
			Disconnect(myMenuID_CustomScript + i, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MyApp::OnScriptMenuSelected), NULL, this);
			omenu->Remove(item);
			delete item;
		}
		//  Create menu item
		menu->Append(myMenuID_CustomScript + i, itemTitle, wxEmptyString, (checkable ? wxITEM_CHECK : wxITEM_NORMAL));
		scriptMenuPositions[i] = menu->GetMenuItemCount() - 1;
		Connect(myMenuID_CustomScript + i, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MyApp::OnScriptMenuSelected), NULL, this);
	}
}

void
MyApp::OnScriptMenuModified(wxCommandEvent& event)
{
	scriptMenuModifiedEventPosted = false;
	UpdateScriptMenu(consoleFrame->GetMenuBar());
	event.Skip();
}

void
MyApp::OnScriptMenuSelected(wxCommandEvent& event)
{
	MainView *mview;
	Molecule *mol;
	int index = event.GetId() - myMenuID_CustomScript;
	if (index < 0 || index >= countScriptMenu)
		return;
	mview = MainViewCallback_activeView();
	if (mview == NULL)
		mol = NULL;
	else mol = mview->mol;
	MolActionCreateAndPerform(NULL, SCRIPT_ACTION("iM"), "lambda { |n, m| $script_menu_commands[n].call(m) }", index, mol);
}

void
MyApp::UpdatePredefinedFragmentMenu(wxMenuBar *mbar)
{
	int i, n;
	wxMenuItem *fmenuItem = mbar->FindItem(myMenuID_PredefinedFragment);
	wxMenu *fmenu = (fmenuItem != NULL ? fmenuItem->GetSubMenu() : NULL);
	if (fmenu == NULL)
		return;
	if (m_NamedFragments == (char **)(-1))
		return;
	
	/*  Rebuild sNamedFragments array  */
	if (m_NamedFragments != NULL) {
		for (i = 0; i < m_CountNamedFragments; i++) {
			free(m_NamedFragments[i * 2]);
			free(m_NamedFragments[i * 2 + 1]);
		}
		free(m_NamedFragments);
	}
	m_NamedFragments = NULL;
	m_CountNamedFragments = 0;
	if (MolActionCreateAndPerform(NULL, SCRIPT_ACTION(";i"), "lambda { $named_fragments.length }", &n) != 0 || n <= 0)
		return;
	m_CountNamedFragments = n;
	m_NamedFragments = (char **)calloc(sizeof(char *), n * 2);
	for (i = 0; i < n; i++) {
		if (MolActionCreateAndPerform(NULL, SCRIPT_ACTION("i;s"), "lambda { |i| $named_fragments[i][0] }", i, &m_NamedFragments[i * 2]) != 0 ||
			MolActionCreateAndPerform(NULL, SCRIPT_ACTION("i;s"), "lambda { |i| $named_fragments[i][1] }", i, &m_NamedFragments[i * 2 + 1]) != 0)
			break;
	}
	if (i < n) {
		for (i = 0; i < m_CountNamedFragments; i++) {
			if (m_NamedFragments[i * 2] != NULL)
				free(m_NamedFragments[i * 2]);
			if (m_NamedFragments[i * 2 + 1] != NULL)
				free(m_NamedFragments[i * 2 + 1]);
		}
		free(m_NamedFragments);
		m_CountNamedFragments = 0;
		m_NamedFragments = NULL;
		return;
	}
	
	wxMenu *predefined_submenu = NULL;
	wxString stitle;
	int sn;
	for (i = sn = 0; i < m_CountNamedFragments; i++) {
		if (strcmp(m_NamedFragments[i * 2 + 1], "-") == 0) {
			if (predefined_submenu != NULL) {
				fmenu->Append(myMenuID_PredefinedFragment + 1 + sn, stitle, predefined_submenu);
			}
			predefined_submenu = new wxMenu;
			stitle = wxString(m_NamedFragments[i * 2], WX_DEFAULT_CONV);
			sn = i;
		} else {
			wxString mtitle(m_NamedFragments[i * 2], WX_DEFAULT_CONV);
			(predefined_submenu != NULL ? predefined_submenu : fmenu)->Append(myMenuID_PredefinedFragment + 1 + i, mtitle);
		}
	}
	if (predefined_submenu != NULL)
		fmenu->Append(myMenuID_PredefinedFragment + 1 + sn, stitle, predefined_submenu);
	Connect(myMenuID_PredefinedFragment + 1, myMenuID_PredefinedFragment + m_CountNamedFragments, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MyApp::OnFragmentMenuSelected), NULL, this);
}

void
MyApp::OnFragmentMenuSelected(wxCommandEvent& event)
{
	int index = event.GetId() - myMenuID_PredefinedFragment - 1;
	if (index < 0 || index >= m_CountNamedFragments)
		return;
	//  Open a predefined fragment as a new file
	char *errbuf;
	Molecule *mol = MoleculeNew();
	char *fullname;
	int result;

	asprintf(&fullname, "%s/Scripts/mbsf/%s", (const char *)(FindResourcePath().mb_str(wxConvFile)), m_NamedFragments[index * 2 + 1]);
	result = MoleculeLoadMbsfFile(mol, fullname, &errbuf);
	if (errbuf != NULL) {
		MyAppCallback_showScriptMessage("On loading %s:\n%s", m_NamedFragments[index * 2 + 1], errbuf);
		free(errbuf);
	}
	if (result != 0) {
		MyAppCallback_errorMessageBox("Cannot open named fragment %s\n(See console for detailed message)", m_NamedFragments[index * 2]);
		free(fullname);
		return;
	}
	free(fullname);
	MyDocument *doc = (MyDocument *)(DocManager()->CreateDocument(wxT(""), wxDOC_NEW));
	wxString title(m_NamedFragments[index * 2], WX_DEFAULT_CONV);
	title = _T("*") + title + _T("*");
	doc->SetMolecule(mol);
	if (mol->natoms > 1000)
		mol->mview->lineMode = 1;
	MainView_resizeToFit(mol->mview);
	
	//  Change the window title
	doc->SetTitle(title);
	//  Propagate the change of the window title
    wxList::compatibility_iterator node = doc->GetViews().GetFirst();
    while (node) {
        wxView *view = (wxView *)node->GetData();
        view->OnChangeFilename();
        node = node->GetNext();
    }
}

void
MyApp::OnUpdateUI(wxUpdateUIEvent& event)
{
	int uid = event.GetId();
	MainView *mview = MainViewCallback_activeView();
	if (uid >= myMenuID_CustomScript && uid < myMenuID_CustomScript + countScriptMenu) {
		//  Check the script menu
		//  If the frontmost window is RubyDialogFrame, then disable any script menu command
		wxWindow *w;
#if defined(__WXMAC__)
		void *MacGetActiveWindow(void);
		w = wxDynamicCast(wxNonOwnedWindow::GetFromWXWindow((WXWindow)MacGetActiveWindow()), wxWindow);
#else
		w = wxGetActiveWindow();
#endif
		if (wxDynamicCast(w, RubyDialogFrame) != NULL) {
			event.Enable(false);
			return;
		}
		Molecule *mol;
		int enabled, checked;
		char *title;
		int index = uid - myMenuID_CustomScript;
		if (mview == NULL)
			mol = NULL;
		else mol = mview->mol;
		checked = -1;
		title = NULL;
		enabled = Ruby_UpdateUI(index, mol, &checked, &title);
		if (checked >= 0)
			event.Check(checked != 0);
		if (title != NULL) {
			wxString wtext(title, WX_DEFAULT_CONV);
			event.SetText(wtext);
			free(title);
		}
		event.Enable(enabled != 0);
	} else if (uid >= myMenuID_PredefinedFragment && uid <= myMenuID_PredefinedFragment + m_CountNamedFragments) {
		event.Enable(true);
#if defined(__WXMAC__)
	} else if (uid >= wxID_OSX_MENU_FIRST && uid <= wxID_OSX_MENU_LAST) {
		event.Enable(true);
#endif
	} else {
		switch (uid) {
			case myMenuID_ExecuteScript:
			case myMenuID_OpenConsoleWindow:
			case myMenuID_EmptyConsoleWindow:
			case myMenuID_ViewParameterFilesList:
			case myMenuID_ViewGlobalParameters:
				event.Enable(true);
				return;
			default:
				if (mview == NULL)
					event.Enable(false);
				else
					event.Skip();
		}
	}
}

void
MyApp::OnExecuteScript(wxCommandEvent &event)
{
	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Choose Script File"), _T(""), _T(""), _T("All files (*.*)|*.*"), wxFD_OPEN | wxFD_CHANGE_DIR | wxFD_FILE_MUST_EXIST);
	if (dialog->ShowModal() == wxID_OK) {
		int status;
		RubyValue retval;
		wxString path = dialog->GetPath();
		
		//  Command line: execute_script('pathname')
		wxString cline(path);
		wxRegEx re(_T("[\\\\']"));   //  A backslash and a single-quote
		re.Replace(&cline, _T("\\\\\\0"));  //  A backslash followed by "\0"
		cline.Prepend(_T("execute_script('"));
		cline += _T("')");
		MyAppCallback_setConsoleColor(3);
		MyAppCallback_showScriptMessage("%s", (const char *)(cline.mb_str(wxConvFile)));
		MyAppCallback_showScriptMessage("\n");
		MyAppCallback_setConsoleColor(0);

		retval = MyAppCallback_executeScriptFromFile((const char *)(path.mb_str(wxConvFile)), &status);
		if (retval == (RubyValue)6 && status == -1)
			MyAppCallback_errorMessageBox("Cannot open Ruby script %s", (const char *)path.mb_str(wxConvFile));
		else if (status != 0)
			Ruby_showError(status);
	}
	dialog->Destroy();
}

void
MyApp::OnActivate(wxActivateEvent &event)
{
#if defined(__WXMAC__) || defined(__WXMSW__)
	MyFrame *frame = GetMainFrame();
	if (frame != NULL)
		frame->Show(false);  /*  Sometimes this "parent" frame gets visible and screw up the menus  */
#endif
	event.Skip();
}

void
MyApp::RequestOpenFilesByEvent(wxString& files)
{
	if (m_pendingFilesToOpen != NULL)
		m_pendingFilesToOpen->Append(files);
	else
		m_pendingFilesToOpen = new wxString(files);
    if (!m_pendingFilesToOpen->EndsWith(wxT("\n")))
        m_pendingFilesToOpen->Append(wxT("\n"));
	wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_openFilesByEvent);
	wxPostEvent(this, myEvent);
}

void
MyApp::OnOpenFilesByEvent(wxCommandEvent& event)
{
	if (m_pendingFilesToOpen == NULL)
		return;
    if (!gInitCompleted) {
        //  Repost this event and try again later
        wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_openFilesByEvent);
        wxPostEvent(this, myEvent);
        return;
    }
	OnOpenFiles(*m_pendingFilesToOpen);
	delete m_pendingFilesToOpen;
	m_pendingFilesToOpen = NULL;
}

bool
MyApp::OnOpenFiles(const wxString &files)
{
	Int start, end;
	bool success = true;
	int status;
	RubyValue retval;
	start = 0;
	while (1) {
		end = files.find(wxT("\n"), start);
		wxString file = files.Mid(start, (end == wxString::npos ? wxString::npos : end - start));
		if (file.Len() == 0)
			break;
		if (file.EndsWith(wxT(".rb")) || file.EndsWith(wxT(".mrb"))) {
			/*  Execute the file as a Ruby script  */
			retval = MyAppCallback_executeScriptFromFile((const char *)file.mb_str(wxConvFile), &status);
			if (status != 0) {
				if (retval == (RubyValue)6 && status == -1)
					MyAppCallback_errorMessageBox("Cannot open Ruby script: %s", (const char *)file.mb_str(wxConvFile));
				else
					Ruby_showError(status);
				return false;
			}
		} else {
			if (NULL == wxGetApp().DocManager()->CreateDocument(file, wxDOC_SILENT))
				success = false;
		}
		if (end == wxString::npos)
			break;
		start = end + 1;
	}
	return success;
}

wxString
MyApp::DefaultSettingsPath()
{
	wxString name = wxStandardPaths::Get().GetUserConfigDir();
	wxChar sep = wxFileName::GetPathSeparator();
	if (name[name.Len() - 1] != sep)
		name += sep;
	name += _T("Molby.settings");
	return name;
}

void
MyApp::LoadDefaultSettings()
{
	wxString name = DefaultSettingsPath();
	m_defaultSettings.clear();
	wxTextFile file(name);
	if (file.Exists() && file.Open()) {
		wxString line;
		int pos;
		for (line = file.GetFirstLine(); ; line = file.GetNextLine()) {
			if (line[0] == '#')
				continue;
			if ((pos = line.Find('=')) != wxNOT_FOUND) {
				wxString key = line.Left(pos);
				wxString value = line.Right(line.Length() - pos - 1);
				SetDefaultSetting(key, value);
			}
			if (file.Eof())
				break;
		}
		file.Close();
	}
}

void
MyApp::SaveDefaultSettings()
{
	wxString name = DefaultSettingsPath();
	wxTextFile file(name);
	if (!file.Exists())
		file.Create();
	else
		file.Open();
	file.Clear();
	MyStringHash::iterator it;
	for (it = m_defaultSettings.begin(); it != m_defaultSettings.end(); it++) {
		wxString key = it->first;
		wxString value = it->second;
		wxString line = key + _T("=") + value;
		file.AddLine(line);
	}
	file.Write();
	file.Close();
}

void
MyApp::SetDefaultSetting(const wxString& key, const wxString& value)
{
	//  TODO: The '=' and '#' characters may need to be escaped
	m_defaultSettings[key] = value;
}

wxString &
MyApp::GetDefaultSetting(const wxString& key)
{
	return m_defaultSettings[key];
}

MyListCtrl *
MyApp::GetGlobalParameterListCtrl()
{
	if (parameterFrame != NULL)
		return parameterFrame->GetListCtrl();
	else return NULL;
}

#define LOG_SUBPROCESS 0
#if LOG_SUBPROCESS
static FILE *fplog;
#endif

#if 0
void
MyApp::OnEndProcess(wxProcessEvent &event)
{
	m_processTerminated = true;
	m_processExitCode = event.GetExitCode();
#if LOG_SUBPROCESS
	if (fplog != NULL)
		fprintf(fplog, "OnEndProcess called\n");
#endif
//	delete m_process;
//	m_process = NULL;
}
#endif

int
MyApp::CallSubProcess(const char *cmdline, const char *procname, int (*callback)(void *), void *callback_data, FILE *fpout, FILE *fperr, int *exitstatus_p, int *pid_p)
{
	int status = 0;
	int callback_result = 0;
	int count = 0;
	bool progress_panel = false;
	char buf[256];
	size_t len, len_total;
	wxString cmdstr(cmdline, WX_DEFAULT_CONV);
	wxLongLong startTime;

	if (m_process != NULL)
		return -1;  //  Another process is already running (CallSubProcess() allows only one subprocess)
	
#if defined(__WXMSW__)
	extern int myKillAllChildren(long pid, wxSignal sig, wxKillError *krc);
#endif
	//  Show progress panel
	if (procname != NULL) {
		snprintf(buf, sizeof buf, "Running %s...", procname);
        ShowProgressPanel(buf);
		progress_panel = true;
	}
	startTime = wxGetUTCTimeMillis();
	
	//  Create log file in the document home directory
#if LOG_SUBPROCESS
    wxDateTime dateTime;
    dateTime.SetToCurrent();
	int nn = 0;
	{
		char *dochome = MyAppCallback_getDocumentHomeDir();
		snprintf(buf, sizeof buf, "%s/molby_subprocess.log", dochome);
		free(dochome);
		fplog = fopen(buf, "a");
		if (fplog == NULL)
			return -1;
	}
#endif

	//  Create proc object and call subprocess
	m_process = new wxBetterProcess(this, -1);
	m_process->Redirect();
	int flag = wxEXEC_ASYNC;
	flag |= wxEXEC_MAKE_GROUP_LEADER;
	long pid = ::wxExecute(cmdstr, flag, m_process);
	if (pid == 0) {
        if (progress_panel)
            HideProgressPanel();
		delete m_process;
#if LOG_SUBPROCESS
		fprintf(fplog, "Cannot start '%s'\n", cmdline);
		fclose(fplog);
#endif
		return -1;
	}
	if (pid_p != NULL)
		*pid_p = pid;
#if LOG_SUBPROCESS
	fprintf(fplog, "%s[DEBUG]pid = %ld\n", (const char *)(dateTime.FormatISOCombined(' ')), pid);
	fflush(fplog);
#endif
	
	//  Wait until process ends or user interrupts
    wxMemoryBuffer memBuffer;  //  Buffer to store standard output
    bool interrupted = false;
    wxString bufstr;
    wxString buferrstr;
    while (1) {
        m_process->GetLine(bufstr);
        if (bufstr.Length() > 0) {
#if LOG_SUBPROCESS
            dateTime.SetToCurrent();
            fprintf(fplog, "%s[STDOUT]%s", (const char *)(dateTime.FormatISOCombined(' ')), (const char *)bufstr);
            fflush(fplog);
#endif
            if (callback == DUMMY_CALLBACK) {
                const char *p = (const char *)bufstr;
                memBuffer.AppendData(p, strlen(bufstr));
            } else if (fpout != NULL && fpout != (FILE *)1) {
                fputs((const char *)bufstr, fpout);
            } else if (fpout == (FILE *)1) {
                MyAppCallback_setConsoleColor(0);
                MyAppCallback_showScriptMessage("%s", (const char *)bufstr);
            }
        }
        m_process->GetErrorLine(buferrstr);
        if (buferrstr.Length() > 0) {
#if LOG_SUBPROCESS
            dateTime.SetToCurrent();
            fprintf(fplog, "%s[STDERR]%s", (const char *)(dateTime.FormatISOCombined(' ')), buf);
            fflush(fplog);
#endif
            if (fperr != NULL && fperr != (FILE *)1) {
                fputs((const char *)buferrstr, fperr);
            } else if (fpout == (FILE *)1) {
                MyAppCallback_setConsoleColor(1);
                MyAppCallback_showScriptMessage("\n%s", (const char *)buferrstr);
                MyAppCallback_setConsoleColor(0);
            }
        }

        ::wxMilliSleep(25);
        if (callback != NULL && callback != DUMMY_CALLBACK) {
            callback_result = (*callback)(callback_data);
            if (callback_result != 0)
                interrupted = true;
        }
        if (progress_panel) {
            SetProgressValue(-1);
            if (IsInterrupted())
                interrupted = true;
        }

#if LOG_SUBPROCESS
        if (++nn >= 10) {
            dateTime.SetToCurrent();
            fprintf(fplog, "%s[DEBUG]pid %ld exists\n", (const char *)(dateTime.FormatISOCombined(' ')), pid);
            fflush(fplog);
            nn = 0;
        }
#endif
        if (m_process->IsTerminated() || !wxProcess::Exists(pid)) {
            /*  The subprocess has terminated  */
            status = m_process->GetStatus();
            break;
        } else if (interrupted) {
            /*  User interrupt  */
            int kflag = wxKILL_CHILDREN;
            wxKillError rc;
            if (
#if __WXMSW__
                myKillAllChildren(pid, wxSIGKILL, &rc) != 0
#else
                ::wxKill(pid, wxSIGTERM, &rc, kflag) != 0
#endif
                ) {
                switch (rc) {
                    case wxKILL_BAD_SIGNAL: status = -3; break; /* No such signal */
                    case wxKILL_ACCESS_DENIED: status = -4; break; /*  Permission denied  */
                    case wxKILL_NO_PROCESS: status = -5; break; /*  No such process  */
                    default: status = -6; break;  /*  unknown error  */
                }
            } else {
                if (callback_result != 0)
                    status = -3;  /*  Interrupt from callback  */
                else
                    status = -2;  /*  User interrupt  */
            }
            m_process->Detach();
            m_process = NULL;
            break;
        }
    }
    
    if (exitstatus_p != NULL)
        *exitstatus_p = status;

    if (progress_panel) {
        HideProgressPanel();
    }
    
	if (m_process != NULL) {
        m_process->Detach();
		m_process = NULL;
	}

	if (callback == DUMMY_CALLBACK) {
        char *membuf = NULL;
        size_t memsize = 0;
        memBuffer.AppendByte(0);
        memsize = memBuffer.GetDataLen();
        membuf = (char *)malloc(memsize);
        if (membuf != NULL) {
            memmove(membuf, memBuffer.GetData(), memsize);
#if __WXMSW__
            {
                /*  Convert "\r\n" to "\n"  */
                char *p, *pend;
                p = pend = membuf + strlen(membuf) + 1;
                while (--p >= membuf) {
                    if (*p == '\r') {
                        memmove(p, p + 1, pend - p);
                        pend--;
                    }
                }
            }
#endif
            *((char **)callback_data) = membuf;
        }
	}

	return status;
}

static int sTimerCount = 0;

void
MyApp::EnableTimerForDocument(MyDocument *doc)
{
	int i;
	if (doc == NULL)
		return;
	for (i = 0; i < m_CountTimerDocs; i++) {
		if (m_TimerDocs[i] == doc)
			return;
	}
	m_TimerDocs = (MyDocument **)realloc(m_TimerDocs, sizeof(MyDocument *) * (m_CountTimerDocs + 1));
	m_TimerDocs[m_CountTimerDocs++] = doc;
	if (m_Timer == NULL)
		m_Timer = new wxTimer(this, -1);
	if (!m_Timer->IsRunning())
		m_Timer->Start(100, wxTIMER_CONTINUOUS);
}

void
MyApp::DisableTimerForDocument(MyDocument *doc)
{
	int i;
	if (doc == NULL)
		return;
	for (i = 0; i < m_CountTimerDocs; i++) {
		if (m_TimerDocs[i] == doc) {
			//  Remove this document from the array
			if (i < m_CountTimerDocs - 1) {
				memmove(&m_TimerDocs[i], &m_TimerDocs[i + 1], sizeof(MyDocument *) * (m_CountTimerDocs - 1 - i));
			}
			m_CountTimerDocs--;
			if (m_CountTimerDocs == 0) {
				free(m_TimerDocs);
				m_TimerDocs = NULL;
				m_Timer->Stop();
			}
			break;
		}
	}
}

void
MyApp::TimerInvoked(wxTimerEvent &event)
{
	int i;
	sTimerCount++;
	for (i = 0; i < m_CountTimerDocs; i++) {
		m_TimerDocs[i]->TimerCallback(sTimerCount);
	}
}

void
MyApp::CheckIfAllWindowsAreGoneHandler(wxCommandEvent &event)
{
	int m = 0;
	wxWindowList::iterator iter;
	wxTopLevelWindow *win;
    for (iter = wxTopLevelWindows.begin(); iter != wxTopLevelWindows.end(); ++iter) {
        win = (wxTopLevelWindow *)(*iter);
		if (win != m_frameToBeDestroyed && win->IsShown())
			m++;
    }
	if (m == 0) {
		const char *p = MyAppCallback_getGlobalSettings(gSettingQuitOnCloseLastWindow);
		int quitFlag = 1;
		if (p != NULL && *p != 0)
			quitFlag = (atoi(p) != 0);
		if (quitFlag || MyAppCallback_messageBox("Do you want to quit Molby?", "Quit Molby", 3, 2)) {
			
			for (iter = wxTopLevelWindows.begin(); iter != wxTopLevelWindows.end(); ++iter) {
				win = (wxTopLevelWindow *)(*iter);
				if (win != m_frameToBeDestroyed)
					win->Destroy();  //  Avoid double destruction
			}
		} else {
			//  Show console window to avoid window-less state
			consoleFrame->Show();
		}
	}
}

void
MyApp::CheckIfAllWindowsAreGone(wxTopLevelWindow *frame)
{
	/*  On Windows, we should avoid the situation where all windows are hidden and
	    still the program is running. So that we check whether all windows are gone
	    and if so ask the user to quit the program. If user chooses not to quit, then
	    the console window is reopened and the program continues to run.  */
	m_frameToBeDestroyed = frame;
	wxCommandEvent myEvent(MyDocumentEvent, myMenuID_Internal_CheckIfAllWindowsAreGone);
	this->AddPendingEvent(myEvent);
}

void
MyApp::OnHelp(wxCommandEvent& WXUNUSED(event) )
{
    static wxString url;
    if (url.IsEmpty()) {
        url = FindResourcePath();
#if defined(__WXMSW__)
        if (url.SubString(0, 1) == wxT("\\\\")) {
            //  Network drive: convert to the drive letter
            wxBetterProcess *process = new wxBetterProcess(this, -1);
            process->Redirect();
            long pid = ::wxExecute("net use", wxEXEC_ASYNC | wxEXEC_MAKE_GROUP_LEADER, process);
            if (pid != 0) {
                wxRegEx re(wxT("^[[:space:]]+([A-Za-z]:)[[:space:]]+(.*)$"));
                wxString bufstr;
                while (process->GetLine(bufstr) >= 0) {
                    bufstr = bufstr.Trim();
                    if (!bufstr.IsEmpty() && re.Matches(bufstr)) {
                        wxString dr = re.GetMatch(bufstr, 1);
                        wxString path = re.GetMatch(bufstr, 2);
                        if (url.Left(path.Length()) == path) {
                            url = dr + url.Mid(path.Length());
                            break;
                        }
                    }
                }
                process->Detach();
            }
        }
#endif
        url.Replace(wxFILE_SEP_PATH, wxT("/"));
        url += "/MolbyDoc/ja/index.html";
    }
    wxLaunchDefaultBrowser(wxT("file:///") + url);
}

int
MyApp::FilterEvent(wxEvent &event)
{
#if 0
    static FILE *fp_eventlog = NULL;
    if (fp_eventlog == NULL) {
        char buf[32];
        int i = 0;
        while (1) {
            snprintf(buf, sizeof buf, "Molby_eventlog_%d.log", i);
            fp_eventlog = fopen(buf, "r");
            if (fp_eventlog == NULL) {
                fp_eventlog = fopen(buf, "wt");
                break;
            } else {
                fclose(fp_eventlog);
                i++;
            }
        }
    }
    if (fp_eventlog != NULL) {
        fprintf(fp_eventlog, "%d %d\n", event.GetEventType(), event.GetId());
        }
        fflush(fp_eventlog);
    }
#endif
    return -1;
}

#pragma mark ====== AboutBox ======

IMPLEMENT_CLASS(AboutDialog, wxDialog)
BEGIN_EVENT_TABLE(AboutDialog, wxDialog)
END_EVENT_TABLE()

AboutDialog::AboutDialog():
    wxDialog(NULL, -1, wxT("About Molby"))
{
    //  vsizer1 --> hsizer1 --> Molby icon
    //          |           |-> vsizer2 --> "Molby"
    //          |                       |-> version strings
    //          |
    //          |-> copyright messages

    char *s1, *s2, *s3;
    s1 = "Molby";
    Molby_getDescription(&s2, &s3);
    wxString str1(s1, WX_DEFAULT_CONV);
    wxString str2(s2, WX_DEFAULT_CONV);
    wxString str3(s3, WX_DEFAULT_CONV);
    free(s2);
    free(s3);
#if defined(__WXMSW__)
    wxFont *textFont0 = new wxFont(12, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD);
    wxFont *textFont1 = new wxFont(10, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
    wxFont *textFont2 = new wxFont(9, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
#else
    wxFont *textFont0 = new wxFont(14, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD);
    wxFont *textFont1 = new wxFont(12, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
    wxFont *textFont2 = new wxFont(11, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
#endif
    wxBoxSizer *vsizer1 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vsizer2 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *hsizer1 = new wxBoxSizer(wxHORIZONTAL);
    wxString tifname = wxGetApp().FindResourcePath() + wxFILE_SEP_PATH + wxT("bitmaps/molby_icon64.png");
    wxBitmap *molbyBitmap = new wxBitmap(tifname, wxBITMAP_TYPE_PNG);
    wxStaticText *stext1 = new wxStaticText(this, -1, wxT("Molby"));
    stext1->SetFont(*textFont0);
    wxStaticText *stext2 = new wxStaticText(this, -1, str2);
    stext2->SetFont(*textFont1);
    wxStaticBitmap *staticBitmap = new wxStaticBitmap(this, -1, *molbyBitmap);
    vsizer2->Add(stext1, 0, wxALL | wxEXPAND, 2);
    vsizer2->Add(stext2, 0, wxALL | wxEXPAND, 2);
    hsizer1->AddSpacer(20);
    hsizer1->Add(staticBitmap, 0, 0, 10);
    hsizer1->AddSpacer(20);
    hsizer1->Add(vsizer2, 0, wxALL | wxEXPAND, 5);
    wxStaticText *stext3 = new wxStaticText(this, -1, str3);
    stext3->SetFont(*textFont2);
    vsizer1->Add(hsizer1, 0, wxALL | wxEXPAND, 5);
    vsizer1->Add(stext3, 0, wxALL | wxEXPAND, 5);
    vsizer1->Add(this->CreateButtonSizer(wxOK), 0, wxALL | wxEXPAND, 10);
    vsizer1->Layout();
    this->SetSizerAndFit(vsizer1);
    this->Centre();
}

#pragma mark ====== MyFrame (top-level window) ======

/*
 * This is the top-level window of the application.
 */
 
IMPLEMENT_CLASS(MyFrame, wxDocParentFrame)
BEGIN_EVENT_TABLE(MyFrame, wxDocParentFrame)
    EVT_MENU(wxID_ABOUT, MyFrame::OnAbout)
END_EVENT_TABLE()

MyFrame::MyFrame(wxDocManager *manager, wxFrame *frame, const wxString& title,
    const wxPoint& pos, const wxSize& size, long type):
  wxDocParentFrame(manager, frame, wxID_ANY, title, pos, size, type, _T("myFrame"))
{
	editMenu = (wxMenu *) NULL;
#if defined(__WXMAC__)
	/*  Avoid this "dummy" top-level window to appear in the window menu.
	    It should not happen because MyApp::OnActivate() tries to hide this window,
	    but this is still here just in case.  */
//	OSStatus sts;
//	sts = ChangeWindowAttributes((WindowRef)m_macWindow, 0, kWindowInWindowMenuAttribute);
/*	printf("m_macWindow = %p, status = %d\n", m_macWindow, (int)sts); */
#endif
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event) )
{
    AboutDialog *d = new AboutDialog();
    if ( d->ShowModal() == wxID_OK )
    d->Destroy();
}

MyFrame *GetMainFrame(void)
{
	return frame;
}

#if 0
#pragma mark ====== Better wxProcess ======
#endif

void
wxBetterProcess::OnTerminate(int pid, int status)
{
    m_terminated = true;
    m_status = status;
}
wxKillError
wxBetterProcess::KillProcess(wxSignal sig, int flags)
{
    wxKillError retval = wxProcess::Kill(this->GetPid(), sig, flags);
    if (retval == wxKILL_OK)
        m_killSignal = sig;
    return retval;
}

int
wxBetterProcess::GetLineSub(wxString &outStr, wxInputStream *stream, wxMemoryBuffer &mbuf)
{
    int err = wxSTREAM_NO_ERROR;
    int trial = 0;
    char *p, *pp;
    long len;
    char buf[1024];
    if (stream == NULL)
        return -3;  //  No stderr stream
    while (1) {
        p = (char *)mbuf.GetData();
        len = mbuf.GetDataLen();
        if (len > 0) {
            pp = (char *)memchr(p, '\n', len);
            if (pp == NULL)
                pp = (char *)memchr(p, '\r', len);
            if (pp == NULL && stream->GetLastError() == wxSTREAM_EOF) {
                //  If EOF, then return all remaining data (without '\n')
                pp = p + mbuf.GetDataLen() - 1;  //  Point to the last char
            }
            if (pp != NULL) {
                //  Return one line and string length
                outStr = wxString(p, wxConvUTF8, pp - p + 1);
                memmove(p, pp + 1, len - (pp - p + 1));
                m_stdout.SetDataLen(len - (pp - p + 1));
                return pp - p + 1;
            }
        }
        if (trial > 0) {
            //  stream->Read() is called only once
            outStr = _T("");
            if (err == wxSTREAM_EOF)
                return -1;  //  EOF and no data left
            return 0;  //  Not EOF, but no data is available at present
        }
        len = 0;
        if (stream->CanRead()) {
            //  We need to read by one character because wxInputStream has
            //  no way to give the available number of bytes
            stream->Read(buf, sizeof buf);
            err = stream->GetLastError();
            if (err != wxSTREAM_NO_ERROR && err != wxSTREAM_EOF)
                return -2;  //  Some read error
            len = stream->LastRead();
        } else err = stream->GetLastError();
        if (len > 0)
            mbuf.AppendData(buf, len);
        trial++;
    }
}

int
wxBetterProcess::GetLine(wxString &outStr)
{
    return GetLineSub(outStr, this->GetInputStream(), m_stdout);
}

int
wxBetterProcess::GetErrorLine(wxString &outStr)
{
    return GetLineSub(outStr, this->GetErrorStream(), m_stderr);
}

int
wxBetterProcess::PutLine(wxString str)
{
    wxOutputStream *stream = this->GetOutputStream();
    if (stream == NULL)
        return -3;  //  No stdin stream
    const char *p = str.utf8_str();
    long len = strlen(p);
    if (len > 0)
        m_stdin.AppendData(p, len);
    char *pp = (char *)m_stdin.GetData();
    len = m_stdin.GetDataLen();
    if (len == 0)
        return 0;
    stream->Write(pp, len);
    long len2 = stream->LastWrite();
    if (len2 > 0) {
        memmove(pp, pp + len2, len - len2);
        m_stdin.SetDataLen(len - len2);
    }
    return len2;
}

void
wxBetterProcess::CloseOutput()
{
    //  We must flush the data in the internal buffer before closing the output
    while (PutLine("") > 0) {}
    wxProcess::CloseOutput();  //  Call the original version
}

#pragma mark ====== Plain-C interface ======

char *
MyAppCallback_getGUIDescriptionString(void)
{
	static char *desc = NULL;
	if (desc == NULL) {
		asprintf(&desc,
			"AmberTools 1.3, http://ambermd.org/\n"
			"  Copyright (C) Junmei Wang, Ross C. Walker, "
			  "Michael F. Crowley, Scott Brozell and David A. Case\n"
			"ORTEP-III, http://web.ornl.gov/sci/ortep/\n"
			"  Michael N. Burnett and Carroll K. Johnson, "
			  "Oak Ridge National Laboratory Report ORNL-6895, "
			  "1996.\n"
			"wxWidgets %d.%d.%d, http://www.wxwidgets.org/\n"
		    "  Copyright (C) 1992-2013 Julian Smart, Vadim "
			  "Zeitlin, Stefan Csomor, Robert Roebling,\n"
            "  and other members of the wxWidgets team\n"
            "  Portions (C) 1996 Artificial Intelligence "
			  "Applications Institute\n",
			wxMAJOR_VERSION, wxMINOR_VERSION, wxRELEASE_NUMBER);
	}
	return desc;
}

void
MyAppCallback_loadGlobalSettings(void)
{
    if (!gUseGUI)
        return;
	wxGetApp().LoadDefaultSettings();
}

void
MyAppCallback_saveGlobalSettings(void)
{
    if (!gUseGUI)
        return;
	wxGetApp().SaveDefaultSettings();
}

/*  Note on the global settings  */
/*  Global settings are stored in a file in the form key="value", where
    the "value" is the 'inspect'-ed representation of Ruby values.
    So that, if the value is a string like 'aaaa', the stored value is "aaaa" (with quotes),
    not aaaa (without quotes). This is convenient for access from Ruby scripts, but it needs
    care for access from C. For C-level access, use MyAppCallback_getGlobalSettingsWithType() and
    MyAppCallback_setGlobalSettingsWithType().  */
char *
MyAppCallback_getGlobalSettings(const char *key)
{
    if (!gUseGUI)
        return NULL;
    wxString wxkey(key, WX_DEFAULT_CONV);
    wxString wxvalue = wxGetApp().GetDefaultSetting(wxkey);
    return strdup(wxvalue.mb_str(WX_DEFAULT_CONV));
}

void
MyAppCallback_setGlobalSettings(const char *key, const char *value)
{
    if (!gUseGUI)
        return;
	wxString wxkey(key, WX_DEFAULT_CONV);
	wxString wxvalue(value, WX_DEFAULT_CONV);
	wxGetApp().SetDefaultSetting(wxkey, wxvalue);
}

int
MyAppCallback_getGlobalSettingsWithType(const char *key, int type, void *ptr)
{
    int retval, temp;
	char *s = MyAppCallback_getGlobalSettings(key);
	char desc[] = SCRIPT_ACTION("s; ");
	desc[sizeof(desc) - 2] = type;
    temp = gMolActionNoErrorDialog;
    gMolActionNoErrorDialog = 1;
	retval = MolActionCreateAndPerform(NULL, desc, "eval", s, ptr);
	free(s);
    gMolActionNoErrorDialog = temp;
	return retval;
}

int
MyAppCallback_setGlobalSettingsWithType(const char *key, int type, const void *ptr)
{
	const char *cmd = "set_global_settings";
	switch (type) {
		case 'i': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("si"), cmd, key, *((const Int *)ptr));
		case 'd': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("sd"), cmd, key, *((const Double *)ptr));
		case 's': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("ss"), cmd, key, (const char *)ptr);
		case 'v': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("sv"), cmd, key, (const Vector *)ptr);
		case 't': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("st"), cmd, key, (const Transform *)ptr);
		default:
			MyAppCallback_errorMessageBox("Internal error: unsupported format '%c' at line %d, file %s", type, __LINE__, __FILE__);
			return -2;
	}
}

int
MyAppCallback_checkInterrupt(void)
{
    if (!gUseGUI)
        return 0;
    return wxGetApp().IsInterrupted();
}

void
MyAppCallback_showProgressPanel(const char *msg)
{
    if (!gUseGUI)
        return;
    wxGetApp().ShowProgressPanel(msg);
}

void
MyAppCallback_hideProgressPanel(void)
{
    if (!gUseGUI)
        return;
    wxGetApp().HideProgressPanel();
}

void
MyAppCallback_setProgressValue(double dval)
{
    if (!gUseGUI)
        return;
    wxGetApp().SetProgressValue(dval);
}

void
MyAppCallback_setProgressMessage(const char *msg)
{
    if (!gUseGUI)
        return;
    wxGetApp().SetProgressMessage(msg);
}

int
MyAppCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize)
{
    if (!gUseGUI) {
        buf[0] = 0;
        return 0;
    }
	wxDialog *dialog = new wxDialog(NULL, -1, _T("Input request"), wxDefaultPosition);
	wxStaticText *stext;
	wxTextCtrl *tctrl;
	int retval;
	wxString pstr(prompt, WX_DEFAULT_CONV);
	wxString defstr(buf, WX_DEFAULT_CONV);
	{	//  Vertical sizer containing [prompt, textbox, buttons]
		wxBoxSizer *sizer1;
		sizer1 = new wxBoxSizer(wxVERTICAL);
		stext = new wxStaticText(dialog, -1, pstr, wxDefaultPosition, wxSize(200, 22));
		sizer1->Add(stext, 0, wxEXPAND | wxALL, 6);
		tctrl = new wxTextCtrl(dialog, -1, defstr, wxDefaultPosition, wxSize(200, 22));
		sizer1->Add(tctrl, 0, wxEXPAND | wxALL, 6);
		wxSizer *bsizer = dialog->CreateButtonSizer(wxOK | wxCANCEL);
		sizer1->Add(bsizer, 0, wxEXPAND | wxALL, 6);
		sizer1->Layout();
		dialog->SetSizerAndFit(sizer1);
		dialog->Centre(wxBOTH);
        tctrl->SelectAll();
		tctrl->SetFocus();
	}
	if (dialog->ShowModal() == wxID_OK) {
		strncpy(buf, (const char *)(tctrl->GetValue().mb_str(WX_DEFAULT_CONV)), bufsize - 1);
		buf[bufsize - 1] = 0;
		retval = 1;
	} else {
		retval = 0;
	}
	dialog->Destroy();
	return retval;
}

/*  Generic message box.  Flags is a bitwise OR of 1 (OK) and 2 (Cancel). Icon is either
    1 (information), 2 (exclamation), or 3 (stop).  */
int
MyAppCallback_messageBox(const char *message, const char *title, int flags, int icon)
{
    if (!gUseGUI) {
        printf("%s\n%s\n", title, message);
        return 1;
    }
    
	int wxflags, wxicon, retval;
	if (!wxGetApp().IsMainLoopRunning()) {
		MyAppCallback_setConsoleColor(1);
		MyAppCallback_showScriptMessage("*** %s ***\n%s\n", message, title);
		MyAppCallback_setConsoleColor(0);
		return 1;
	}
	if (flags == 0)
		flags = 1;
	wxflags = ((flags & 1) ? wxOK : 0) | ((flags & 2) ? wxCANCEL : 0);
	switch (icon) {
		case 3: wxicon = wxICON_ERROR; break;
		case 2: wxicon = wxICON_EXCLAMATION; break;
		default: wxicon = wxICON_INFORMATION; break;
	}
	wxString wxmessage(message, WX_DEFAULT_CONV);
	wxString wxtitle(title, WX_DEFAULT_CONV);
	retval = ::wxMessageBox(wxmessage, wxtitle, wxflags | wxicon);
	return (retval == wxOK ? 1 : 0);
}

void
MyAppCallback_errorMessageBox(const char *fmt, ...)
{
    if (!gUseGUI) {
        va_list ap;
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        return;
    }
	char *s;
	int need_free = 0;
	va_list ap;
	va_start(ap, fmt);
	if (strchr(fmt, '%') == 0) {
		s = (char *)fmt;
	} else if (strcmp(fmt, "%s") == 0) {
		s = va_arg(ap, char *);
	} else {
		vasprintf(&s, fmt, ap);
		need_free = 1;
	}
	MyAppCallback_messageBox(s, "Error", 0, 3);
	if (need_free)
		free(s);
}
	
char *
MyAppCallback_getHomeDir(void)
{
	char *s;
#if __WXMSW__
	/*  wxFileName::GetHomeDir() may return unexpected value under MSYS  */
	s = getenv("USERPROFILE");
#else
	s = getenv("HOME");
#endif
	return (s == NULL ? NULL : strdup(s));
}

#if __WXMSW__
#include <Shlobj.h>
#endif

char *
MyAppCallback_getDocumentHomeDir(void)
{
#if __WXMSW__
	char appData[MAX_PATH * 2];
	HRESULT hResult;
	hResult = SHGetFolderPathA(NULL, CSIDL_PERSONAL, NULL, 0, appData);
	if (hResult == S_OK) {
		return strdup(appData);
	} else {
		return MyAppCallback_getHomeDir();
	}
#else
	char *s;
	s = getenv("HOME");
	return (s == NULL ? NULL : strdup(s));
#endif
}

int
MyAppCallback_registerScriptMenu(const char *title)
{
    if (!gUseGUI)
        return -1;
	return wxGetApp().RegisterScriptMenu(title);
}

int
MyAppCallback_lookupScriptMenu(const char *title)
{
    if (!gUseGUI)
        return 0;
	return wxGetApp().LookupScriptMenu(title);
}

RubyValue
MyAppCallback_executeScriptFromFile(const char *cpath, int *status)
{
    if (!gUseGUI) {
        return 0;
    }
    
	RubyValue retval;
	wxString cwd = wxFileName::GetCwd();
	wxString path(cpath, wxConvFile);
	char *p = strdup(cpath);
	char sep = wxFileName::GetPathSeparator();
	char *pp, *script = NULL;
	if ((pp = strrchr(p, sep)) != NULL) {
		*pp++ = 0;
		wxString dirname(p, wxConvFile);
		wxFileName::SetCwd(dirname);
	} else pp = p;
	
	/*  Read the content of the file  */
	FILE *fp = fopen(cpath, "rb");
	if (fp != NULL) {
		off_t len;
		fseek(fp, 0, SEEK_END);
		len = ftell(fp);
		fseek(fp, 0, SEEK_SET);
		script = (char *)malloc(len + 1);
		if (script!= NULL) {
			fread(script, 1, len, fp);
			script[len] = 0;
		}
		fclose(fp);
	}

	if (script == NULL) {
		*status = -1;
		return (RubyValue)6;  /*  Cannot open file  */
	}
	
	/*  Check the encoding specification, and if present convert it to default encoding  */
	if (0) {
		char *lp = script, *eolp;
		int n = 0;
		while (n < 2) {  /*  Check the first two lines  */
			while (*lp && isspace(*lp))
				lp++;
			if (*lp == 0 || *lp++ != '#')  /*  Check only the comment line  */
				break;
			if (*lp == 0)
				break;
			if (*lp == '!') { /*  Shebang line  */
				while (*lp && *lp != '\n')
					lp++;  /*  Skip until end of line  */
				n++;
				lp++;
				continue;
			}
			for (eolp = lp; *eolp && *eolp != '\n'; eolp++);
			if (*eolp != '\n')
				break;
			*eolp = 0;  /*  Limit the search area  */
			lp = strstr(lp, "coding:");
			*eolp = '\n';  /*  Restore original string  */
			if (lp != NULL) {
				lp += 7;
				while (*lp && isspace(*lp))
					lp++;
				if (strncasecmp(lp, "shift-jis", 9) == 0) {
					wxString s(script, wxCSConv(wxT("cp932")));
					free(script);
					script = strdup(s.mb_str(WX_DEFAULT_CONV));
				} else if (strncasecmp(lp, "utf-8", 5) == 0) {
					wxString s(script, wxConvUTF8);
					free(script);
					script = strdup(s.mb_str(WX_DEFAULT_CONV));
				}
				break;
			}
			lp = eolp + 1;
			n++;
		}
	}
	
	retval = Molby_evalRubyScriptOnMolecule(script, MoleculeCallback_currentMolecule(), pp, status);
	
	free(script);
	free(p);
	wxFileName::SetCwd(cwd);
	return retval;
}

void MyAppCallback_beginUndoGrouping(void)
{
    if (!gUseGUI)
        return;
	wxList &doclist = wxGetApp().DocManager()->GetDocuments();
	wxList::iterator iter;
	for (iter = doclist.begin(); iter != doclist.end(); ++iter) {
		((MyDocument *)(*iter))->BeginUndoGrouping();
	}
}

void MyAppCallback_endUndoGrouping(void)
{
    if (!gUseGUI)
        return;
	wxList &doclist = wxGetApp().DocManager()->GetDocuments();
	wxList::iterator iter;
	for (iter = doclist.begin(); iter != doclist.end(); ++iter) {
		((MyDocument *)(*iter))->EndUndoGrouping();
	}
}

int MyAppCallback_callSubProcess(const char *cmdline, const char *procname, int (*callback)(void *), void *callback_data, FILE *output, FILE *errout, int *exitstatus_p, int *pid_p)
{
    if (!gUseGUI)
        return system(cmdline);
	return wxGetApp().CallSubProcess(cmdline, procname, callback, callback_data, output, errout, exitstatus_p, pid_p);
}

void MyAppCallback_showConsoleWindow(void)
{
    if (!gUseGUI)
        return;
	ConsoleFrame *frame = wxGetApp().GetConsoleFrame();
	frame->Show(true);
	frame->Raise();
}

void MyAppCallback_hideConsoleWindow(void)
{
    if (!gUseGUI)
        return;
	ConsoleFrame *frame = wxGetApp().GetConsoleFrame();
	frame->Hide();
}

void MyAppCallback_bell(void)
{
    if (!gUseGUI)
        return;
    wxBell();
}

int MyAppCallback_playSound(const char *filename, int flag)
{
    if (!gUseGUI)
        return 0;
	unsigned uflag = wxSOUND_SYNC;
	if (flag == 1)
		uflag = wxSOUND_ASYNC;
	else if (flag == 3)
		uflag = wxSOUND_ASYNC | wxSOUND_LOOP;
	wxString fnamestr(filename, wxConvFile);
	bool retval = wxSound::Play(fnamestr, uflag);
	return retval;
}

void MyAppCallback_stopSound(void)
{
    if (!gUseGUI)
        return;
	wxSound::Stop();
}

void MyAppCallback_initImageHandlers(void)
{
	static bool handlers_init = false;
	if (!handlers_init) {
		wxInitAllImageHandlers();
		handlers_init = true;
	}
}
