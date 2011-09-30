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

#include "MyApp.h"
#include "MyDocument.h"
#include "MoleculeView.h"
#include "ConsoleFrame.h"
#include "ProgressFrame.h"
#include "GlobalParameterFrame.h"
#include "GlobalParameterFilesFrame.h"
#include "MyMBConv.h"

#include "../MolLib/MolLib.h"
#include "../MolLib/Ruby_bind/Molby.h"
#include "../MolLib/Missing.h"

#include <wchar.h>
#include <stdio.h>

#if defined(__WXMAC__)
#include <CoreFoundation/CoreFoundation.h>
#undef T_DATA
#include <Carbon/Carbon.h>
#include <sys/wait.h>  /* for waitpid()  */
#endif

#pragma mark ====== MyApp ======

static char *sLastBuildString = "";

MyFrame *frame = (MyFrame *) NULL;

IMPLEMENT_APP(MyApp)

//IMPLEMENT_CLASS(MyApp, wxApp)

BEGIN_EVENT_TABLE(MyApp, wxApp)
	//EVT_KEY_DOWN(MyApp::OnChar)
	//EVT_MOUSE_EVENTS(MyApp::OnMouseEvent)
	EVT_COMMAND(MyDocumentEvent_scriptMenuModified, MyDocumentEvent, MyApp::OnScriptMenuModified)
	EVT_UPDATE_UI_RANGE(myMenuID_MyFirstMenuItem, myMenuID_MyLastMenuItem, MyApp::OnUpdateUI)
	EVT_MENU(myMenuID_ExecuteScript, MyApp::OnExecuteScript)
	EVT_MENU(myMenuID_OpenConsoleWindow, MyApp::OnOpenConsoleWindow)
//	EVT_MENU(myMenuID_ReadParameters, MyApp::OnReadParameters)
	EVT_MENU(myMenuID_ViewGlobalParameters, MyApp::OnViewGlobalParameters)
	EVT_MENU(myMenuID_ViewParameterFilesList, MyApp::OnViewParameterFilesList)
	EVT_MENU(myMenuID_ImportAmberLib, MyApp::OnImportAmberLib)
#if defined(__WXMAC__)
	EVT_ACTIVATE(MyApp::OnActivate)
#endif
	EVT_END_PROCESS(-1, MyApp::OnEndProcess)
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
	m_progressFrame = NULL;
	m_processTerminated = false;
	m_processExitCode = 0;
	countScriptMenu = 0;
	scriptMenuCommands = NULL;
	scriptMenuTitles = NULL;
	scriptMenuModifiedEventPosted = false;
	parameterFrame = NULL;
	parameterFilesFrame = NULL;
	consoleFrame = NULL;
}

bool MyApp::OnInit(void)
{
	//  Set defaults
#ifdef __WXMAC__
	wxSystemOptions::SetOption(wxT("mac.listctrl.always_use_generic"), 1);
#endif

#if __WXMSW__
	{
		//  Check if the same application is already running
		char *buf;
		asprintf(&buf, "Molby-%s", (const char *)wxGetUserId().mb_str(WX_DEFAULT_CONV));
		wxString name(buf, WX_DEFAULT_CONV);
		malloc(16);
		free(buf);
		m_checker = new wxSingleInstanceChecker(name);
		if (m_checker->IsAnotherRunning()) {
			wxLogError(_T("Molby is already running."));
			return false;
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
	new wxDocTemplate(m_docManager, _T("Gaussian Output File"), _T("*.out"), _T(""), _T("out"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Gaussian Checkpoint File"), _T("*.fchk;*.fch"), _T(""), _T("fchk"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("GAMESS Input File"), _T("*.inp"), _T(""), _T("inp"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("GAMESS Output File"), _T("*.log"), _T(""), _T("log"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("GAMESS DAT File"), _T("*.dat"), _T(""), _T("dat"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("ORTEP Input File"), _T("*.tep"), _T(""), _T("tep"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("SHELX Input File"), _T("*.ins;*.res"), _T(""), _T("ins"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Crystallographic Information File"), _T("*.cif"), _T(""), _T("cif"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Cartesian"), _T("*.xyz"), _T(""), _T("xyz"), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));
	new wxDocTemplate(m_docManager, _T("Any Molecule"), _T("*.*"), _T(""), _T(""), _T("Molecule Doc"), _T("Molecule View"), CLASSINFO(MyDocument), CLASSINFO(MoleculeView));

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

#if defined(__WXMAC__)
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
		static const char fname[] = "startup.rb";
		wxString dirname = FindResourcePath();
		char *wbuf;
	
		dirname += wxFILE_SEP_PATH;
		dirname += wxT("Scripts");
		wxString cwd = wxGetCwd();
		wxSetWorkingDirectory(dirname);

		/*  Read build information (for About dialog)  */
		{
			FILE *fp = fopen("../buildInfo.txt", "r");
			if (fp != NULL) {
				char buf[200];
				if (fgets(buf, sizeof(buf), fp) != NULL) {
					char *p1 = strchr(buf, '\"');
					char *p2 = strrchr(buf, '\"');
					if (p1 != NULL && p2 != NULL && p2 - p1 > 1) {
						memmove(buf, p1 + 1, p2 - p1 - 1);
						buf[p2 - p1 - 1] = 0;
						asprintf(&sLastBuildString, "Last compile: %s\n", buf);
					}
				}
				fclose(fp);
			}
		}
		
		/*  Read atom display parameters  */
		if (ElementParameterInitialize("element.par", &wbuf) != 0) {
			SetConsoleColor(1);
			AppendConsoleMessage(wbuf);
			SetConsoleColor(0);
			free(wbuf);
		}
		
		/*  Read default parameters  */
		ParameterReadFromFile(gBuiltinParameters, "default.par", &wbuf, NULL);
		if (wbuf != NULL) {
			SetConsoleColor(1);
			AppendConsoleMessage(wbuf);
			SetConsoleColor(0);
			free(wbuf);
		}

		wxString fnamestr(fname, wxConvFile);
		Molby_startup(wxFileExists(fnamestr) ? fname : NULL, (const char *)dirname.mb_str(wxConvFile));
		
		wxSetWorkingDirectory(cwd);
		MyAppCallback_showScriptMessage("%% ");
	}
	
	/*  Open given files as MyDocument  */
	if (argc == 1) {
#if __WXMSW__
		m_docManager->CreateDocument(wxEmptyString, wxDOC_NEW);
#endif
	} else {
		while (argc > 1) {
			wxString file(argv[1]);
			m_docManager->CreateDocument(file, wxDOC_SILENT);
			argc--;
			argv++;
		}
	}
	
	return true;
}

//  Create Menu Bars
//  kind == 0: main menu
//  kind == 1: molecule window
//  kind == 2: console window
wxMenuBar *
MyApp::CreateMenuBar(int kind, wxMenu **out_file_history_menu, wxMenu **out_edit_menu)
{
	
	//// Make a menubar
	wxMenu *file_menu = new wxMenu;

	file_menu->Append(wxID_NEW, _T("&New...\tCtrl-N"));
	file_menu->Append(wxID_OPEN, _T("&Open...\tCtrl-O"));
	if (out_file_history_menu != NULL) {
		*out_file_history_menu = new wxMenu;
		file_menu->AppendSubMenu(*out_file_history_menu, _T("Open Recent"));
		m_docManager->FileHistoryAddFilesToMenu(*out_file_history_menu);
		m_docManager->FileHistoryUseMenu(*out_file_history_menu);  //  Should be removed when menu is discarded
	}

	file_menu->AppendSeparator();
	file_menu->Append(wxID_CLOSE, _T("&Close\tCtrl-W"));
	file_menu->Append(wxID_SAVE, _T("&Save\tCtrl-S"));
	file_menu->Append(wxID_SAVEAS, _T("Save &As..."));	
	
	file_menu->AppendSeparator();
	file_menu->Append(myMenuID_Import, _T("Import..."));	
	file_menu->Append(myMenuID_Export, _T("Export..."));	
	
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
	
	wxMenu *show_menu = new wxMenu;
	show_menu->Append(myMenuID_FitToScreen, _T("Fit To Screen\tCtrl-T"));
	show_menu->Append(myMenuID_CenterSelection, _T("Center Selection"));
	show_menu->AppendSeparator();
	show_menu->Append(myMenuID_ShowUnitCell, _T("Show Unit Cell"), _T(""), wxITEM_CHECK);
/*	show_menu->Append(myMenuID_ShowPeriodicBox, _T("Show Periodic Box"), _T(""), wxITEM_CHECK); */
	show_menu->Append(myMenuID_ShowHydrogens, _T("Show Hydrogen Atoms"), _T(""), wxITEM_CHECK);
	show_menu->Append(myMenuID_ShowDummyAtoms, _T("Show Dummy Atoms"), _T(""), wxITEM_CHECK);
	show_menu->Append(myMenuID_ShowExpandedAtoms, _T("Show Expanded Atoms"), _T(""), wxITEM_CHECK);
	show_menu->Append(myMenuID_ShowEllipsoids, _T("Show Ellipsoids"), _T(""), wxITEM_CHECK);
	show_menu->Append(myMenuID_ShowRotationCenter, _T("Show Rotation Center"), _T(""), wxITEM_CHECK);
	show_menu->AppendSeparator();
	show_menu->Append(myMenuID_HideSelected, _T("Hide Selected"), _T(""));
	show_menu->Append(myMenuID_HideUnselected, _T("Hide Unselected"), _T(""));
	show_menu->Append(myMenuID_HideReverse, _T("Hide Reverse"), _T(""));
	show_menu->Append(myMenuID_ShowAllAtoms, _T("Show All Atoms"), _T(""));
	show_menu->AppendSeparator();
	show_menu->Append(myMenuID_ShowGraphite, _T("Show Graphite..."));
	show_menu->AppendSeparator();
	show_menu->Append(myMenuID_LineMode, _T("Line Mode"), _T(""), wxITEM_CHECK);

	wxMenu *md_menu = new wxMenu;
	md_menu->Append(myMenuID_MolecularDynamics, _T("Molecular Dynamics..."));
	md_menu->Append(myMenuID_Minimize, _T("Minimize..."));
	md_menu->Append(myMenuID_StopMDRun, _T("Stop\tCtrl-."));
	md_menu->AppendSeparator();
//	md_menu->Append(myMenuID_ReadParameters, _T("Read Parameters..."));	
	md_menu->Append(myMenuID_ViewGlobalParameters, _T("View Global Parameters..."));
	md_menu->Append(myMenuID_ViewParameterFilesList, _T("Load/Unload Global Parameters..."));
	md_menu->AppendSeparator();
	md_menu->Append(myMenuID_DefinePeriodicBox, _T("Define Unit Cell..."));
	md_menu->Append(myMenuID_ShowPeriodicImage, _T("Show Periodic Image..."));
	md_menu->Append(myMenuID_PressureControl, _T("Pressure Control..."));
/*	md_menu->Append(myMenuID_DefineSymmetry, _T("Define Symmetry Operations..."));
	md_menu->Append(myMenuID_ExpandBySymmetry, _T("Expand by Symmetry...")); */
	md_menu->AppendSeparator();
	wxMenu *md_tools_menu = new wxMenu;
	md_tools_menu->Append(myMenuID_RunAntechamber, _T("Antechamber/parmchk..."));
	md_tools_menu->Append(myMenuID_RunResp, _T("GAMESS/RESP..."));
	md_tools_menu->Append(myMenuID_CreateSanderInput, _T("Create SANDER input..."));
	md_tools_menu->Append(myMenuID_ImportAmberLib, _T("Import AMBER Lib..."));
	md_tools_menu->Append(myMenuID_ImportAmberFrcmod, _T("Import AMBER Frcmod..."));
	md_menu->Append(myMenuID_MDTools, _T("Tools"), md_tools_menu);

	wxMenu *qc_menu = new wxMenu;
	qc_menu->Append(myMenuID_CreateGamessInput, _T("Create GAMESS input..."));
	qc_menu->Append(myMenuID_CreateMOCube, _T("Create MO cube..."));
	
	wxMenu *script_menu = new wxMenu;
	script_menu->Append(myMenuID_ExecuteScript, _T("Execute Script..."));
	script_menu->Append(myMenuID_OpenConsoleWindow, _T("Open Console Window..."));
	script_menu->AppendSeparator();
	countNonCustomScriptMenu = script_menu->GetMenuItemCount();

	wxMenu *help_menu = new wxMenu;
	help_menu->Append(wxID_ABOUT, _T("&About...\tF1"));
	
	wxMenuBar *menu_bar = new wxMenuBar;
	
	menu_bar->Append(file_menu, _T("&File"));
	menu_bar->Append(edit_menu, _T("&Edit"));
	menu_bar->Append(show_menu, _T("Show"));
	menu_bar->Append(md_menu, _T("MM/MD"));
	menu_bar->Append(qc_menu, _T("QChem"));
	menu_bar->Append(script_menu, _T("&Script"));
	menu_bar->Append(help_menu, _T("&Help"));
	
	UpdateScriptMenu(menu_bar);
	
	return menu_bar;
}

#if __WXMAC__
/*  When the application is launched without any documents, an empty document is opened.
    This should be implemented by overriding this special method; parsing argc/argv does
    not work, because the list of files is passed through an Apple Event.  */
void
MyApp::MacNewFile()
{
	m_docManager->CreateDocument(_T(""), wxDOC_NEW);
}
#endif

int MyApp::OnExit(void)
{
	SaveDefaultSettings();
    delete m_docManager;
#if __WXMSW__
	delete m_checker;
#endif
    return 0;
}

int
MyApp::AppendConsoleMessage(const char *mes)
{
	wxTextCtrl *textCtrl;
	if (consoleFrame != NULL && (textCtrl = consoleFrame->textCtrl) != NULL) {
		wxString string(mes, WX_DEFAULT_CONV);
		textCtrl->AppendText(string);
		return string.Len();
	} else return 0;
}

void
MyApp::FlushConsoleMessage()
{
	wxTextCtrl *textCtrl = consoleFrame->textCtrl;
	textCtrl->Refresh();
	textCtrl->Update();
}

void
MyApp::SetConsoleColor(int color)
{
	wxTextCtrl *textCtrl = consoleFrame->textCtrl;
	static wxTextAttr *col[4];
	if (col[0] == NULL) {
		col[0] = new wxTextAttr(*wxBLACK);
		col[1] = new wxTextAttr(*wxRED);
		col[2] = new wxTextAttr(*wxGREEN);
		col[3] = new wxTextAttr(*wxBLUE);
	}
	textCtrl->SetDefaultStyle(*(col[color % 4]));
}

void
MyApp::ShowProgressPanel(const char *mes)
{
	wxString string((mes ? mes : ""), WX_DEFAULT_CONV);
	if (m_progressFrame == NULL) {
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
	}
}

void
MyApp::HideProgressPanel()
{
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
}

void
MyApp::SetProgressValue(double dval)
{
	if (m_progressFrame != NULL) {
		m_progressFrame->SetProgressValue(dval);
	}
}

void
MyApp::SetProgressMessage(const char *mes)
{
	if (m_progressFrame != NULL) {
		wxString string((mes ? mes : ""), WX_DEFAULT_CONV);
		m_progressFrame->SetProgressMessage(string);
	}
}

int
MyApp::IsInterrupted()
{
	if (m_progressFrame != NULL)
		return m_progressFrame->CheckInterrupt();
	else {
		if (::wxGetKeyState(WXK_ESCAPE))
			return 1;
		else return 0;		
	}
}

/*
#warning "TODO: Move this to MyDocument and 'import parameters' "

 void
MyApp::OnReadParameters(wxCommandEvent& event)
{
	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Choose Parameter File"), _T(""), _T(""), _T("All files (*.*)|*.*"), wxFD_OPEN | wxFD_CHANGE_DIR | wxFD_FILE_MUST_EXIST);
	if (dialog->ShowModal() == wxID_OK) {
		char *p = strdup((const char *)(dialog->GetPath().mb_str(wxConvFile)));
		char *wbuf;
		ParameterReadFromFile(NULL, p, &wbuf, NULL);
		if (wbuf != NULL) {
			SetConsoleColor(1);
			AppendConsoleMessage(wbuf);
			SetConsoleColor(0);
			free(wbuf);
		}
		free(p);
	}
	dialog->Destroy();
}
*/

void
MyApp::OnOpenConsoleWindow(wxCommandEvent& event)
{
	consoleFrame->Show(true);
	consoleFrame->Raise();
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
MyApp::OnImportAmberLib(wxCommandEvent &event)
{
	MolActionCreateAndPerform(NULL, SCRIPT_ACTION(""), "cmd_import_amberlib");
}

void
MyApp::RegisterScriptMenu(const char *cmd, const char *title)
{
	int i;
	if (cmd[0] == 0 && title[0] == 0)
		i = countScriptMenu;  /*  A sepearator */
	else {
		for (i = 0; i < countScriptMenu; i++) {
			if (strcmp(cmd, scriptMenuCommands[i]) == 0) {
				free(scriptMenuTitles[i]);
				scriptMenuTitles[i] = strdup(title);
				break;
			} else if (strcmp(title, scriptMenuTitles[i]) == 0) {
				free(scriptMenuCommands[i]);
				scriptMenuCommands[i] = strdup(cmd);
				break;
			}
		}
	}
	if (i >= countScriptMenu) {
		if (countScriptMenu == 0) {
			scriptMenuTitles = (char **)malloc(sizeof(char *));
			scriptMenuCommands = (char **)malloc(sizeof(char *));
		} else {
			scriptMenuTitles = (char **)realloc(scriptMenuTitles, sizeof(char *) * (countScriptMenu + 1));
			scriptMenuCommands = (char **)realloc(scriptMenuCommands, sizeof(char *) * (countScriptMenu + 1));
		}
		scriptMenuTitles[countScriptMenu] = strdup(title);
		scriptMenuCommands[countScriptMenu] = strdup(cmd);
		countScriptMenu++;
	}
	if (!scriptMenuModifiedEventPosted) {
		wxCommandEvent myEvent(MyDocumentEvent, MyDocumentEvent_scriptMenuModified);
		wxPostEvent(this, myEvent);	
		scriptMenuModifiedEventPosted = true;
	}
}

void
MyApp::UpdateScriptMenu(wxMenuBar *mbar)
{
	int i;

	wxMenu *smenu = mbar->GetMenu(myMenuIndex_Script);
	if (smenu == NULL)
		return;
	
	//  Remove all custom items
	for (i = smenu->GetMenuItemCount() - 1; i >= countNonCustomScriptMenu; i--) {
		wxMenuItem *item = smenu->FindItemByPosition(i);
		if (!item->IsSeparator()) {
			int n = item->GetId();
			Disconnect(n, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MyApp::OnScriptMenuSelected), NULL, this);
		}
		smenu->Remove(item);
		delete item;
	}
	
	//  Build script menu from internal array
	for (i = 0; i < countScriptMenu; i++) {
		const char *title = scriptMenuTitles[i];
		if (title == NULL || title[0] == 0) {
			smenu->AppendSeparator();
		} else {
			wxString stitle(scriptMenuTitles[i], WX_DEFAULT_CONV);
			wxMenuItem *item = new wxMenuItem(smenu, myMenuID_CustomScript + i, stitle);
			smenu->Append(item);
		}
	}
	Connect(myMenuID_CustomScript, myMenuID_CustomScript + countScriptMenu - 1, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MyApp::OnScriptMenuSelected), NULL, this);
}

void
MyApp::OnScriptMenuModified(wxCommandEvent& event)
{
	scriptMenuModifiedEventPosted = false;
	UpdateScriptMenu(GetMainFrame()->GetMenuBar());
	UpdateScriptMenu(consoleFrame->GetMenuBar());
	event.Skip();
}

void
MyApp::OnScriptMenuSelected(wxCommandEvent& event)
{
	char *cmd;
	int methodType;
	MainView *mview;
	Molecule *mol;
	int index = event.GetId() - myMenuID_CustomScript;
	if (index < 0 || index >= countScriptMenu)
		return;
	cmd = scriptMenuCommands[index];
	methodType = Ruby_methodType("Molecule", cmd);
	if (methodType == 0)
		return;
	mview = MainViewCallback_activeView();
	if (mview == NULL)
		mol = NULL;
	else mol = mview->mol;
	if (methodType == 1 && mol != NULL)  /*  Instance method (with no arguments)  */
		MolActionCreateAndPerform(mol, SCRIPT_ACTION(""), cmd);
	else if (methodType == 2)  /*  Class method (with molecule as an only argument)  */
		MolActionCreateAndPerform(NULL, SCRIPT_ACTION("M"), cmd, mol);
	else return;
}

void
MyApp::OnUpdateUI(wxUpdateUIEvent& event)
{
	int uid = event.GetId();
	MainView *mview = MainViewCallback_activeView();
	if (uid >= myMenuID_CustomScript && uid < myMenuID_CustomScript + countScriptMenu) {
		//  Check the script menu
		char *cmd;
		int methodType;
		Molecule *mol;
		int index = uid - myMenuID_CustomScript;
		cmd = scriptMenuCommands[index];
		methodType = Ruby_methodType("Molecule", cmd);
		event.Enable(false);
		if (methodType != 0) {
			if (mview == NULL)
				mol = NULL;
			else mol = mview->mol;
			if (methodType == 1 && mol != NULL)  /*  Instance method (with no arguments)  */
				event.Enable(true);
			else if (methodType == 2)  /*  Class method (with molecule as an only argument)  */
				event.Enable(true);
		}
	} else {
		switch (uid) {
			case myMenuID_ExecuteScript:
			case myMenuID_OpenConsoleWindow:
			case myMenuID_ViewParameterFilesList:
			case myMenuID_ViewGlobalParameters:
			case myMenuID_MDTools:
			case myMenuID_ImportAmberLib:
				event.Enable(true);
				return;
			default:
				if (mview == NULL)
					event.Enable(false);
				else event.Skip();
		}
	}
}

void
MyApp::OnExecuteScript(wxCommandEvent &event)
{
	wxFileDialog *dialog = new wxFileDialog(NULL, _T("Choose Script File"), _T(""), _T(""), _T("All files (*.*)|*.*"), wxFD_OPEN | wxFD_CHANGE_DIR | wxFD_FILE_MUST_EXIST);
	if (dialog->ShowModal() == wxID_OK) {
		int status;
		wxString path = dialog->GetPath();
		
		//  Command line: execute_script('pathname')
		wxString cline(path);
		wxRegEx re(_T("[\\\\']"));   //  A backslash and a single-quote
		re.Replace(&cline, _T("\\\\\\0"));  //  A backslash followed by "\0"
		cline.Prepend(_T("execute_script('"));
		cline += _T("')");
		MyAppCallback_setConsoleColor(3);
		wxGetApp().AppendConsoleMessage((const char *)(cline.mb_str(wxConvFile)));
		wxGetApp().AppendConsoleMessage("\n");
		MyAppCallback_setConsoleColor(0);

		MyAppCallback_executeScriptFromFile((const char *)(path.mb_str(wxConvFile)), &status);
		if (status != 0)
			Molby_showError(status);
	}
	dialog->Destroy();
}

void
MyApp::OnActivate(wxActivateEvent &event)
{
#if defined(__WXMAC__)
	MyFrame *frame = GetMainFrame();
	frame->Show(false);  /*  Sometimes this "parent" frame gets visible and screw up the menus  */
#endif
	event.Skip();
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

void
MyApp::OnEndProcess(wxProcessEvent &event)
{
	m_processTerminated = true;
	m_processExitCode = event.GetExitCode();
#if LOG_SUBPROCESS
	if (fplog != NULL)
		fprintf(fplog, "OnEndProcess called\n");
#endif
}

int
MyApp::CallSubProcess(const char *cmdline, const char *procname)
{
	const int sEndProcessMessageID = 2;
	int status = 0;
	char buf[256];
	size_t len, len_total;
	wxString cmdstr(cmdline, WX_DEFAULT_CONV);
#if defined(__WXMSW__)
	extern int myKillAllChildren(long pid, wxSignal sig, wxKillError *krc);
#endif
	//  Show progress panel
	if (procname == NULL)
		procname = "subprocess";
	snprintf(buf, sizeof buf, "Running %s...", procname);
	ShowProgressPanel(buf);
	
	//  Create log file in the document home directory
#if LOG_SUBPROCESS
	int nn = 0;
	{
		char *dochome = MyAppCallback_getDocumentHomeDir();
		snprintf(buf, sizeof buf, "%s/%s.log", dochome, procname);
		free(dochome);
		fplog = fopen(buf, "w");
		if (fplog == NULL)
			return -1;
	}
#endif

	//  Create proc object and call subprocess
	wxProcess *proc = new wxProcess(wxGetApp().GetProgressFrame(), sEndProcessMessageID);
	proc->Redirect();
	int flag = wxEXEC_ASYNC;
//#if !__WXMSW__
	flag |= wxEXEC_MAKE_GROUP_LEADER;
//#endif
	m_processTerminated = false;
	m_processExitCode = 0;
	long pid = ::wxExecute(cmdstr, flag, proc);
	if (pid == 0) {
		MyAppCallback_errorMessageBox("Cannot start %s", procname);
		proc->Detach();
		HideProgressPanel();
#if LOG_SUBPROCESS
		fclose(fplog);
#endif
		return -1;
	}
#if LOG_SUBPROCESS
	fprintf(fplog, "[DEBUG]pid = %ld\n", pid);
#endif
	
	//  Wait until process ends or user interrupts
	wxInputStream *in = proc->GetInputStream();
	wxInputStream *err = proc->GetErrorStream();
	len_total = 0;
	while (1) {
		if (m_processTerminated || !wxProcess::Exists(pid)) {
			if (m_processExitCode != 0) {
				/*  Error from subprocess  */
				MyAppCallback_errorMessageBox("%s failed with exit code %d.", procname, m_processExitCode);
				status = m_processExitCode;
			} else status = 0;
			break;
		}
#if defined(__WXMAC__)
		if (waitpid(pid, &status, WNOHANG) != 0) {
			/*  Already finished, although not detected by wxProcess  */
			/*  This sometimes happens on ppc Mac  */
			proc->Detach();
			status = WEXITSTATUS(status);
			break;
		}
#endif
#if LOG_SUBPROCESS
		if (++nn >= 100) {
			fprintf(fplog, "[DEBUG]pid %ld exists\n", pid);
			nn = 0;
		}
#endif
		if (wxGetApp().IsInterrupted()) {
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
				const char *emsg;
				switch (rc) {
					case wxKILL_BAD_SIGNAL: emsg = "no such signal"; break;
					case wxKILL_ACCESS_DENIED: emsg = "permission denied"; break;
					case wxKILL_NO_PROCESS: emsg = "no such process"; break;
					default: emsg = "unknown error"; break;
				}
				MyAppCallback_errorMessageBox("Cannot kill subprocess: %s", emsg);
			}
			proc->Detach();
			status = -2;
			break;
		}
		while (in->CanRead()) {
			in->Read(buf, sizeof buf - 1);
			if ((len = in->LastRead()) > 0) {
				buf[len] = 0;
				len_total += len;
#if LOG_SUBPROCESS
				fprintf(fplog, "%s", buf);
#endif
				MyAppCallback_setConsoleColor(0);
				MyAppCallback_showScriptMessage("%s", buf);
			}
		}
		while (err->CanRead()) {
			err->Read(buf, sizeof buf - 1);
			if ((len = err->LastRead()) > 0) {
				buf[len] = 0;
				len_total += len;
#if LOG_SUBPROCESS
				fprintf(fplog, "%s", buf);
#endif
				MyAppCallback_setConsoleColor(1);
				MyAppCallback_showScriptMessage("\n%s", buf);
				MyAppCallback_setConsoleColor(0); 
			}
		}
	}
#if LOG_SUBPROCESS
	fclose(fplog);
#endif

	HideProgressPanel();
/*	if (len_total > 0)
		MyAppCallback_showRubyPrompt(); */
	return status;
}

#pragma mark ====== MyFrame (top-level window) ======

/*
 * This is the top-level window of the application.
 */
 
IMPLEMENT_CLASS(MyFrame, wxDocMDIParentFrame)
BEGIN_EVENT_TABLE(MyFrame, wxDocMDIParentFrame)
    EVT_MENU(wxID_ABOUT, MyFrame::OnAbout)
END_EVENT_TABLE()

MyFrame::MyFrame(wxDocManager *manager, wxFrame *frame, const wxString& title,
    const wxPoint& pos, const wxSize& size, long type):
  wxDocMDIParentFrame(manager, frame, wxID_ANY, title, pos, size, type, _T("myFrame"))
{
	editMenu = (wxMenu *) NULL;
#if defined(__WXMAC__)
	/*  Avoid this "dummy" top-level window to appear in the window menu.
	    It should not happen because MyApp::OnActivate() tries to hide this window,
	    but this is still here just in case.  */
	OSStatus sts;
	sts = ChangeWindowAttributes((WindowRef)m_macWindow, 0, kWindowInWindowMenuAttribute);
/*	printf("m_macWindow = %p, status = %d\n", m_macWindow, (int)sts); */
#endif
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event) )
{
	extern const char *gVersionString, *gCopyrightString;
	char *s;
	asprintf(&s, 
			 "Molby %s\n%s\n%s\n"
			 "Including:\n"
			 "AmberTools 1.3, http://ambermd.org/\n"
			 "  Copyright (c) Junmei Wang, Ross C. Walker, \n"
			 "  Michael F. Crowley, Scott Brozell and David A. Case\n"
			 "wxWidgets %d.%d.%d, Copyright (c) 1992-2008 Julian Smart, \n"
			 "  Robert Roebling, Vadim Zeitlin and other members of the \n"
			 "  wxWidgets team\n"
			 "  Portions (c) 1996 Artificial Intelligence Applications Institute\n"
			 "ruby %s\n%s",
			 gVersionString, gCopyrightString, sLastBuildString,
			 wxMAJOR_VERSION, wxMINOR_VERSION, wxRELEASE_NUMBER,
			 gRubyVersion, gRubyCopyright);
	wxString str(s, WX_DEFAULT_CONV);
    (void)wxMessageBox(str, _T("Molby"));
}

MyFrame *GetMainFrame(void)
{
	return frame;
}

#pragma mark ====== Plain-C interface ======

void
MyAppCallback_loadGlobalSettings(void)
{
	wxGetApp().LoadDefaultSettings();
}

void
MyAppCallback_saveGlobalSettings(void)
{
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
	wxString wxkey(key, WX_DEFAULT_CONV);
	wxString wxvalue = wxGetApp().GetDefaultSetting(wxkey);
	return strdup(wxvalue.mb_str(WX_DEFAULT_CONV));
}

void
MyAppCallback_setGlobalSettings(const char *key, const char *value)
{
	wxString wxkey(key, WX_DEFAULT_CONV);
	wxString wxvalue(value, WX_DEFAULT_CONV);
	wxGetApp().SetDefaultSetting(wxkey, wxvalue);
}

int
MyAppCallback_getGlobalSettingsWithType(const char *key, int type, void *ptr)
{
	int retval;
	char *s = MyAppCallback_getGlobalSettings(key);
	char desc[] = SCRIPT_ACTION("s; ");
	desc[sizeof(desc) - 2] = type;
	retval = MolActionCreateAndPerform(NULL, desc, "eval", s, ptr);
	free(s);
	return retval;
}

int
MyAppCallback_setGlobalSettingsWithType(const char *key, int type, const void *ptr)
{
	const char *cmd = "set_global_settings";
	switch (type) {
		case 'i': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("i"), cmd, *((const Int *)ptr));
		case 'd': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("d"), cmd, *((const Double *)ptr));
		case 's': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("s"), cmd, (const char *)ptr);
		case 'v': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("v"), cmd, (const Vector *)ptr);
		case 't': return MolActionCreateAndPerform(NULL, SCRIPT_ACTION("t"), cmd, (const Transform *)ptr);
		default:
			MyAppCallback_errorMessageBox("Internal error: unsupported format '%c' at line %d, file %s", type, __LINE__, __FILE__);
			return -2;
	}
}

int
MyAppCallback_showScriptMessage(const char *fmt, ...)
{
	if (fmt != NULL) {
		char *p;
		va_list ap;
		int retval;
		va_start(ap, fmt);
		if (strchr(fmt, '%') == NULL) {
			/*  No format characters  */
			return wxGetApp().AppendConsoleMessage(fmt);
		} else if (strcmp(fmt, "%s") == 0) {
			/*  Direct output of one string  */
			p = va_arg(ap, char *);
			return wxGetApp().AppendConsoleMessage(p);
		}
#if 1
		vasprintf(&p, fmt, ap);
#else
		/*  Use safe wxString method  */
		/*  Not necessary any longer; vasprintf() is implemented in Missing.c  */
		{
			wxString str;
			str.PrintfV(wxString::FromUTF8(fmt).GetData(), ap);
			p = strdup((const char *)str.mb_str(WX_DEFAULT_CONV));
		}
#endif
		if (p != NULL) {
			retval = wxGetApp().AppendConsoleMessage(p);
			free(p);
			return retval;
		} else return 0;
	} else {
		wxGetApp().FlushConsoleMessage();
		return 0;
	}
  return 0;
}

void
MyAppCallback_setConsoleColor(int color)
{
	wxGetApp().SetConsoleColor(color);
}

void
MyAppCallback_showRubyPrompt(void)
{
	MyAppCallback_setConsoleColor(0);
	MyAppCallback_showScriptMessage("%% ");
}

int
MyAppCallback_checkInterrupt(void)
{
	return wxGetApp().IsInterrupted();
}

void
MyAppCallback_showProgressPanel(const char *msg)
{
	wxGetApp().ShowProgressPanel(msg);
}

void
MyAppCallback_hideProgressPanel(void)
{
	wxGetApp().HideProgressPanel();
}

void
MyAppCallback_setProgressValue(double dval)
{
	wxGetApp().SetProgressValue(dval);
}

void
MyAppCallback_setProgressMessage(const char *msg)
{
	wxGetApp().SetProgressMessage(msg);
}

int
MyAppCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize)
{
	wxDialog *dialog = new wxDialog(NULL, -1, _T("Input request"), wxDefaultPosition);
	wxStaticText *stext;
	wxTextCtrl *tctrl;
	int retval;
	wxString pstr(prompt, WX_DEFAULT_CONV);
	{	//  Vertical sizer containing [prompt, textbox, buttons]
		wxBoxSizer *sizer1;
		sizer1 = new wxBoxSizer(wxVERTICAL);
		stext = new wxStaticText(dialog, -1, pstr, wxDefaultPosition, wxSize(200, 22));
		sizer1->Add(stext, 0, wxEXPAND | wxALL, 6);
		tctrl = new wxTextCtrl(dialog, -1, _T(""), wxDefaultPosition, wxSize(200, 22));
		sizer1->Add(tctrl, 0, wxEXPAND | wxALL, 6);
		wxSizer *bsizer = dialog->CreateButtonSizer(wxOK | wxCANCEL);
		sizer1->Add(bsizer, 0, wxEXPAND | wxALL, 6);
		sizer1->Layout();
		dialog->SetSizerAndFit(sizer1);
		dialog->Centre(wxBOTH);
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
	int wxflags, wxicon, retval;
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

char *
MyAppCallback_getDocumentHomeDir(void)
{
	char *s;
#if __WXMSW__
	char *ss;
	s = getenv("USERPROFILE");
	asprintf(&ss, "%s\\My Documents", s);
	return ss;
#else
	s = getenv("HOME");
	return (s == NULL ? NULL : strdup(s));
#endif
}

void
MyAppCallback_registerScriptMenu(const char *cmd, const char *title)
{
	wxGetApp().RegisterScriptMenu(cmd, title);
}


RubyValue
MyAppCallback_executeScriptFromFile(const char *cpath, int *status)
{
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
	wxFile file;
	if (file.Open((const wxChar *)path, wxFile::read)) {
		wxFileOffset len = file.Length();
		script = (char *)malloc(len + 1);
		if (script != NULL) {
			file.Read(script, len);
			script[len] = 0;
		}
	}
	file.Close();
	
	retval = Molby_evalRubyScriptOnMolecule(script, MoleculeCallback_currentMolecule(), status);
	free(script);
	free(p);
	wxFileName::SetCwd(cwd);
	return retval;
}

void MyAppCallback_beginUndoGrouping(void)
{
	wxList &doclist = wxGetApp().DocManager()->GetDocuments();
	wxList::iterator iter;
	for (iter = doclist.begin(); iter != doclist.end(); ++iter) {
		((MyDocument *)(*iter))->BeginUndoGrouping();
	}
}

void MyAppCallback_endUndoGrouping(void)
{
	wxList &doclist = wxGetApp().DocManager()->GetDocuments();
	wxList::iterator iter;
	for (iter = doclist.begin(); iter != doclist.end(); ++iter) {
		((MyDocument *)(*iter))->EndUndoGrouping();
	}
}

int MyAppCallback_callSubProcess(const char *cmdline, const char *procname)
{
	return wxGetApp().CallSubProcess(cmdline, procname);
}
