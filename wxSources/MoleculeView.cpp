/*
 *  MoleculeView.cpp
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

#include "MoleculeView.h"

#include "MyApp.h"
#include "MyDocument.h"
#include "MyGLCanvas.h"
#include "MyCommand.h"
#include "MySlider.h"
#include "MyListCtrl.h"
#include "../MolLib/Missing.h"
#include "../MolLib/Ruby_bind/Molby_extern.h"
#include "MyMBConv.h"
#include "MyProgressIndicator.h"

#include "wx/tglbtn.h"
#include "wx/listctrl.h"
#include "wx/splitter.h"
#include "wx/choice.h"
#include "wx/font.h"

#if defined(__WXMSW__)
#include "OpenGL_extensions.h"
#endif

//#include "../MolLib/Ruby_bind/Molby_extern.h"

enum {
	myID_RotButton = 500,
	myID_TransButton,
	myID_ScaleButton,
	myID_SelectButton,
	myID_BondButton,
	myID_EraseButton,
	myID_RotateBondSlider,
	myID_RotateXSlider,
	myID_RotateYSlider,
	myID_FrameControlPanel,
	myID_FrameSlider,
	myID_JumpToStartButton,
	myID_PlayBackwardButton,
	myID_FrameText,
	myID_PlayForwardButton,
	myID_JumpToEndButton,
	myID_Table,
	myID_TableMenu,
	myID_StopProgressButton
};

IMPLEMENT_DYNAMIC_CLASS(MoleculeView, wxView)

BEGIN_EVENT_TABLE(MoleculeView, wxView)
	EVT_TOGGLEBUTTON(myID_RotButton, MoleculeView::OnButtonPressed)
	EVT_TOGGLEBUTTON(myID_TransButton, MoleculeView::OnButtonPressed)
	EVT_TOGGLEBUTTON(myID_ScaleButton, MoleculeView::OnButtonPressed)
	EVT_TOGGLEBUTTON(myID_SelectButton, MoleculeView::OnButtonPressed)
	EVT_TOGGLEBUTTON(myID_BondButton, MoleculeView::OnButtonPressed)
	EVT_TOGGLEBUTTON(myID_EraseButton, MoleculeView::OnButtonPressed)
	EVT_BUTTON(myID_StopProgressButton, MoleculeView::OnStopProgressPressed)
	EVT_COMMAND(wxID_ANY, MySliderEvent, MoleculeView::OnSliderAction)
	EVT_COMMAND_SCROLL(myID_FrameSlider, MoleculeView::OnFrameSliderAction)
	EVT_TEXT_ENTER(myID_FrameText, MoleculeView::OnFrameTextAction)
	EVT_CHOICE(myID_TableMenu, MoleculeView::OnSelectTable)
	EVT_ACTIVATE(MoleculeView::OnActivate)
END_EVENT_TABLE()
#define ConnectMouseDownEvents(src, func, target) \
	(src->Connect(wxEVT_LEFT_DOWN, wxMouseEventHandler(func), NULL, target), \
	src->Connect(wxEVT_LEFT_DCLICK, wxMouseEventHandler(func), NULL, target))
#define DisconnectMouseDownEvents(src, func, target) \
  (src->Disconnect(wxEVT_LEFT_DOWN, wxMouseEventHandler(func), NULL, target), \
  src->Disconnect(wxEVT_LEFT_DCLICK, wxMouseEventHandler(func), NULL, target))

WX_DEFINE_ARRAY_PTR(MoleculeView *, ArrayOfMoleculeViews);

ArrayOfMoleculeViews sActiveViews;

MoleculeView::~MoleculeView()
{
  /*  Unlink the connection with the widgets (otherwise the widgets may try to
   access to this view after deletion  */
  canvas->Disconnect(-1, wxEVT_CHAR, wxKeyEventHandler(MoleculeView::OnChar), NULL, this);
  DisconnectMouseDownEvents(bbuttons[0], MoleculeView::OnFrameButtonAction, this);
  DisconnectMouseDownEvents(bbuttons[1], MoleculeView::OnFrameButtonAction, this);
  DisconnectMouseDownEvents(bbuttons[2], MoleculeView::OnFrameButtonAction, this);
  DisconnectMouseDownEvents(bbuttons[3], MoleculeView::OnFrameButtonAction, this);

  //  Disconnect the notification handler
  MolDocument()->Disconnect(MyDocumentEvent_documentModified, MyDocumentEvent, wxCommandEventHandler(MoleculeView::OnDocumentModified), NULL, this);
  wxGetApp().Disconnect(MyDocumentEvent_scriptMenuModified, MyDocumentEvent, wxCommandEventHandler(MoleculeView::OnScriptMenuModified), NULL, this);

  //  Unbind the double-click handler of MyListCtrl
  listctrl->GetScrolledWindow()->Unbind(wxEVT_LEFT_DCLICK, &MoleculeView::OnLeftDClickInListCtrl, this);

  //  Unset data source for the list control
  listctrl->SetDataSource(NULL);

}

bool
MoleculeView::OnCreate(wxDocument *doc, long WXUNUSED(flags) )
{
	int i;

	const wxFont *ctrlFont;
#if __WXOSX_COCOA__
	ctrlFont = wxSMALL_FONT;
#else
	ctrlFont = NULL;
#endif
	
	// Make a document frame
	frame = new wxDocChildFrame(doc, this, GetMainFrame(), wxID_ANY, _T("New Molby Document"),
						   wxDefaultPosition, wxDefaultSize,
						   wxDEFAULT_FRAME_STYLE |
						   wxNO_FULL_REPAINT_ON_RESIZE);
    frame->SetPosition(FromFrameDIP(frame, wxPoint(10, 24)));
    frame->SetClientSize(FromFrameDIP(frame, wxSize(680, 400)));
	canvas = NULL;
	mview = NULL;
	listmenu = NULL;
	listctrl = NULL;
	file_history_menu = NULL;
	edit_menu = NULL;
	memset(tbuttons, 0, sizeof(tbuttons));
	infotext = NULL;
	frameControlPanel = NULL;
	frameSlider = NULL;
	frameText = NULL;
	isRebuildingTable = false;

	Molecule *mol = ((MyDocument *)doc)->GetMolecule();
	if (mol != NULL) {
		mview = mol->mview;
		MainView_setViewObject(mview, this);
	}

	wxMenuBar *menu_bar = wxGetApp().CreateMenuBar(1, &file_history_menu, &edit_menu);
	
	// Associate the menu bar with the frame
	frame->SetMenuBar(menu_bar);

	// Associate the edit menu with the command processor
	doc->GetCommandProcessor()->SetEditMenu(edit_menu);
	
	// Create the window content
	
	//  A splitter window embraces a grid (left) and main screen (right)
	wxSplitterWindow *splitter = new wxSplitterWindow(frame, -1, wxDefaultPosition, wxDefaultSize, wxSP_3D | wxSP_LIVE_UPDATE);
	splitter->SetMinimumPaneSize(1);
	
	//  Create the left half
	//  A panel containing a popup menu and a list window
	wxPanel *panel0 = new wxPanel(splitter);
	{
		char buf[32];
		wxBoxSizer *sizer0;
		sizer0 = new wxBoxSizer(wxVERTICAL);
		wxArrayString choiceItems;
		for (i = 0; ; i++) {
			MainView_tableTitleForIndex(mview, i, buf, sizeof buf);
			if (buf[0] == 0)
				break;
			wxString itemTitle(buf, WX_DEFAULT_CONV);
			choiceItems.Add(itemTitle);
		}
		
		listmenu = new wxChoice(panel0, myID_TableMenu, wxDefaultPosition, wxDefaultSize, choiceItems);
		sizer0->Add(listmenu, 0, wxALL, 0);
		
		listctrl = new MyListCtrl();
		listctrl->Create(panel0, myID_Table, wxDefaultPosition, wxDefaultSize);
		sizer0->Add(listctrl, 1, wxALL | wxEXPAND, 0);
		panel0->SetSizer(sizer0);		
	}
		
	//  Create the right half
	//  A panel containing MyGLCanvas, buttons, sliders, etc.
	wxPanel *panel1 = new wxPanel(splitter);
	
	{	//  Vertical sizer containing [sizer2, sizer3, sizer4]
		wxBoxSizer *sizer1;
		sizer1 = new wxBoxSizer(wxVERTICAL);
		
		{	//  Horizontal sizer containing [button0, ..., button5, infotext, progress_indicator]
			wxBoxSizer *sizer2;
			sizer2 = new wxBoxSizer(wxHORIZONTAL);
			
			{	// Button0..5 (Rot/Trans/Scale/Select/Bond/Erase)
				wxString labels[] = {
					wxT("Rot"), wxT("Trans"), wxT("Scale"), wxT("Select"), wxT("Bond"), wxT("Erase")
				};
				wxWindowID ids[] = {
					myID_RotButton, myID_TransButton, myID_ScaleButton, 
					myID_SelectButton, myID_BondButton, myID_EraseButton
				};
				for (i = 0; i < 6; i++) {
					tbuttons[i] = new MyToggleButton(panel1, ids[i], labels[i], wxDefaultPosition, FromFrameDIP(frame, wxSize(40, 32)), wxTOGGLEBUTTON_STYLE);
					if (ctrlFont)
						tbuttons[i]->SetFont(*ctrlFont);
					sizer2->Add(tbuttons[i], 0, wxALL | wxEXPAND, FromFrameDIP(frame, 3));
				}
				tbuttons[0]->SetValue(true);
			}
			{	// Information text
				infotext = new wxStaticText(panel1, -1, wxT(""), wxDefaultPosition, FromFrameDIP(frame, wxSize(80, 32)), wxST_NO_AUTORESIZE | wxBORDER_SUNKEN);
				infotext->SetMinSize(FromFrameDIP(frame, wxSize(80, 32)));
#if defined(__WXMSW__)
				infotext->SetFont(*wxSMALL_FONT);
#else
				if (ctrlFont)
					infotext->SetFont(*ctrlFont);
#endif
				sizer2->Add(infotext, 1, wxALL | wxEXPAND, 3);   // Can expand horizontally
			}
			{	// Custom progress indicator
				progress = new MyProgressIndicator(panel1, myID_StopProgressButton, wxDefaultPosition, FromFrameDIP(frame, wxSize(12, 24)));
				sizer2->Add(progress, 0, wxALL | wxEXPAND, 3);
			}
			sizer1->Add(sizer2, 0, wxALL | wxEXPAND, 0);
		}
		
		{	// Horizontal sizer containing [sizer31, sizer32, sizer33]
			wxBoxSizer *sizer3 = new wxBoxSizer(wxHORIZONTAL);
			
			{	// Vertical sizer containing [button, mySlider]
				wxBoxSizer *sizer31 = new wxBoxSizer(wxVERTICAL);
				{	// "Rotate bond" button and mySlider
					#include "../bitmaps/rotate_bond.xpm"
					wxBitmap bmp1(rotate_bond_xpm);
					wxBitmapButton *button1 = new wxBitmapButton(panel1, -1, bmp1, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, 21)), wxTOGGLEBUTTON_STYLE);
					sizer31->Add(button1, 0, 0, 0);
					button1->Disable();
					MySlider *slider1 = new MySlider(panel1, myID_RotateBondSlider, wxVERTICAL, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, 21)));
					sizer31->Add(slider1, 1, wxEXPAND);
				}
				sizer3->Add(sizer31, 0, wxALL | wxEXPAND, 0);
			}
			
			{	// Vertical sizer containing [Canvas, [button, mySlider]]
				wxBoxSizer *sizer32 = new wxBoxSizer(wxVERTICAL);
				{
					canvas = new MyGLCanvas(this, panel1, wxDefaultPosition, FromFrameDIP(frame, wxSize(100, 100)));
					sizer32->Add(canvas, 1, wxALL | wxEXPAND, 0);
					
					//  Let the MyGLCanvas pass the keyboard event to this
					canvas->Connect(-1, wxEVT_CHAR, wxKeyEventHandler(MoleculeView::OnChar), NULL, this);
                }
				{
					wxBoxSizer *sizer321 = new wxBoxSizer(wxHORIZONTAL);
					{
						#include "../bitmaps/rotate_y.xpm"
						wxBitmap bmp2(rotate_y_xpm);
						wxBitmapButton *button2 = new wxBitmapButton(panel1, -1, bmp2, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, 21)), wxTOGGLEBUTTON_STYLE);
						sizer321->Add(button2, 0, 0, 0);
						button2->Disable();
						MySlider *slider2 = new MySlider(panel1, myID_RotateYSlider, wxHORIZONTAL, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, 21)));
						sizer321->Add(slider2, 1, wxEXPAND);
					}
					sizer32->Add(sizer321, 0, wxEXPAND);
				}
				sizer3->Add(sizer32, 1, wxEXPAND);
			}

			{	// Vertical sizer containing [button, mySlider]
				wxBoxSizer *sizer33 = new wxBoxSizer(wxVERTICAL);
				{	// "Rotate bond" button and mySlider
					#include "../bitmaps/rotate_x.xpm"
					wxBitmap bmp3(rotate_x_xpm);
					wxBitmapButton *button3 = new wxBitmapButton(panel1, -1, bmp3, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, 21)), wxTOGGLEBUTTON_STYLE);
					button3->Disable();
					sizer33->Add(button3, 0, 0, 0);
					
					MySlider *slider3 = new MySlider(panel1, myID_RotateXSlider, wxVERTICAL, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, 21)));
					sizer33->Add(slider3, 1, wxEXPAND);
				}
				sizer3->Add(sizer33, 0, wxALL | wxEXPAND, 0);
			}
			
			sizer1->Add(sizer3, 1, wxALL | wxEXPAND, 0);
		}
		
		{	//  Horizontal sizer containing frame controls
			
			const int height = 18;
			frameControlPanel = new wxPanel(panel1, myID_FrameControlPanel, wxDefaultPosition, FromFrameDIP(frame, wxSize(200, height)));
			wxBoxSizer *sizer4 = new wxBoxSizer(wxHORIZONTAL);
			{
				frameSlider = new wxSlider(frameControlPanel, myID_FrameSlider, 0, 0, 1, wxDefaultPosition, FromFrameDIP(frame, wxSize(40, height - 2)));
				frameSlider->Enable(false);
				sizer4->Add(frameSlider, 1, wxALL | wxEXPAND, 1);
			
				#include "../bitmaps/jump_to_start.xpm"
				wxBitmap bmp41(jump_to_start_xpm);
				wxBitmapButton *button41 = new wxBitmapButton(frameControlPanel, myID_JumpToStartButton, bmp41, wxDefaultPosition, FromFrameDIP(frame, wxSize(16, height)), wxTOGGLEBUTTON_STYLE);
				sizer4->Add(button41, 0, wxEXPAND);
        bbuttons[0] = button41;
				ConnectMouseDownEvents(bbuttons[0], MoleculeView::OnFrameButtonAction, this);

				#include "../bitmaps/play_backward.xpm"
				wxBitmap bmp42(play_backward_xpm);
				wxBitmapButton *button42 = new wxBitmapButton(frameControlPanel, myID_PlayBackwardButton, bmp42, wxDefaultPosition, FromFrameDIP(frame, wxSize(16, height)), wxTOGGLEBUTTON_STYLE);
				sizer4->Add(button42, 0, wxEXPAND);
        bbuttons[1] = button42;
				ConnectMouseDownEvents(bbuttons[1], MoleculeView::OnFrameButtonAction, this);
				
				{
					frameText = new wxTextCtrl(frameControlPanel, myID_FrameText, wxT(""), wxDefaultPosition, FromFrameDIP(frame, wxSize(40, height)));
					if (ctrlFont) {
						wxTextAttr attr(*wxBLACK, wxNullColour, *ctrlFont);
						frameText->SetDefaultStyle(attr);
						frameText->SetFont(*ctrlFont);
					}
					sizer4->Add(frameText, 0, wxEXPAND);
				}
			
				#include "../bitmaps/play_forward.xpm"
				wxBitmap bmp43(play_forward_xpm);
				wxBitmapButton *button43 = new wxBitmapButton(frameControlPanel, myID_PlayForwardButton, bmp43, wxDefaultPosition, FromFrameDIP(frame, wxSize(16, height)), wxTOGGLEBUTTON_STYLE);
				sizer4->Add(button43, 0, wxEXPAND);
        bbuttons[2] = button43;
				ConnectMouseDownEvents(bbuttons[2], MoleculeView::OnFrameButtonAction, this);

				#include "../bitmaps/jump_to_end.xpm"
				wxBitmap bmp44(jump_to_end_xpm);
				wxBitmapButton *button44 = new wxBitmapButton(frameControlPanel, myID_JumpToEndButton, bmp44, wxDefaultPosition, FromFrameDIP(frame, wxSize(16, height)), wxTOGGLEBUTTON_STYLE);
				sizer4->Add(button44, 0, wxEXPAND);
        bbuttons[3] = button44;
				ConnectMouseDownEvents(bbuttons[3], MoleculeView::OnFrameButtonAction, this);
				
				wxPanel *spacer = new wxPanel(frameControlPanel, -1, wxDefaultPosition, FromFrameDIP(frame, wxSize(21, height)));
				sizer4->Add(spacer, 0, wxEXPAND);
			}
			frameControlPanel->SetSizer(sizer4);
			sizer1->Add(frameControlPanel, 0, wxALL | wxEXPAND, 0);
		//	controls->Disable();
		}

		panel1->SetSizer(sizer1);
	}

	splitter->SplitVertically(panel0, panel1);

	wxBoxSizer *mainsizer = new wxBoxSizer(wxHORIZONTAL);
	mainsizer->Add(splitter, 1, wxEXPAND);
	frame->SetSizer(mainsizer);

	mainsizer->Layout();
	splitter->SetSashPosition(FromFrameDIP(frame, 240), true);

	//  Initialize table view
	MainView_createColumnsForTableAtIndex(mview, 0);

	//  Select table popup
	listmenu->SetSelection(0);
	
#if defined(__X__) || defined(__WXMAC__)
    // X seems to require a forced resize
    int x, y;
    frame->GetSize(&x, &y);
    frame->SetSize(wxDefaultCoord, wxDefaultCoord, x, y);
#endif
    frame->Show(true);
    Activate(true);

	//  Initial keyboard focus is on the GL canvas (to accept 'S' for scale, etc.)
	canvas->SetFocus();
	
	//  Connect the notification handler
	doc->Connect(MyDocumentEvent_documentModified, MyDocumentEvent, wxCommandEventHandler(MoleculeView::OnDocumentModified), NULL, this);

	wxGetApp().Connect(MyDocumentEvent_scriptMenuModified, MyDocumentEvent, wxCommandEventHandler(MoleculeView::OnScriptMenuModified), NULL, this);

	//  Intercept the double-click handler of MyListCtrl
	listctrl->GetScrolledWindow()->Bind(wxEVT_LEFT_DCLICK, &MoleculeView::OnLeftDClickInListCtrl, this);

	//  Set data source for the list control
	listctrl->SetDataSource(this);
				
#if defined(__WXMSW__)
	//  Initialize OpenGL extension
	{
		static int openGLExtension_inited = 0;
		if (openGLExtension_inited == 0) {
			if (InitializeOpenGLExtensions() != 0) {
				MyAppCallback_errorMessageBox("Fatal internal error: cannot initialize OpenGL extensions");
				openGLExtension_inited = -1;
			} else openGLExtension_inited = 1;
		}
	}
	
#endif
	
    return true;
}

void
MoleculeView::OnDraw(wxDC *dc)
{
	if (mview != NULL) {
		if (mview->isPrinting) {
			float scale = 4.0;
			wxImage *img = CaptureGLCanvas(scale);
			if (img != NULL) {
				wxBitmap bitmap(*img);
				double sx, sy;
				dc->GetUserScale(&sx, &sy);
				dc->SetUserScale(sx / scale, sy / scale);
				dc->DrawBitmap(bitmap, 0, 0);
				delete img;
			}
			mview->isPrinting = 0;
		} else {
			if (mview->mol != NULL) {
				MoleculeLock(mview->mol);
				MainView_drawModel(mview);
				MoleculeUnlock(mview->mol);
			}
		}
	}
}

wxImage *
MoleculeView::CaptureGLCanvas(float scale, int bg_color, int width, int height)
{
	if (canvas && mview->mol != NULL) {
		int x, y, cwidth, cheight;
		float bgcol[4], rx, ry;
		
		canvas->SetCurrent();
		canvas->GetClientSize(&cwidth, &cheight);
		if (width <= 0)
			width = cwidth;
		if (height <= 0)
			height = cheight;
		rx = (float)width / cwidth;
		ry = (float)height / cheight;
		if (rx > ry)
			rx = ry;
		width *= scale;
		height *= scale;

		//  Create OpenGL offscreen buffer
		GLuint frame_buf, render_buf, depth_buf;
		frame_buf = render_buf = depth_buf = 0;
		glGenFramebuffersEXT(1, &frame_buf);
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frame_buf);
		
		glGenRenderbuffersEXT(1, &render_buf);
		glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, render_buf);
		glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA, width, height);

		glGenRenderbuffersEXT(1, &depth_buf);
		glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depth_buf);
		glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, width, height);
		
		glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, render_buf);
		glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depth_buf);

		MainView_initializeOpenGL();

		mview->offline_scale = 1.0;  /*  Scale is handled in offline_width and offline_height  */
		mview->offline_width = width;
		mview->offline_height = height;
		
		for (x = 0; x < 4; x++) {
			bgcol[x] = mview->background_color[x];
		}
		if (bg_color == 0) {
			mview->background_color[0] = 0;
			mview->background_color[1] = 0;
			mview->background_color[2] = 0;
			mview->background_color[3] = 0;
		} else if (bg_color == 1) {
			mview->background_color[0] = 0;
			mview->background_color[1] = 0;
			mview->background_color[2] = 0;
			mview->background_color[3] = 1;
		} else if (bg_color == 2) {
			mview->background_color[0] = 1;
			mview->background_color[1] = 1;
			mview->background_color[2] = 1;
			mview->background_color[3] = 1;
		}
		MoleculeLock(mview->mol);
		MainView_drawModel(mview);
		MoleculeUnlock(mview->mol);
		for (x = 0; x < 4; x++) {
			mview->background_color[x] = bgcol[x];
		}
		mview->offline_scale = 0.0;
		mview->offline_width = 0;
		mview->offline_height = 0;

		glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
		unsigned char *glBitmapData = (unsigned char *)malloc(4 * width * height);
		unsigned char *glRGBData = (unsigned char *)malloc(3 * width * height);
		unsigned char *glAlphaData = (unsigned char *)malloc(width * height);
		glReadPixels((GLint)0, (GLint)0, (GLint)width, (GLint)height, GL_RGBA, GL_UNSIGNED_BYTE, glBitmapData);
		for (y = 0; y < height; y++) {
			unsigned char *p1 = glBitmapData + 4 * y * width;
			unsigned char *p2 = glRGBData + 3 * (height - 1 - y) * width;
			unsigned char *p3 = glAlphaData + (height - 1 - y) * width;
			for (x = 0; x < width; x++) {
				//  Copy RGB data and Alpha data separately, with flipping vertically
				*p2++ = *p1++;
				*p2++ = *p1++;
				*p2++ = *p1++;
				*p3++ = *p1++;
			}
		}
		wxImage *img = new wxImage(width, height, glRGBData, glAlphaData);
		free(glBitmapData);  // glRGBData and glAlphaData are deallocated within wxImage constructor
		
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		glDeleteFramebuffersEXT(1, &frame_buf);
		glDeleteRenderbuffersEXT(1, &render_buf);
		
		canvas->Refresh();
	
		return img;
	}
	return NULL;
}

int
MoleculeView::DoExportGraphic(wxString& fname, float scale, int bg_color, int width, int height)
{
	wxImage *img = CaptureGLCanvas(scale, bg_color, width, height);
	if (img == NULL)
		return -1;
	wxString ext = fname.AfterLast('.');
	wxBitmapType type = wxBITMAP_TYPE_PNG;
	if (ext.CmpNoCase(_T("tif")) == 0)
		type = wxBITMAP_TYPE_TIF;
	img->SaveFile(fname, type);
	delete img;
	return 0;
}

void
MoleculeView::OnUpdate(wxView *WXUNUSED(sender), wxObject *WXUNUSED(hint))
{
	if (canvas)
		canvas->Refresh();

/*  Maybe necessary in some platforms (not in MacOSX and MSW)  */
#if 0
  if (canvas) {
      wxClientDC dc(canvas);
      dc.Clear();
      OnDraw(&dc);
  }
#endif
}

bool
MoleculeView::OnClose(bool deleteWindow)
{
//#if !defined(__WXMAC__) && !defined(__WXOSX__)
	//  On wxOSX, this causes invocation of MyDocument::Close() twice, which
	//  apprently is not very good. However, on wxMSW this is not the case.
	//  So we need to keep this code for wxMSW but not for wxOSX.
	if (!GetDocument()->Close())
		return false;
//#endif

    wxDocChildFrame *saveFrame = frame;

	//  Dispose relationship between this and Molecule (MainView)
	MainView_setViewObject(mview, NULL);
	mview = NULL;
    
	//  Remove this from the active view list
	sActiveViews.Remove(this);

	//  Dispose Connection between DocManager and file history menu
	wxGetApp().DocManager()->FileHistoryRemoveMenu(file_history_menu);

	wxGetApp().Disconnect(MyDocumentEvent_scriptMenuModified, MyDocumentEvent, wxCommandEventHandler(MoleculeView::OnScriptMenuModified), NULL, this);

	if (deleteWindow) {
		saveFrame->Destroy();
    } else saveFrame->Hide();
	
	//  Check if all windows are gone
	wxGetApp().CheckIfAllWindowsAreGone(saveFrame);

	return true;
}

void
MoleculeView::Activate(bool activate)
{
	if (activate) {
		int i, n = sActiveViews.GetCount();
		for (i = 0; i < n; i++) {
			if (sActiveViews.Item(i) == this)
				break;
		}
		if (i < n) {
			/*  Remove this if existing  */
			sActiveViews.RemoveAt(i);
		}
		/*  Insert self as the first item  */
		sActiveViews.Insert(this, 0);
	}
	wxView::Activate(activate);
    if (frame != NULL)
        frame->Refresh();
}

wxPrintout *
MoleculeView::OnCreatePrintout()
{
	if (mview != NULL)
		mview->isPrinting = 1;
	return wxView::OnCreatePrintout();
}

void
MoleculeView::UpdateFrameControlValues()
{
	if (mview == NULL || mview->mol == NULL)
		return;
	wxString str;
	int cframe;
	cframe = mview->mol->cframe;
	str.Printf(_T("%d"), cframe);
//	frameText->SetValue(str);
	frameText->Clear();
	frameText->AppendText(str);
	frameSlider->SetValue(cframe);
}

void
MoleculeView::UpdateFrameControls()
{
	int nframes;
	bool enabled = false;
	if (mview != NULL && mview->mol != NULL) {
		MoleculeLock(mview->mol);
		nframes = MoleculeGetNumberOfFrames(mview->mol);
		MoleculeUnlock(mview->mol);
		if (nframes > 1)
			enabled = true;
	}
	
	frameControlPanel->Enable(enabled);
	if (enabled) {
		frameSlider->Enable(true);
		frameSlider->SetRange(0, nframes - 1);
		UpdateFrameControlValues();
	} else {
		frameSlider->Enable(false);
		frameSlider->SetRange(0, 1);
		frameText->SetValue(_T(""));
	}
	frameControlPanel->Update();
}

void
MoleculeView::InvalidateProgressIndicator()
{
	progress->Refresh();
}

void
MoleculeView::ProceedProgressIndicator()
{
	progress->ProceedIndicatorState();
}

void
MoleculeView::SelectButtonForMode(int mode)
{
	int i;
	if (mode >= 1 && mode <= 6) {
		for (i = 0; i < 6; i++) {
			tbuttons[i]->SetValue((i == mode - 1));
		}
	}
	MainViewCallback_setKeyboardFocus(mview);
}

void
MoleculeView::OnButtonPressed(wxCommandEvent& event)
{
	int eventId = event.GetId();
	int mode = eventId - myID_RotButton + kTrackballRotateMode;
	MainView_setMode(mview, mode);
	SelectButtonForMode(mode);
}

void
MoleculeView::OnSliderAction(wxCommandEvent& event)
{
	int eventId = event.GetId();
	int mode = eventId - myID_RotateBondSlider + 1;
    char buf[256];
	MySlider *sender = (MySlider *)event.GetEventObject();
	float angle = sender->GetFloatValue();
    float angledeg = angle * 360.0;
    float fixangledeg = floor(angledeg / 15.0 + 0.5) * 15.0;  //  Nearest multiple of 15 deg
    if (angledeg > fixangledeg - 0.9 && angledeg < fixangledeg + 0.9) {
        angledeg = fixangledeg;  //  Snap to the nearest multiple of 15 deg
        angle = angledeg / 360.0;
    }
	int mouseStatus = sender->GetMouseStatus();
    
	MoleculeLock(mview->mol);
	MainView_rotateBySlider(mview, angle * 3.1415927 * 2, mode, mouseStatus, MainViewCallback_modifierFlags(NULL));

    //  We do not use here MainViewCallback_drawInfoText(), because the character
    //  "°" here is written in UTF-8 encoding, whereas MainViewCallback_drawInfoText()
    //  accepts a C string in native encoding (shift-JIS in Windows)
    snprintf(buf, sizeof buf, "%.1f°", angledeg);
    wxString labelstr(buf, wxConvUTF8);  //  force UTF-8 encoding here
    infotext->SetLabel(labelstr);

    MainViewCallback_updateCanvas(mview);
	MoleculeUnlock(mview->mol);
}

void
MoleculeView::OnFrameButtonAction(wxMouseEvent &event)
{
	int ival, nframes, bid;
	if (mview == NULL || mview->mol == NULL)
		goto skip;
	nframes = MoleculeGetNumberOfFrames(mview->mol);
	if (nframes == 0)
		goto skip;
	bid = event.GetId();
	if (bid == myID_JumpToStartButton) {
		ival = 0;
	} else if (bid == myID_JumpToEndButton) {
		ival = nframes - 1;
	} else if (bid == myID_PlayForwardButton) {
		ival = (mview->mol->cframe + 1) % nframes;
	} else if (bid == myID_PlayBackwardButton) {
		ival = (mview->mol->cframe + nframes - 1) % nframes;
	}
	//  TODO: implement continuous move
	if (ival >= 0 && ival < nframes) {
		MoleculeLock(mview->mol);
		MoleculeSelectFrame(mview->mol, ival, 1);
		MoleculeUnlock(mview->mol);
		MainViewCallback_display(mview);
		UpdateFrameControlValues();
	}
	
skip:
	event.Skip();
}

void
MoleculeView::OnFrameSliderAction(wxScrollEvent &event)
{
	int ival, nframes;
	if (mview == NULL || mview->mol == NULL)
		return;
	MoleculeLock(mview->mol);
	nframes = MoleculeGetNumberOfFrames(mview->mol);
	if (nframes != 0) {
		ival = frameSlider->GetValue();
		if (ival >= 0 && ival < nframes) {
			MoleculeSelectFrame(mview->mol, ival, 1);
			MoleculeUnlock(mview->mol);
			MainViewCallback_display(mview);
			UpdateFrameControlValues();
			return;
		}
	}
	MoleculeUnlock(mview->mol);
}

void
MoleculeView::OnFrameTextAction(wxCommandEvent &event)
{
	int ival, nframes;
	wxString str;
	if (mview == NULL || mview->mol == NULL)
		return;
	MoleculeLock(mview->mol);
	nframes = MoleculeGetNumberOfFrames(mview->mol);
	if (nframes != 0) {
		str = frameText->GetValue();
		ival = atoi((const char *)str.mb_str(WX_DEFAULT_CONV));
		if (ival >= 0 && ival < nframes) {
			MoleculeSelectFrame(mview->mol, ival, 1);
			MoleculeUnlock(mview->mol);
			MainViewCallback_display(mview);
			UpdateFrameControlValues();
			return;
		}
	}
	MoleculeUnlock(mview->mol);
}

void
MoleculeView::OnStopProgressPressed(wxCommandEvent& event)
{
	mview->mol->requestAbortThread = 1;
}

void
MoleculeView::OnDocumentModified(wxCommandEvent& event)
{
    int xpos, ypos, rowindex, newypos;
	if (!mview->freezeScreen) {
		if (canvas)
			canvas->Refresh();
		UpdateFrameControls();
		MoleculeLock(mview->mol);
        
        //  Get the current row/column value; if the table index is the same, keep the
        //  displayed position as close as possible
        listctrl->GetScrollPosition(&xpos, &ypos);
        if (mview->cachedTableIndex != mview->tableIndex) {
            xpos = ypos = 0;  /*  Different table is shown; should reset the view  */
        } else {
            if (mview->tableIndex >= kMainViewAtomTableIndex && mview->tableIndex <= kMainViewImproperTableIndex) {
                //  Atom, Bond, Angle, Dihedral, Improper; use the cache
                int count = IntGroupGetCount(mview->tableCache);
                if (ypos < count)
                    rowindex = IntGroupGetNthPoint(mview->tableCache, ypos);
                else if (count > 0)
                    rowindex = IntGroupGetNthPoint(mview->tableCache, count - 1);
                else rowindex = 0;
            } else rowindex = ypos;
        }
		MainView_refreshTable(mview);
        if (mview->tableIndex >= kMainViewAtomTableIndex && mview->tableIndex <= kMainViewImproperTableIndex) {
            int ic, is, ie, i;
            ic = IntGroupGetIntervalCount(mview->tableCache);
            newypos = 0;
            for (i = 0; i < ic; i++) {
                is = IntGroupGetStartPoint(mview->tableCache, i);
                ie = IntGroupGetEndPoint(mview->tableCache, i);
                if (rowindex < ie) {
                    if (rowindex >= is)
                        newypos += (rowindex - is);
                    break;
                }
                newypos += ie - is;
            }
        } else newypos = rowindex;
        listctrl->SetScrollPosition(xpos, newypos);
		MoleculeUnlock(mview->mol);
	}
	
	if (mview->tableIndex == kMainViewParameterTableIndex && mview->mol->parameterTableSelectionNeedsClear) {
		/*  Clear parameter selection if necessary  */
		MainViewCallback_setTableSelection(mview, NULL);
		mview->mol->parameterTableSelectionNeedsClear = 0;
	}
	
/*	printf("MoleculeView::OnDocumentModified invoked\n"); */
	event.Skip();  /*  Continue processing of the notification  */
}

void
MoleculeView::OnScriptMenuModified(wxCommandEvent& event)
{
	wxGetApp().UpdateScriptMenu(frame->GetMenuBar());
	event.Skip();
}

void
MoleculeView::OnChar(wxKeyEvent &event)
{
	int code = event.GetKeyCode();
	int mode = 0;
	bool noMod = ((event.GetModifiers() | wxMOD_SHIFT) == wxMOD_SHIFT);
//	MyAppCallback_showScriptMessage("MoleculeView::OnChar invoked\n");
	if (code == WXK_BACK || code == WXK_DELETE || code == 0x7f || code == 8) {
		MoleculeLock(mview->mol);
		MainView_delete(mview);
		MoleculeUnlock(mview->mol);
	} else if (noMod && (code == 'r' || code == 'R'))
		mode = kTrackballRotateMode;
	else if (noMod && (code == 't' || code == 'T'))
		mode = kTrackballTranslateMode;
	else if (noMod && (code == 's' || code == 'S'))
		mode = kTrackballScaleMode;
	else if (noMod && (code == ' '))
		mode = kTrackballSelectionMode;
	else if (noMod && (code == 'b' || code == 'B'))
		mode = kTrackballCreateMode;
	else if (noMod && (code == 'e' || code == 'E'))
		mode = kTrackballEraseMode;
	else {
		event.Skip();
		return;
	}
	if (mode > 0) {
		MainView_setMode(mview, mode);
		MainViewCallback_selectMatrixCellForMode(mview, mode);
	}
}

void
MoleculeView::SelectTable(int idx)
{
	if (idx >= 0 && idx < listmenu->GetCount() && idx != mview->tableIndex) {
		isRebuildingTable = true;
		if (MainView_createColumnsForTableAtIndex(mview, idx) == 0) {
			/*  Invalid menu item  */
			listmenu->SetSelection(mview->tableIndex);
			isRebuildingTable = false;
			return;
		}
		listmenu->SetSelection(idx);
		isRebuildingTable = false;
		MoleculeLock(mview->mol);
		MainView_refreshTable(mview);
		MoleculeUnlock(mview->mol);
		if (idx >= kMainViewParameterTableIndex) {
			MainViewCallback_setTableSelection(mview, NULL);
			if (idx == kMainViewParameterTableIndex) {
				/*  Not sure whether it is appropriate rebuilding MDArena 
				    *every time* the parameter table is opened...  */
				MoleculePrepareMDArena(mview->mol, 1, NULL);
			}
		}
	}
}

void
MoleculeView::OnSelectTable(wxCommandEvent &event)
{
	if (!isRebuildingTable) {
		MoleculeLock(mview->mol);
		SelectTable(listmenu->GetSelection());
		MoleculeUnlock(mview->mol);
	}
}

void
MoleculeView::OnLeftDClickInListCtrl(wxMouseEvent &event)
{
	listctrl->OnLeftDClick(event);
	if (mview->tableIndex >= kMainViewBondTableIndex && mview->tableIndex <= kMainViewImproperTableIndex /* && mview->mol->par != NULL */ ) {
		int row, col, i;
		char indices[64], names[64], types[64], value[20], partypes[64], params[3][20];
        const char *ptype;
        char *parstr;
		wxPoint pos = event.GetPosition();
		if (!listctrl->FindItemAtPosition(pos, &row, &col) || col < 4)
			return;
		/*  Start editing the local parameter; open a separate dialog  */
		MainView_valueForTable(mview, 0, row, indices, sizeof indices);
		MainView_valueForTable(mview, 1, row, names, sizeof names);
		MainView_valueForTable(mview, 2, row, types, sizeof types);
		MainView_valueForTable(mview, 3, row, value, sizeof value);
		MainView_valueForTable(mview, 4, row, partypes, sizeof partypes);
		for (i = 0; i < 3; i++) {
			MainView_valueForTable(mview, 5 + i, row, params[i], sizeof(params[0]));			
		}
		switch (mview->tableIndex) {
			case kMainViewBondTableIndex: ptype = "bond"; break;
			case kMainViewAngleTableIndex: ptype = "angle"; break;
			case kMainViewDihedralTableIndex: ptype = "dihedral"; break;
			case kMainViewImproperTableIndex: ptype = "improper"; break;
			default: return;
		}
		asprintf(&parstr, "%s %s %s", params[0], params[1], params[2]);
		MolActionCreateAndPerform(mview->mol, SCRIPT_ACTION("sssssss"), "cmd_edit_local_parameter_in_mainview", ptype, indices, names, types, value, partypes, parstr);
        free(parstr);
	}
}

void
MoleculeView::OnActivate(wxActivateEvent &event)
{
	if (!event.GetActive()) {
		if (listctrl != NULL)
			listctrl->EndEditText(true);
	}
	event.Skip();
}

void
MoleculeView::OnMoleculeReplaced()
{
	Molecule *mol = ((MyDocument *)GetDocument())->GetMolecule();
	if (mol == NULL && mview == NULL)
		return;
	if (mol != NULL && mol->mview == mview) {
		/*  Clear internal cache  */
		MainView_refreshCachedInfo(mol->mview);
		return;
	}
	if (mview != NULL)
		MainView_setViewObject(mview, NULL);
	if (mol != NULL) {
		int tableIndex;
		mview = mol->mview;
		MainView_setViewObject(mview, this);
		/*  Force updating the table  */
		tableIndex = mview->tableIndex;
		mview->tableIndex = -1;
		SelectTable(tableIndex);
	}
}

void
MoleculeView::EnableProgressIndicator(bool flag)
{
	progress->SetEnabled(flag);
}

#pragma mark ====== MyListCtrl data source ======

int
MoleculeView::GetItemCount(MyListCtrl *ctrl)
{
	return MainView_numberOfRowsInTable(mview);
}

wxString
MoleculeView::GetItemText(MyListCtrl *ctrl, long row, long column) const
{
	char buf[128];
	MainView_valueForTable(mview, column, row, buf, sizeof buf);
	wxString *str = new wxString(buf, WX_DEFAULT_CONV);
	return *str;
}

int
MoleculeView::SetItemText(MyListCtrl *ctrl, long row, long column, const wxString &value)
{
	MainView_setValueForTable(mview, column, row, value.mb_str(WX_DEFAULT_CONV));
	return 0;
}

void
MoleculeView::DragSelectionToRow(MyListCtrl *ctrl, long row)
{
	MainView_dragTableSelectionToRow(mview, row);
}

bool
MoleculeView::IsItemEditable(MyListCtrl *ctrl, long row, long column)
{
	return MainView_isTableItemEditable(mview, column, row);
}

bool
MoleculeView::IsDragAndDropEnabled(MyListCtrl *ctrl, long row)
{
	/*  Only enabled for the atom table  */
	return (MainView_tableType(mview) == 0);
}

void
MoleculeView::OnSelectionChanged(MyListCtrl *ctrl)
{
	MainView_setSelectionFromTable(mview);
}

int
MoleculeView::SetItemColor(MyListCtrl *ctrl, long row, long col, float *fg, float *bg)
{
	if (mview != NULL && mview->mol != NULL) {
		return MainView_setColorForTable(mview, col, row, fg, bg);
	/*
		if (mview->tableIndex == kMainViewParameterTableIndex && col == -1) {
			int src = ParameterTableGetItemSource(mview->mol->par, row);
			if (src == -2) {  // separator line 
				bg[0] = bg[1] = bg[2] = 0.6;
				return 2;
			} else if (src == -1) { //  undefined parameters
				bg[0] = 1.0;
				bg[1] = bg[2] = 0.2;
				return 2;
			} else if (src == 0) {  //  local parameter
				bg[0] = bg[1] = 1.0;
				bg[2] = 0.6;
				return 2;
			}
		} else if (mview->tableIndex > 0 && mview->tableIndex < 5) {
			return MainView_setColorForTable(mview, col, row, fg, bg);
		}
	*/
	}
	return 0;
}

bool
MoleculeView::IsRowSelectable(MyListCtrl *ctrl, long row)
{
    return MainView_isRowSelectable(mview, row);
}

#pragma mark ====== Plain C interface ======

int
MainViewCallback_modifierFlags(void *eventRef)
{
	int flags = 0;
	unsigned modifiers;
	wxMouseState state = ::wxGetMouseState();
	modifiers = 0;
	if (state.ShiftDown())
	  flags |= kShiftKeyMask;
	if (state.AltDown())
	  flags |= kAltKeyMask;
	return flags;
}

int
MainViewCallback_clickCount(void *eventRef)
{
  wxMouseEvent *mevent = (wxMouseEvent *)eventRef;
  if (mevent != NULL) {
    if (mevent->LeftDClick())
      return 2;
    else if (mevent->LeftDown() || mevent->LeftUp())
      return 1;
    else return 0;
  }
  return 0;
}

void
MainViewCallback_lockFocus(MainView *mview)
{
	if (mview != NULL && mview->ref != NULL)
	  ((MoleculeView *)(mview->ref))->canvas->SetCurrent();
}

void
MainViewCallback_unlockFocus(MainView *mview)
{
  //	if (mview != NULL && mview->ref != NULL)
  //  [[(MyWindowController *)(mview->ref) myOpenGLView] unlockFocus];
}

void
MainViewCallback_frame(MainView *mview, float *rect)
{
	if (mview != NULL && mview->ref != NULL) {
	  int width, height;
	  ((MoleculeView *)(mview->ref))->canvas->GetClientSize(&width, &height);
	  rect[0] = rect[1] = 0.0;
	  rect[2] = width;
	  rect[3] = height;
	} else {
		rect[0] = rect[1] = rect[2] = rect[3] = 0.0;
	}
}

float
MainViewCallback_getContentScaleFactor(MainView *mview)
{
    if (gUseGUI && mview != NULL && mview->ref != NULL) {
        return ((MoleculeView *)(mview->ref))->canvas->GetContentScaleFactor();
    } else return 1.0;
}


void
MainViewCallback_display(MainView *mview)
{
	if (gUseGUI && mview != NULL && mview->ref != NULL) {
		wxWindow *canvas = ((MoleculeView *)(mview->ref))->canvas;
		canvas->Refresh();
		canvas->Update();
	}
}

void
MainViewCallback_makeFront(MainView *mview)
{
	if (gUseGUI && mview != NULL && mview->ref != NULL) {
		((MoleculeView *)(mview->ref))->GetFrame()->Raise();
	}
}

void
MainViewCallback_setNeedsDisplay(MainView *mview, int flag)
{
  if (gUseGUI && mview != NULL && mview->ref != NULL) {
    if (flag)
      ((MoleculeView *)(mview->ref))->canvas->Refresh();
  }
}

void
MainViewCallback_updateCanvas(MainView *mview)
{
    if (gUseGUI && mview != NULL && mview->ref != NULL) {
        ((MoleculeView *)(mview->ref))->canvas->Update();
    }
}

void
MainViewCallback_setKeyboardFocus(MainView *mview)
{
	if (gUseGUI && mview != NULL && mview->ref != NULL) {
		((MoleculeView *)(mview->ref))->canvas->SetFocus();
	}
}

void
MainViewCallback_enableToggleButton(MainView *mview, int mode, int flag)
{
	if (mode >= kTrackballRotateMode && mode <= kTrackballEraseMode)
		((MoleculeView *)(mview->ref))->GetToggleButtonAtIndex(mode - kTrackballRotateMode)->Enable(flag);
}

void
MainViewCallback_clearLabels(MainView *mview)
{
	return;
/*
	if (mview != NULL && mview->ref != NULL) {
		id view = [(MyWindowController *)(mview->ref) myOverlayView];
		NSRect bounds = [view bounds];
		[view lockFocus];
		NSEraseRect(bounds);
		[[NSColor cyanColor] set];
		bounds.origin.x = bounds.size.width / 2;
		NSFrameRect(bounds);
		[view unlockFocus];
		[view setNeedsDisplay: YES];
	}
*/
}

void
MainViewCallback_drawLabel(MainView *mview, const float *pos, const char *label)
{
	return;
/*
	if (mview != NULL && mview->ref != NULL) {
		id view = [(MyWindowController *)(mview->ref) myOverlayView];
		NSString *s = [NSString stringWithUTF8String: label];
		NSDictionary *attr = [(MyWindowController *)(mview->ref) labelAttributes];
		[view lockFocus];
		[s drawAtPoint: NSMakePoint(pos[0], pos[1]) withAttributes: attr];
		[view unlockFocus];
		[view setNeedsDisplay: YES];
	}
*/
}

void
MainViewCallback_drawInfoText(MainView *mview, const char *label)
{
	if (gUseGUI && mview != NULL && mview->ref != NULL) {
		wxString labelstr(label, WX_DEFAULT_CONV);
		((MoleculeView *)(mview->ref))->infotext->SetLabel(labelstr);
	}
}

int
MainViewCallback_mouseCheck(MainView *mview)
{
	return 0;
}

void
MainViewCallback_selectMatrixCellForMode(MainView *mview, int mode)
{
  	if (mview != NULL && mview->ref != NULL)
		((MoleculeView *)(mview->ref))->SelectButtonForMode(mode);
}

//int
//MainViewCallback_getTag(MainView *mview)
//{
  //	if (mview != NULL && mview->ref != NULL)
  //	return [(MyWindowController *)(mview->ref) myTag];
  //else return -1;
//}

MainView *
MainViewCallback_viewWithTag(int tag)
{
	wxList &doclist = wxGetApp().DocManager()->GetDocuments();
	wxList::iterator iter;
	int i = 0;
	for (i = 0, iter = doclist.begin(); iter != doclist.end(); ++i, ++iter) {
		if (i == tag) {
			return ((MoleculeView *)(((MyDocument *)(*iter))->GetFirstView()))->mview;
		}
	}
	return NULL;
  //wxList::compatibility_iterator iter = doclist.Item(tag);
  //if (iter != NULL)
  //  return (MyDocument *)(*iter)->GetFirstView()->mview;
  //else
//    return NULL;
  //	int i;
  //	if (sMyWindowControllers == nil)
  //		return NULL;
  //	for (i = [sMyWindowControllers count] - 1; i >= 0; i--) {
  //		id obj = [sMyWindowControllers objectAtIndex: i];
  //	if ([obj myTag] == tag)
  //			return [obj mainView];
  //	}
  //	return NULL;
}

MainView *
MainViewCallback_activeView(void)
{
	if (sActiveViews.GetCount() == 0)
		return NULL;
	else {
		MoleculeView *aview = sActiveViews.Item(0);
		return aview->mview;
	}
//	MoleculeView *cview = (MoleculeView *)(wxGetApp().DocManager()->GetCurrentView());
//	if (cview == NULL)
//		return NULL;
//	else
//		return cview->mview;
}

/*
MainView *
MainViewCallback_newFromFile(const char *fname)
{
  wxDocument *doc;
  wxDocManager *manager = wxGetApp().DocManager();
  if (fname == NULL || *fname == 0) {
    doc = manager->CreateDocument(wxT(""), wxDOC_NEW);
  } else {
    wxString fnamestr(fname, wxConvFile);
    doc = manager->CreateDocument(fnamestr, wxDOC_SILENT);
  }
  return MainViewCallback_activeView();
}
*/

int
MainViewCallback_importFromFile(MainView *mview, const char *fname)
{
  MyDocument *doc;
  if (mview != NULL && mview->ref != NULL && (doc = (((MoleculeView *)(mview->ref))->MolDocument())) != NULL) {
    wxString fnamestr(fname, wxConvFile);
    // doc->importFromFile(fnamestr);
    MainViewCallback_setNeedsDisplay(mview, 1);
    return 1;
  }
  return 0;
}

void
MainViewCallback_getFilename(MainView *mview, char *buf, int bufsize)
{
  MyDocument *doc;
  if (mview != NULL && mview->ref != NULL && (doc = (((MoleculeView *)(mview->ref))->MolDocument())) != NULL) {
    wxString fname;
    fname = doc->GetFilename();
    strncpy(buf, (const char*)fname.mb_str(wxConvFile), bufsize - 1);
    buf[bufsize - 1] = 0;
  } else {
    buf[0] = 0;
  }
}

/*
void
MainViewCallback_moleculeReplaced(MainView *mview, struct Molecule *mol)
{
	if (mview != NULL && mview->ref != NULL) {
		MyDocument *doc = ((MoleculeView *)(mview->ref))->MolDocument();
		if (doc != NULL)
			doc->SetMolecule(mol);
		MyListCtrl *listctrl = ((MoleculeView *)(mview->ref))->GetListCtrl();
		if (listctrl != NULL)
			listctrl->SetDataSource((MoleculeView *)(mview->ref));
	}
}
*/

typedef struct Label {
  //	StringTexture *tex;
} Label;

struct Label *
MainViewCallback_newLabel(MainView *mview, const char *message, float fontsize, const float *forecolor, const float *backcolor)
{
  /*
	Label *label;
	NSDictionary *attr;
	NSColor *textColor, *boxColor;
	label = (Label *)malloc(sizeof(Label));
	if (label == NULL)
		return NULL;
	memset(label, 0, sizeof(Label));
//	MainViewCallback_lockFocus(mview);
	if (forecolor != NULL)
		textColor = [NSColor colorWithDeviceRed: forecolor[0] green: forecolor[1] blue: forecolor[2] alpha: forecolor[3]];
	else
		textColor = [NSColor whiteColor];
	if (backcolor != NULL)
		boxColor = [NSColor colorWithDeviceRed: backcolor[0] green: backcolor[1] blue: backcolor[2] alpha: backcolor[3]];
	else
		boxColor = [NSColor colorWithDeviceRed:1.0f green:1.0f blue:1.0f alpha:0.0f];
	attr = [NSDictionary dictionaryWithObjectsAndKeys: [NSFont userFontOfSize: fontsize], NSFontAttributeName, textColor, NSForegroundColorAttributeName, nil];
	label->tex = [[StringTexture alloc] 
		initWithString: [NSString stringWithUTF8String: message] 
		withAttributes: attr
		withTextColor: textColor
		withBoxColor: boxColor
		withBorderColor: [NSColor colorWithDeviceRed:1.0f green:1.0f blue:1.0f alpha:0.0f]];
//	MainViewCallback_unlockFocus(mview);
	return label;
  */
	return NULL;
}

void
MainViewCallback_releaseLabel(struct Label *label)
{
  /*	if (label != NULL) {
		[label->tex release];
		free(label);
	}
  */
}

void
MainViewCallback_drawLabelAtPoint(struct Label *label, const float *pos)
{
  /*
	if (label != NULL && pos != NULL) {
	//	GLint matrixMode;
		glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	//	glGetIntegerv (GL_MATRIX_MODE, &matrixMode);
	//	glMatrixMode(GL_MODELVIEW);
	//	glPushMatrix();
	//	glTranslatef(0.0f, 0.0f, pos[3]);
		[label->tex drawAtPoint: NSMakePoint(pos[0], pos[1]) withDepth: pos[2]];
	//	glPopMatrix();
	//	glMatrixMode (matrixMode);
	}
  */
}

void
MainViewCallback_labelSize(struct Label *label, float *outSize)
{
  /*	if (label != NULL) {
		NSSize size = [label->tex frameSize];
		if (outSize != NULL) {
			outSize[0] = size.width;
			outSize[1] = size.height;
		}
	}
  */
}

int
MainViewCallback_exportGraphic(MainView *mview, const char *fname, float scale, int bg_color, int width, int height)
{
    if (!gUseGUI)
        return 0;
	if (mview != NULL && mview->ref != NULL && ((MoleculeView *)(mview->ref))->MolDocument() != NULL) {
		wxString fnamestr(fname, wxConvFile);
		return ((MoleculeView *)(mview->ref))->DoExportGraphic(fnamestr, scale, bg_color, width, height);
	}
	return -100;
}

#pragma mark ====== Plain C Interface (MyListCtrl) ======

/*  These interface functions are also used for accessing MyListCtrl of GlobalParameterFrame (when mview == NULL) */

static MyListCtrl *
s_MyListCtrlFromMainView(MainView *mview)
{
	if (mview == NULL)
		return wxGetApp().GetGlobalParameterListCtrl();
	if (mview != NULL && mview->ref != NULL)
		return ((MoleculeView *)(mview->ref))->GetListCtrl();
	else return NULL;
}

void
MainViewCallback_selectTable(MainView *mview, int idx)
{
	if (mview != NULL && mview->ref != NULL)
		((MoleculeView *)(mview->ref))->SelectTable(idx);
}

int
MainViewCallback_numberOfTableColumns(MainView *mview)
{
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl != NULL)
		return listctrl->GetColumnCount();
	else return 0;
}

int
MainViewCallback_addTableColumn(MainView *mview, const char *name, int width, int editable)
{
	int idx;
	wxString nstr(name, WX_DEFAULT_CONV);
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl == NULL)
		return 0;
	idx = listctrl->GetColumnCount();
	listctrl->InsertColumn(idx, nstr, wxLIST_FORMAT_LEFT);
	listctrl->SetColumnWidth(idx, FromFrameDIP(listctrl, width * 10));
	return idx;
}

int
MainViewCallback_removeTableColumnAtIndex(MainView *mview, int idx)
{
	int ncolumns;
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl == NULL)
		return 0;
	ncolumns = listctrl->GetColumnCount();
	if (idx >= 0 && idx < ncolumns) {
		listctrl->DeleteColumn(idx);
		ncolumns--;
	}
	return ncolumns;
}

void
MainViewCallback_reloadTableData(MainView *mview)
{
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl != NULL)
		listctrl->RefreshTable();
/*	
	int nrows = MainView_numberOfRowsInTable(mview);
	listctrl->SetItemCount(nrows);
	if (nrows > 0)
		listctrl->RefreshItems(0, nrows - 1);
*/
}

void
MainViewCallback_setTableSelection(MainView *mview, IntGroup *selection)
{
	int i, n, n1, n2;
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl == NULL)
		return;
	n = 0;
	listctrl->EnableSelectionChangeNotification(false);
	if (selection != NULL) {
		for (i = 0; (n1 = IntGroupGetStartPoint(selection, i)) >= 0; i++) {
			n2 = IntGroupGetEndPoint(selection, i);
			while (n < n1) {
				listctrl->SetItemState(n, 0, wxLIST_STATE_SELECTED);
				n++;
			}
			while (n < n2) {
				listctrl->SetItemState(n, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED);
				n++;
			}
		}
	}
	listctrl->RefreshTable();
	n1 = MainView_numberOfRowsInTable(mview);
	while (n < n1) {
		listctrl->SetItemState(n, 0, wxLIST_STATE_SELECTED);
		n++;
	}
//	listctrl->RefreshItems(0, n1 - 1);
	{
		//  EVT_LIST_ITEM_SELECTED is sent by wxPostEvent rather than ProcessEvent,
		//  so re-enable of table selection event should also be sent by wxPostEvent.
		//  Otherwise, the enable flag is set to true _before_ EVT_LIST_ITEM_SELECTED is sent.
		wxCommandEvent myEvent(MyListCtrlEvent, MyListCtrlEvent_enableTableSelectionNotification);
		wxPostEvent(listctrl, myEvent);
	}
}

IntGroup *
MainViewCallback_getTableSelection(MainView *mview)
{
	int i, n;
	IntGroup *ig;
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl == NULL)
		return NULL;
	ig = IntGroupNew();
	n = MainView_numberOfRowsInTable(mview);
	for (i = 0; i < n; i++) {
		if (listctrl->GetItemState(i, wxLIST_STATE_SELECTED) != 0)
			IntGroupAdd(ig, i, 1);
	}
	return ig;
}

void
MainViewCallback_showTable(MainView *mview)
{
  //	if (mview != NULL && mview->ref != NULL) {
  //	MyDocument *doc = (MyDocument *)[(id)(mview->ref) document];
  //	[doc showTable: doc];
  //}
}

void
MainViewCallback_hideTable(MainView *mview)
{
  //	if (mview != NULL && mview->ref != NULL) {
  //	MyDocument *doc = (MyDocument *)[(id)(mview->ref) document];
  //		[doc hideTable: doc];
  //	}
}

void
MainViewCallback_ensureVisible(MainView *mview, int row)
{
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl != NULL)
		listctrl->EnsureVisible(row);
}

void
MainViewCallback_startEditText(MainView *mview, int row, int column)
{
	MyListCtrl *listctrl = s_MyListCtrlFromMainView(mview);
	if (listctrl != NULL)
		listctrl->StartEditText(row, column);
}

