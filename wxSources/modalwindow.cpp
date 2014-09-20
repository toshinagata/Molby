/**********************************************************************************************
 *
 * Filename  : modalwindow.cpp
 * Purpose   : Allow a modalwindow like wxDialog but allowing menus and such.
 * Author    : John A. Mason
 * Created   : 8/27/2008 07:54:12 AM
 * Copyright : Released under wxWidgets original license.
 *
**********************************************************************************************/

#include "modalwindow.h"

IMPLEMENT_DYNAMIC_CLASS(wxModalWindow, wxFrame)

wxModalWindow::wxModalWindow()
: wxFrame() {
   Init();
}

wxModalWindow::wxModalWindow(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style, const wxString& name)
: wxFrame(parent, id, title, pos, size, style, name) {
   Init();
}


wxModalWindow::~wxModalWindow() {
   delete m_eventLoop;
}


void wxModalWindow::Init()
{
    m_returnCode = 0;
    m_windowDisabler = NULL;
    m_eventLoop = NULL;
    m_isShowingModal = false;
}

bool wxModalWindow::Show(bool show)
{
    if ( !show )
    {
        // if we had disabled other app windows, reenable them back now because
        // if they stay disabled Windows will activate another window (one
        // which is enabled, anyhow) and we will lose activation
        if ( m_windowDisabler )
        {
            delete m_windowDisabler;
            m_windowDisabler = NULL;
        }

        if ( IsModal() )
            EndModal(wxID_CANCEL);
    }

    bool ret = wxFrame::Show(show);

   //I don't think we need this. Since it is a wxFrame that we are extending,
   // we don't need wxEVT_INIT_DIALOG firing off - that's what InitDialog does...
   // and this would only make sense if we have a wxDialog and validators
//    if ( show )
        //InitDialog();

    return ret;
}

bool wxModalWindow::IsModal() const {
   return m_isShowingModal;
}

int wxModalWindow::ShowModal() {
  if ( IsModal() )
    {
       wxFAIL_MSG( wxT("wxModalWindow:ShowModal called twice") );
       return GetReturnCode();
    }

    // use the apps top level window as parent if none given unless explicitly
    // forbidden
    if ( !GetParent() )
    {
        wxWindow *parent = wxTheApp->GetTopWindow();
        if ( parent && parent != this )
        {
            m_parent = parent;
        }
    }

    Show(true);

    m_isShowingModal = true;

    wxASSERT_MSG( !m_windowDisabler, _T("disabling windows twice?") );

#if defined(__WXGTK__) || defined(__WXMGL__)
    wxBusyCursorSuspender suspender;
    // FIXME (FIXME_MGL) - make sure busy cursor disappears under MSW too
#endif

/*  20140919 Toshi Nagata  */
#if defined(__WXMAC__)
	/*  --- OSX specific --- */
	{
		extern void MacRunModalForWindow(void *w);  /*  Defined in modalwindow_osx.mm  */
		WXWindow w = this->GetWXWindow();
		MacRunModalForWindow(w);
	}
	
#else
	/*  --- The original code ---  */

	m_windowDisabler = new wxWindowDisabler(this);
    if ( !m_eventLoop )
        m_eventLoop = new wxEventLoop;

    m_eventLoop->Run();

#endif
/*  End Toshi Nagata  */

    return GetReturnCode();
}



void wxModalWindow::EndModal(int retCode) {

	/*  20140919 Toshi Nagata  */

#if defined(__WXMAC__)
	/*  --- OSX specific ---  */
	{
		extern void MacStopModal(void);  /*  Defined in modalwindow_osx.m  */
		SetReturnCode(retCode);
		m_isShowingModal = false;
		MacStopModal();
	}
#else
    wxASSERT_MSG( m_eventLoop, _T("wxModalWindow is not modal") );

    SetReturnCode(retCode);

    if ( !IsModal() )
    {
        wxFAIL_MSG( wxT("wxModalWindow:EndModal called twice") );
        return;
    }

    m_isShowingModal = false;

	m_eventLoop->Exit();
#endif
/*  End Toshi Nagata  */
	
    Show(false);
}


void wxModalWindow::SetReturnCode(int retCode) {
   m_returnCode=retCode;
}


int wxModalWindow::GetReturnCode() const {
   return m_returnCode;
}

