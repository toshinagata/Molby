/*
 *  MyIPCSupport.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 12/10/10.
 *  Copyright 2012 Toshi Nagata. All rights reserved.
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

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif


#if defined(__WXMSW__)

#include "MyIPCSupport.h"
#include "MyApp.h"

wxString *gIPCServiceName = NULL;

bool
MyClientConnection::OnDisconnect()
{
	wxGetApp().m_client->Disconnect();
}

MyClient::MyClient()
{
	m_clientConnection = NULL;
}

MyClient::~MyClient()
{
	Disconnect();
}
	
void
MyClient::Disconnect()
{
	if (m_clientConnection != NULL) {
		m_clientConnection->Disconnect();
		m_clientConnection = NULL;
	}
}

wxConnectionBase *
MyClient::OnMakeConnection()
{
	if (m_clientConnection == NULL)
		m_clientConnection = new MyClientConnection;
	return m_clientConnection;
}

bool
MyServerConnection::OnDisconnect()
{
	wxGetApp().m_server->Disconnect();
}

bool
MyServerConnection::OnExecute(const wxString& topic, wxChar* data, int size, wxIPCFormat format)
{
	if (topic == MOLBY_IPC_TOPIC) {
		wxString files(data);
		wxGetApp().RequestOpenFilesByIPC(files);
		return true;
	} else return false;
}

MyServer::MyServer()
{
	m_serverConnection = NULL;
}

MyServer::~MyServer()
{
	Disconnect();
}

void
MyServer::Disconnect()
{
	if (m_serverConnection != NULL) {
		m_serverConnection->Disconnect();
		m_serverConnection = NULL;
	}
}

wxConnectionBase *
MyServer::OnAcceptConnection(const wxString &topic)
{
    if (topic == MOLBY_IPC_TOPIC) {
		if (m_serverConnection == NULL)
			m_serverConnection = new MyServerConnection();
        return m_serverConnection;
    }
    return NULL;
}

#endif  // defined(__WXMSW__)