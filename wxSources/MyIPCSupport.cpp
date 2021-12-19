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
	return wxGetApp().m_client->Disconnect();
}

MyClient::MyClient()
{
	m_clientConnection = NULL;
}

MyClient::~MyClient()
{
	Disconnect();
}
	
bool
MyClient::Disconnect()
{
	if (m_clientConnection != NULL) {
        if (m_clientConnection->Disconnect()) {
            m_clientConnection = NULL;
            return true;
        } else return false;
    } else return true;
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
	return wxGetApp().m_server->Disconnect();
}

bool
MyServerConnection::OnExecute(const wxString& topic, const void *data, size_t size, wxIPCFormat format)
{
	if (topic == MOLBY_IPC_TOPIC) {
		wxString files((wxChar *)data);
		wxGetApp().RequestOpenFilesByEvent(files);
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

bool
MyServer::Disconnect()
{
	if (m_serverConnection != NULL) {
        if (m_serverConnection->Disconnect()) {
            m_serverConnection = NULL;
            return true;
        } else return false;
    } else return true;
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
