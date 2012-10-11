/*
 *  MyIPCSupport.h
 *  Molby
 *
 *  Created by Toshi Nagata on 12/10/10.
 *  Copyright 2010 Toshi Nagata. All rights reserved.
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

#if defined(__WXMSW__)

#include "wx/ipc.h"

#define MOLBY_IPC_TOPIC wxT("MOLBY_IPC_TOPIC")

extern wxString *gIPCServiceName;

class MyClientConnection: public wxConnection
{
public:
	virtual bool OnDisconnect();
};

class MyClient: public wxClient
{
public:
	MyClient();
	~MyClient();
	void Disconnect();
    wxConnectionBase *OnMakeConnection();
	MyClientConnection *m_clientConnection;
};

class MyServerConnection: public wxConnection
{
public:
	virtual bool OnDisconnect();
	virtual bool OnExecute(const wxString& topic, wxChar* data, int size, wxIPCFormat format);
};

class MyServer: public wxServer
{
public:
	MyServer();
	~MyServer();
	void Disconnect();
    wxConnectionBase *OnAcceptConnection(const wxString& topic);
	MyServerConnection *m_serverConnection;
};

#endif // defined(__WXMSW__)
