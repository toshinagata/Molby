/*
 *  MyClipboardData.h
 *  Molby
 *
 *  Created by Toshi Nagata on 08/11/26.
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

#ifndef __MyClipboardData_h__
#define __MyClipboardData_h__

#include "wx/dataobj.h"

class MyClipboardData: public wxDataObjectSimple {
public:
	char *buffer;
	int length;
	
	MyClipboardData(const char *type);
	virtual ~MyClipboardData();

	virtual size_t GetDataSize() const;
	virtual bool GetDataHere(void *buf) const;
	virtual bool SetData(size_t len, const void *buf);
};

#endif /* __MyClipboardData_h__ */
