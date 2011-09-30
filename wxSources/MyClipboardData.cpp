/*
 *  MyClipboardData.cpp
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

#include "MyClipboardData.h"
#include "MyMBConv.h"

MyClipboardData::MyClipboardData(const char *type):
	wxDataObjectSimple()
{
	wxString ftype(type, WX_DEFAULT_CONV);
	wxDataFormat customFormat((const wxChar *)ftype);
	buffer = NULL;
	length = 0;
	SetFormat(customFormat);
}

MyClipboardData::~MyClipboardData()
{
	if (buffer != NULL)
		free(buffer);
}

size_t
MyClipboardData::GetDataSize() const
{
	return length;
}

bool
MyClipboardData::GetDataHere(void *buf) const
{
	if (buffer != NULL) {
		memmove(buf, buffer, length);
		return true;
	} else return false;
}

bool
MyClipboardData::SetData(size_t len, const void *buf)
{
	char *p;
	if (buffer == NULL)
		p = (char *)malloc(len);
	else
		p = (char *)realloc(buffer, len);
	if (p == NULL)
		return false;
	memmove(p, buf, len);
	buffer = p;
	length = len;
	return true;
}
