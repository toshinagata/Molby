/*
 *  Object.c
 *
 *  Created by Toshi Nagata on 06/03/11.
 *  Copyright 2006-2008 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "Object.h"
#include "Types.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

void
ObjectInit(Object *obj, Object **objRootPtr, const char *name)
{
	obj->refCount = 1;
	obj->next = *objRootPtr;
	*objRootPtr = obj;
	ObjectSetName(obj, name, *objRootPtr);
}

int
ObjectIncrRefCount(Object *obj)
{
	if (obj != NULL)
		return ++(obj->refCount);
	else return -1;
}

int
ObjectDecrRefCount(Object *obj)
{
	if (obj != NULL)
		return --(obj->refCount);
	else return -1;
}

void
ObjectDealloc(Object *obj, Object **objRootPtr)
{
	Object **objp;
	free((void *)(obj->name));
	for (objp = objRootPtr; *objp != NULL; objp = &((*objp)->next)) {
		if (*objp == obj) {
			*objp = obj->next;
			break;
		}
	}
	free(obj);
}

void
ObjectSetName(Object *obj, const char *name, Object *objRoot)
{
	obj->name = strdup(name);
	return;
#if 0
	char *buf, *p;
	int n;
	size_t size;
	/* Allocate a copy of string with some space to add a number to the end */
	size = strlen(name) + 6;
	buf = (char *)malloc(size);
	if (buf == NULL)
		Panic("Cannot set object name");
	strcpy(buf, name);
	/* Determine the position to place the number; if name ends with space + number,
	   then the number will be replaced. Otherwise, the number is appended at the
	   end of the string with one space */
	n = 0;
	p = buf + strlen(name) - 1;
	while (p >= buf && *p == ' ')
		p--;
	if (p >= buf && isdigit(*p)) {
		char *p1 = p;
		while (p1 >= buf && isdigit(*p1))
			p1--;
		if (p1 >= buf && *p1 == ' ') {
			n = atoi(p1 + 1);
			p = p1;
		}
	}
	p++;
	while (n < 1000) {
		Object *op;
		for (op = objRoot; op != NULL; op = op->next) {
			if (op == obj)
				continue;
			if (strcmp(op->name, buf) == 0)
				break;
		}
		if (op == NULL)
			break;
		if (p > buf && p[-1] != ' ')
			*p++ = ' ';
		snprintf(p, size - (p - buf), "%d", ++n);
	}
	if (n >= 1000)
		Panic("Cannot set object name");
	if (obj->name != NULL)
		free((void *)(obj->name));
	obj->name = buf;
#endif
}

const char *
ObjectGetName(Object *obj)
{
	if (obj == NULL)
		return NULL;
	else return obj->name;
}

Object *
ObjectWithName(const char *name, Object *objRoot)
{
	Object *obj;
	for (obj = objRoot; obj != NULL; obj = obj->next) {
		if (strcmp(obj->name, name) == 0)
			return obj;
	}
	return NULL;
}

