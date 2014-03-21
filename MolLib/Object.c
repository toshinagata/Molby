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
	if (obj->name != NULL)
		free((void *)(obj->name));
	obj->name = strdup(name);
	return;
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

