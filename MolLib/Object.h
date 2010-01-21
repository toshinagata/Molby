/*
 *  Object.h
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

#ifndef __Object_h__
#define __Object_h__

#ifdef __cplusplus
extern "C" {
#endif
	
/*  A "base object type" to implement reference counting
    and management by unique names  */
typedef struct Object {
	int    refCount;
	struct Object *next;
	const char *name;
} Object;

void ObjectInit(Object *obj, Object **objRootPtr, const char *name);
int ObjectIncrRefCount(Object *obj);
int ObjectDecrRefCount(Object *obj);
void ObjectDealloc(Object *obj, Object **objRootPtr);
void ObjectSetName(Object *obj, const char *name, Object *objRoot);
const char *ObjectGetName(Object *obj);
Object *ObjectWithName(const char *name, Object *objRoot);

#ifdef __cplusplus
}
#endif
		
#endif /* __Object_h__ */
