/*
 *  MDGraphite.h
 *  Molby
 *
 *  Created by Toshi Nagata on 07/10/07.
 *  Copyright 2007 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __MDGRAPHITE_H__
#define __MDGRAPHITE_H__

#include "MDCore.h"

#ifdef __cplusplus
extern "C" {
#endif
	
typedef struct MDGraphiteArena MDGraphiteArena;

MDGraphiteArena *graphite_new(void);
void graphite_release(MDGraphiteArena *graphite);
void graphite_set_origin(MDGraphiteArena *graphite, const Vector *vp);
void graphite_force(MDGraphiteArena *graphite, MDArena *arena, Double *energy, Vector *force);
void graphite_get_axes(MDGraphiteArena *graphite, Vector *op, Vector *xp, Vector *yp, Vector *zp, Double *rp);
void graphite_set_needs_update(MDGraphiteArena *graphite, int flag);
	
#ifdef __cplusplus
}
#endif
		
#endif /* __MDGRAPHITE_H__ */
