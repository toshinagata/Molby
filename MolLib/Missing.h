/*
 *  Missing.h
 *  Molby
 *
 *  Created by Toshi Nagata on 08/11/06.
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

#ifndef __MISSING_H__
#define __MISSING_H__

#ifdef __WXMSW__
#define MISSING_STRSEP 1
#define MISSING_ASPRINTF 1
#define MISSING_STRTOK_R 1
#define MISSING_MERGESORT 1
#endif

#ifdef __cplusplus
extern "C" {
#endif

void translate_char(char *p, int from, int to);

#if MISSING_STRSEP
char *strpbrk(const char *cs, const char *ct);
char *strsep(char **stringp, const char *delim);
#endif

#if MISSING_ASPRINTF
#include <stdarg.h>
int asprintf(char **ret, const char *format, ...);
int vasprintf(char **ret, const char *format, va_list ap);
#endif /* MISSING_ASPRINTF */

#if MISSING_STRTOK_R
char *strtok_r(char *str, const char *sep, char **lasts);
#endif /* MISSING_STRTOK_R */

#if MISSING_MERGESORT
#include <stdlib.h>
int mergesort(void *base, size_t nel, size_t width, int (*compar)(const void *, const void *));
#endif /* MISSING_MERGESORT */

#ifdef __cplusplus
}
#endif

#endif /* __MISSING_H__ */
