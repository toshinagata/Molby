/*
 *  OpenGL_extensions.c
 *  Molby
 *
 *  Created by Toshi Nagata on 2014/03/17.
 *  Copyright 2008-2014 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

/*  Implement OpenGL extentions  */

#include "OpenGL_extensions.h"

/*  Function pointers  */
typedef GLboolean (*t_glIsRenderbufferEXT)(GLuint renderbuffer);
typedef void (*t_glBindRenderbufferEXT)(GLenum target, GLuint renderbuffer);
typedef void (*t_glDeleteRenderbuffersEXT)(GLsizei n, const GLuint * renderbuffers);
typedef void (*t_glGenRenderbuffersEXT)(GLsizei n, GLuint * renderbuffers);
typedef void (*t_glRenderbufferStorageEXT)(GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
typedef void (*t_glGetRenderbufferParameterivEXT)(GLenum target, GLenum pname, GLint * params);
typedef GLboolean (*t_glIsFramebufferEXT)(GLuint framebuffer);
typedef void (*t_glBindFramebufferEXT)(GLenum target, GLuint framebuffer);
typedef void (*t_glDeleteFramebuffersEXT)(GLsizei n, const GLuint * framebuffers);
typedef void (*t_glGenFramebuffersEXT)(GLsizei n, GLuint * framebuffers);
typedef GLenum (*t_glCheckFramebufferStatusEXT)(GLenum target);
typedef void (*t_glFramebufferTexture1DEXT)(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
typedef void (*t_glFramebufferTexture2DEXT)(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
typedef void (*t_glFramebufferTexture3DEXT)(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
typedef void (*t_glFramebufferRenderbufferEXT)(GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
typedef void (*t_glGetFramebufferAttachmentParameterivEXT)(GLenum target, GLenum attachment, GLenum pname, GLint * params);
typedef void (*t_glGenerateMipmapEXT)(GLenum target);

static t_glIsRenderbufferEXT p_glIsRenderbufferEXT;
static t_glBindRenderbufferEXT p_glBindRenderbufferEXT;
static t_glDeleteRenderbuffersEXT p_glDeleteRenderbuffersEXT;
static t_glGenRenderbuffersEXT p_glGenRenderbuffersEXT;
static t_glRenderbufferStorageEXT p_glRenderbufferStorageEXT;
static t_glGetRenderbufferParameterivEXT p_glGetRenderbufferParameterivEXT;
static t_glIsFramebufferEXT p_glIsFramebufferEXT;
static t_glBindFramebufferEXT p_glBindFramebufferEXT;
static t_glDeleteFramebuffersEXT p_glDeleteFramebuffersEXT;
static t_glGenFramebuffersEXT p_glGenFramebuffersEXT;
static t_glCheckFramebufferStatusEXT p_glCheckFramebufferStatusEXT;
static t_glFramebufferTexture1DEXT p_glFramebufferTexture1DEXT;
static t_glFramebufferTexture2DEXT p_glFramebufferTexture2DEXT;
static t_glFramebufferTexture3DEXT p_glFramebufferTexture3DEXT;
static t_glFramebufferRenderbufferEXT p_glFramebufferRenderbufferEXT;
static t_glGetFramebufferAttachmentParameterivEXT p_glGetFramebufferAttachmentParameterivEXT;
static t_glGenerateMipmapEXT p_glGenerateMipmapEXT;

static int p_inited = 0;

int
InitializeOpenGLExtensions(void)
{
    int n = 0;
    if (p_inited != 0)
        return 0;
    p_glIsRenderbufferEXT = (t_glIsRenderbufferEXT)wglGetProcAddress("glIsRenderbufferEXT");
    if (p_glIsRenderbufferEXT == NULL)
        n++;
    p_glBindRenderbufferEXT = (t_glBindRenderbufferEXT)wglGetProcAddress("glBindRenderbufferEXT");
    if (p_glBindRenderbufferEXT == NULL)
        n++;
    p_glDeleteRenderbuffersEXT = (t_glDeleteRenderbuffersEXT)wglGetProcAddress("glDeleteRenderbuffersEXT");
    if (p_glDeleteRenderbuffersEXT == NULL)
        n++;
    p_glGenRenderbuffersEXT = (t_glGenRenderbuffersEXT)wglGetProcAddress("glGenRenderbuffersEXT");
    if (p_glGenRenderbuffersEXT == NULL)
        n++;
    p_glRenderbufferStorageEXT = (t_glRenderbufferStorageEXT)wglGetProcAddress("glRenderbufferStorageEXT");
    if (p_glRenderbufferStorageEXT == NULL)
        n++;
    p_glGetRenderbufferParameterivEXT = (t_glGetRenderbufferParameterivEXT)wglGetProcAddress("glGetRenderbufferParameterivEXT");
    if (p_glGetRenderbufferParameterivEXT == NULL)
        n++;
    p_glIsFramebufferEXT = (t_glIsFramebufferEXT)wglGetProcAddress("glIsFramebufferEXT");
    if (p_glIsFramebufferEXT == NULL)
        n++;
    p_glBindFramebufferEXT = (t_glBindFramebufferEXT)wglGetProcAddress("glBindFramebufferEXT");
    if (p_glBindFramebufferEXT == NULL)
        n++;
    p_glDeleteFramebuffersEXT = (t_glDeleteFramebuffersEXT)wglGetProcAddress("glDeleteFramebuffersEXT");
    if (p_glDeleteFramebuffersEXT == NULL)
        n++;
    p_glGenFramebuffersEXT = (t_glGenFramebuffersEXT)wglGetProcAddress("glGenFramebuffersEXT");
    if (p_glGenFramebuffersEXT == NULL)
        n++;
    p_glCheckFramebufferStatusEXT = (t_glCheckFramebufferStatusEXT)wglGetProcAddress("glCheckFramebufferStatusEXT");
    if (p_glCheckFramebufferStatusEXT == NULL)
        n++;
    p_glFramebufferTexture1DEXT = (t_glFramebufferTexture1DEXT)wglGetProcAddress("glFramebufferTexture1DEXT");
    if (p_glFramebufferTexture1DEXT == NULL)
        n++;
    p_glFramebufferTexture2DEXT = (t_glFramebufferTexture2DEXT)wglGetProcAddress("glFramebufferTexture2DEXT");
    if (p_glFramebufferTexture2DEXT == NULL)
        n++;
    p_glFramebufferTexture3DEXT = (t_glFramebufferTexture3DEXT)wglGetProcAddress("glFramebufferTexture3DEXT");
    if (p_glFramebufferTexture3DEXT == NULL)
        n++;
    p_glFramebufferRenderbufferEXT = (t_glFramebufferRenderbufferEXT)wglGetProcAddress("glFramebufferRenderbufferEXT");
    if (p_glFramebufferRenderbufferEXT == NULL)
        n++;
    p_glGetFramebufferAttachmentParameterivEXT = (t_glGetFramebufferAttachmentParameterivEXT)wglGetProcAddress("glGetFramebufferAttachmentParameterivEXT");
    if (p_glGetFramebufferAttachmentParameterivEXT == NULL)
        n++;
    p_glGenerateMipmapEXT = (t_glGenerateMipmapEXT)wglGetProcAddress("glGenerateMipmapEXT");
    if (p_glGenerateMipmapEXT == NULL)
        n++;
    p_inited = 1;
    return n;
}

GLboolean glIsRenderbufferEXT(GLuint renderbuffer)
{
	return (*p_glIsRenderbufferEXT)(renderbuffer);
}

void glBindRenderbufferEXT(GLenum target, GLuint renderbuffer)
{
	(*p_glBindRenderbufferEXT)(target, renderbuffer);
}

void glDeleteRenderbuffersEXT(GLsizei n, const GLuint * renderbuffers)
{
	(*p_glDeleteRenderbuffersEXT)(n, renderbuffers);
}

void glGenRenderbuffersEXT(GLsizei n, GLuint * renderbuffers)
{
	(*p_glGenRenderbuffersEXT)(n, renderbuffers);
}

void glRenderbufferStorageEXT(GLenum target, GLenum internalformat, GLsizei width, GLsizei height)
{
	(*p_glRenderbufferStorageEXT)(target, internalformat, width, height);
}

void glGetRenderbufferParameterivEXT(GLenum target, GLenum pname, GLint * params)
{
	(*p_glGetRenderbufferParameterivEXT)(target, pname, params);
}

GLboolean glIsFramebufferEXT(GLuint framebuffer)
{
	(*p_glIsFramebufferEXT)(framebuffer);
}

void glBindFramebufferEXT(GLenum target, GLuint framebuffer)
{
	(*p_glBindFramebufferEXT)(target, framebuffer);
}

void glDeleteFramebuffersEXT(GLsizei n, const GLuint * framebuffers)
{
	(*p_glDeleteFramebuffersEXT)(n, framebuffers);
}

void glGenFramebuffersEXT(GLsizei n, GLuint * framebuffers)
{
	(*p_glGenFramebuffersEXT)(n, framebuffers);
}

GLenum glCheckFramebufferStatusEXT(GLenum target)
{
	return (*p_glCheckFramebufferStatusEXT)(target);
}

void glFramebufferTexture1DEXT(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level)
{
	(*p_glFramebufferTexture1DEXT)(target, attachment, textarget, texture, level);
}

void glFramebufferTexture2DEXT(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level)
{
	(*p_glFramebufferTexture2DEXT)(target, attachment, textarget, texture, level);
}

void glFramebufferTexture3DEXT(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset)
{
	(*p_glFramebufferTexture3DEXT)(target, attachment, textarget, texture, level, zoffset);
}

void glFramebufferRenderbufferEXT(GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer)
{
	(*p_glFramebufferRenderbufferEXT)(target, attachment, renderbuffertarget, renderbuffer);
}

void glGetFramebufferAttachmentParameterivEXT(GLenum target, GLenum attachment, GLenum pname, GLint * params)
{
	(*p_glGetFramebufferAttachmentParameterivEXT)(target, attachment, pname, params);
}

void glGenerateMipmapEXT(GLenum target)
{
	(*p_glGenerateMipmapEXT)(target);
}
