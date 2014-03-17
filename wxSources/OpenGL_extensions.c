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
static GLboolean (*p_glIsRenderbufferEXT)(GLuint renderbuffer);
static void (*p_glBindRenderbufferEXT)(GLenum target, GLuint renderbuffer);
static void (*p_glDeleteRenderbuffersEXT)(GLsizei n, const GLuint * renderbuffers);
static void (*p_glGenRenderbuffersEXT)(GLsizei n, GLuint * renderbuffers);
static void (*p_glRenderbufferStorageEXT)(GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
static void (*p_glGetRenderbufferParameterivEXT)(GLenum target, GLenum pname, GLint * params);
GLboolean (*p_glIsFramebufferEXT)(GLuint framebuffer);
static void (*p_glBindFramebufferEXT)(GLenum target, GLuint framebuffer);
static void (*p_glDeleteFramebuffersEXT)(GLsizei n, const GLuint * framebuffers);
static void (*p_glGenFramebuffersEXT)(GLsizei n, GLuint * framebuffers);
static GLenum (*p_glCheckFramebufferStatusEXT)(GLenum target);
static void (*p_glFramebufferTexture1DEXT)(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
static void (*p_glFramebufferTexture2DEXT)(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
static void (*p_glFramebufferTexture3DEXT)(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
static void (*p_glFramebufferRenderbufferEXT)(GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
static void (*p_glGetFramebufferAttachmentParameterivEXT)(GLenum target, GLenum attachment, GLenum pname, GLint * params);
static void (*p_glGenerateMipmapEXT)(GLenum target);

static int p_inited = 0;

int
InitializeOpenGLExtensions(void)
{
    int n = 0;
    p_glIsRenderbufferEXT = (void *)wglGetProcAddress("glIsRenderbufferEXT");
    if (p_glIsRenderbufferEXT == NULL)
        n++;
    p_glBindRenderbufferEXT = (void *)wglGetProcAddress("glBindRenderbufferEXT");
    if (p_glBindRenderbufferEXT == NULL)
        n++;
    p_glDeleteRenderbuffersEXT = (void *)wglGetProcAddress("glDeleteRenderbuffersEXT");
    if (p_glDeleteRenderbuffersEXT == NULL)
        n++;
    p_glGenRenderbuffersEXT = (void *)wglGetProcAddress("glGenRenderbuffersEXT");
    if (p_glGenRenderbuffersEXT == NULL)
        n++;
    p_glRenderbufferStorageEXT = (void *)wglGetProcAddress("glRenderbufferStorageEXT");
    if (p_glRenderbufferStorageEXT == NULL)
        n++;
    p_glGetRenderbufferParameterivEXT = (void *)wglGetProcAddress("glGetRenderbufferParameterivEXT");
    if (p_glGetRenderbufferParameterivEXT == NULL)
        n++;
    p_glIsFramebufferEXT = (void *)wglGetProcAddress("glIsFramebufferEXT");
    if (p_glIsFramebufferEXT == NULL)
        n++;
    p_glBindFramebufferEXT = (void *)wglGetProcAddress("glBindFramebufferEXT");
    if (p_glBindFramebufferEXT == NULL)
        n++;
    p_glDeleteFramebuffersEXT = (void *)wglGetProcAddress("glDeleteFramebuffersEXT");
    if (p_glDeleteFramebuffersEXT == NULL)
        n++;
    p_glGenFramebuffersEXT = (void *)wglGetProcAddress("glGenFramebuffersEXT");
    if (p_glGenFramebuffersEXT == NULL)
        n++;
    p_glCheckFramebufferStatusEXT = (void *)wglGetProcAddress("glCheckFramebufferStatusEXT");
    if (p_glCheckFramebufferStatusEXT == NULL)
        n++;
    p_glFramebufferTexture1DEXT = (void *)wglGetProcAddress("glFramebufferTexture1DEXT");
    if (p_glFramebufferTexture1DEXT == NULL)
        n++;
    p_glFramebufferTexture2DEXT = (void *)wglGetProcAddress("glFramebufferTexture2DEXT");
    if (p_glFramebufferTexture2DEXT == NULL)
        n++;
    p_glFramebufferTexture3DEXT = (void *)wglGetProcAddress("glFramebufferTexture3DEXT");
    if (p_glFramebufferTexture3DEXT == NULL)
        n++;
    p_glFramebufferRenderbufferEXT = (void *)wglGetProcAddress("glFramebufferRenderbufferEXT");
    if (p_glFramebufferRenderbufferEXT == NULL)
        n++;
    p_glGetFramebufferAttachmentParameterivEXT = (void *)wglGetProcAddress("glGetFramebufferAttachmentParameterivEXT");
    if (p_glGetFramebufferAttachmentParameterivEXT == NULL)
        n++;
    p_glGenerateMipmapEXT = (void *)wglGetProcAddress("glGenerateMipmapEXT");
    if (p_glGenerateMipmapEXT == NULL)
        n++;
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
