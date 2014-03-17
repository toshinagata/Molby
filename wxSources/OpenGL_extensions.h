/*
 *  OpenGL_extensions.h
 *  Molby
 *
 *  Created by Toshi Nagata on 2014/03/17.
 *  Copyright 2014 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

/*  Implement OpenGL extensions  */
/*  (For Windows)  */

#ifndef __OpenGL_extensions_h__
#define __OpenGL_extensions_h__

#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#ifdef __cplusplus
extern "C" {
#endif

/*  Constants  */
#define GL_INVALID_FRAMEBUFFER_OPERATION_EXT               0x0506
#define GL_MAX_RENDERBUFFER_SIZE_EXT                       0x84E8
#define GL_FRAMEBUFFER_BINDING_EXT                         0x8CA6
#define GL_RENDERBUFFER_BINDING_EXT                        0x8CA7
#define GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT          0x8CD0
#define GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT          0x8CD1
#define GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL_EXT        0x8CD2
#define GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE_EXT 0x8CD3
#define GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_3D_ZOFFSET_EXT   0x8CD4
#define GL_FRAMEBUFFER_COMPLETE_EXT                        0x8CD5
#define GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT           0x8CD6
#define GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT   0x8CD7
#define GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT           0x8CD9
#define GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT              0x8CDA
#define GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT          0x8CDB
#define GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT          0x8CDC
#define GL_FRAMEBUFFER_UNSUPPORTED_EXT                     0x8CDD
#define GL_MAX_COLOR_ATTACHMENTS_EXT                       0x8CDF
#define GL_COLOR_ATTACHMENT0_EXT                           0x8CE0
#define GL_COLOR_ATTACHMENT1_EXT                           0x8CE1
#define GL_COLOR_ATTACHMENT2_EXT                           0x8CE2
#define GL_COLOR_ATTACHMENT3_EXT                           0x8CE3
#define GL_COLOR_ATTACHMENT4_EXT                           0x8CE4
#define GL_COLOR_ATTACHMENT5_EXT                           0x8CE5
#define GL_COLOR_ATTACHMENT6_EXT                           0x8CE6
#define GL_COLOR_ATTACHMENT7_EXT                           0x8CE7
#define GL_COLOR_ATTACHMENT8_EXT                           0x8CE8
#define GL_COLOR_ATTACHMENT9_EXT                           0x8CE9
#define GL_COLOR_ATTACHMENT10_EXT                          0x8CEA
#define GL_COLOR_ATTACHMENT11_EXT                          0x8CEB
#define GL_COLOR_ATTACHMENT12_EXT                          0x8CEC
#define GL_COLOR_ATTACHMENT13_EXT                          0x8CED
#define GL_COLOR_ATTACHMENT14_EXT                          0x8CEE
#define GL_COLOR_ATTACHMENT15_EXT                          0x8CEF
#define GL_DEPTH_ATTACHMENT_EXT                            0x8D00
#define GL_STENCIL_ATTACHMENT_EXT                          0x8D20
#define GL_FRAMEBUFFER_EXT                                 0x8D40
#define GL_RENDERBUFFER_EXT                                0x8D41
#define GL_RENDERBUFFER_WIDTH_EXT                          0x8D42
#define GL_RENDERBUFFER_HEIGHT_EXT                         0x8D43
#define GL_RENDERBUFFER_INTERNAL_FORMAT_EXT                0x8D44
#define GL_STENCIL_INDEX1_EXT                              0x8D46
#define GL_STENCIL_INDEX4_EXT                              0x8D47
#define GL_STENCIL_INDEX8_EXT                              0x8D48
#define GL_STENCIL_INDEX16_EXT                             0x8D49
#define GL_RENDERBUFFER_RED_SIZE_EXT                       0x8D50
#define GL_RENDERBUFFER_GREEN_SIZE_EXT                     0x8D51
#define GL_RENDERBUFFER_BLUE_SIZE_EXT                      0x8D52
#define GL_RENDERBUFFER_ALPHA_SIZE_EXT                     0x8D53
#define GL_RENDERBUFFER_DEPTH_SIZE_EXT                     0x8D54
#define GL_RENDERBUFFER_STENCIL_SIZE_EXT                   0x8D55

/*  Functions  */
extern GLboolean glIsRenderbufferEXT(GLuint renderbuffer);
extern void glBindRenderbufferEXT(GLenum target, GLuint renderbuffer);
extern void glDeleteRenderbuffersEXT(GLsizei n, const GLuint * renderbuffers);
extern void glGenRenderbuffersEXT(GLsizei n, GLuint * renderbuffers);
extern void glRenderbufferStorageEXT(GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
extern void glGetRenderbufferParameterivEXT(GLenum target, GLenum pname, GLint * params);
extern GLboolean glIsFramebufferEXT(GLuint framebuffer);
extern void glBindFramebufferEXT(GLenum target, GLuint framebuffer);
extern void glDeleteFramebuffersEXT(GLsizei n, const GLuint * framebuffers);
extern void glGenFramebuffersEXT(GLsizei n, GLuint * framebuffers);
extern GLenum glCheckFramebufferStatusEXT(GLenum target);
extern void glFramebufferTexture1DEXT(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
extern void glFramebufferTexture2DEXT(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
extern void glFramebufferTexture3DEXT(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
extern void glFramebufferRenderbufferEXT(GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
extern void glGetFramebufferAttachmentParameterivEXT(GLenum target, GLenum attachment, GLenum pname, GLint * params);
extern void glGenerateMipmapEXT(GLenum target);

extern int InitializeOpenGLExtensions(void);

#ifdef __cplusplus
}
#endif

#endif  /* __OpenGL_extensions_h__ */

