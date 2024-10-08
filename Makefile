ifeq ($(TARGET_PLATFORM),MSW)
 ifeq ($(TARGET_ARCH),x86_64)
  TOOL_PREFIX = x86_64-w64-mingw32-
  MSW_BUILD = build-win
  LIB_SUFFIX = -3.2-x86_64-w64-mingw32
  WINE_PATH=/Applications/Wine\ Stable.app/Contents/Resources/wine/bin
  PRODUCT_SUFFIX = 64
  TARGET_ARCH_DEFINE = -DTARGET_ARCH=64
  SETUP_NAME = SetupMolbyWin
 else
  TOOL_PREFIX = i686-w64-mingw32-
  MSW_BUILD = build-win32
  LIB_SUFFIX = -3.2-i686-w64-mingw32
  WINE_PATH=/Applications/Wine\ Stable.app/Contents/Resources/wine/bin
  PRODUCT_SUFFIX = 32
  TARGET_ARCH_DEFINE = -DTARGET_ARCH=32
  SETUP_NAME = SetupMolbyWin32
 endif
 WX_DIR = $(PWD)/../../wxWidgets-3.2.5
 WX_LIB_DIR = $(WX_DIR)/$(MSW_BUILD)/lib
 WX_ARCH_DIR = $(WX_LIB_DIR)/wx/include/$(TOOL_PREFIX)msw-unicode-static-3.2
 WX_CPPFLAGS = -isystem $(WX_ARCH_DIR) -isystem $(WX_DIR)/include -D_LARGEFIILE_SOURCE=unknown -D__WXMSW__ $(TARGET_ARCH_DEFINE)
 WX_LDFLAGS = -L$(WX_LIB_DIR) -Wl,--subsystem,windows -mwindows $(WX_LIB_DIR)/libwx_mswu_gl$(LIB_SUFFIX).a -lopengl32 -lglu32 $(WX_LIB_DIR)/libwx_mswu$(LIB_SUFFIX).a -limm32 -lwxtiff$(LIB_SUFFIX) -lwxjpeg$(LIB_SUFFIX) -lwxpng$(LIB_SUFFIX) -lwxregexu$(LIB_SUFFIX) -lwxscintilla$(LIB_SUFFIX) -lwxexpat$(LIB_SUFFIX) -lwxzlib$(LIB_SUFFIX) -lrpcrt4 -loleaut32 -lole32 -luuid -luxtheme -lwinspool -lwinmm -lshell32 -lshlwapi -lcomctl32 -lcomdlg32 -ladvapi32 -lversion -lws2_32 -lgdi32 -loleacc -lwinhttp
 CPP_EXTRA_FLAGS = -isystem $(PWD)/../../CLAPACK-3.1.1.1-mingw/INCLUDE -isystem $(PWD)/../../fftw-3.3.2/$(MSW_BUILD)/include -I$(PWD)/../MolLib
 LD_EXTRA_FLAGS = -L$(PWD)/../../CLAPACK-3.1.1.1-mingw/$(MSW_BUILD)/lib -L$(PWD)/../../fftw-3.3.2/$(MSW_BUILD)/lib -llapackMinGW -lblasMinGW -lf2c_nomain -lfftw3 -static-libgcc -static-libstdc++ -Wl,-Bstatic,-lpthread
 RUBY_DIR = $(PWD)/../../ruby-2.0.0-p353
 ifeq ($(TARGET_ARCH),x86_64)
  RUBY_CFLAGS = -isystem $(RUBY_DIR)/$(MSW_BUILD)/include/ruby-2.0.0 -I$(RUBY_DIR) -I$(RUBY_DIR)/$(MSW_BUILD)/include/ruby-2.0.0/x64-mingw32
  RUBY_LDFLAGS = -L$(RUBY_DIR)/$(MSW_BUILD)/lib -lx64-msvcrt-ruby200-static -lws2_32 -lshlwapi -limagehlp -lenc -ltrans
 else
  RUBY_CFLAGS = -isystem $(RUBY_DIR)/$(MSW_BUILD)/include/ruby-2.0.0 -I$(RUBY_DIR) -I$(RUBY_DIR)/$(MSW_BUILD)/include/ruby-2.0.0/i386-mingw32
  RUBY_LDFLAGS = -L$(RUBY_DIR)/$(MSW_BUILD)/lib -lmsvcrt-ruby200-static -lws2_32 -lshlwapi -limagehlp -lenc -ltrans
 endif
 EXECUTABLE = _Molby.exe_
 FINAL_EXECUTABLE = Molby$(FINAL_EXECUTABLE_SUFFIX).exe
 EXE_SUFFIX = .exe
endif

WXLIB_LIST = core,base,gl,adv
OBJECTS = ConsoleFrame.o GlobalParameterFrame.o GlobalParameterFilesFrame.o MoleculeView.o MyApp.o MyCommand.o MyDocument.o MyGLCanvas.o MySlider.o MyClipboardData.o ProgressFrame.o MyListCtrl.o MyDocManager.o wxKillAddition.o RubyDialogFrame.o MyIPCSupport.o MyVersion.o MyThread.o MyProgressIndicator.o MyToggleButton.o modalwindow.o MyTextCtrl.o
LIBS = MolLib.a Ruby_bind.a
RUBY_EXTLIB = scanf.rb

ifeq ($(TARGET_PLATFORM),MAC)
PRODUCT = Molby.app
else
PRODUCT_DIR = Molby
PRODUCT = $(PRODUCT_DIR)/$(FINAL_EXECUTABLE)
endif

CPP = $(TOOL_PREFIX)g++
CC = $(TOOL_PREFIX)gcc
AR = $(TOOL_PREFIX)ar
RANLIB = $(TOOL_PREFIX)ranlib

ifeq ($(MAKECMDGOALS),debug)
 DEBUG = 1
endif

ifeq ($(DEBUG),1)
 DESTPREFIX = build/debug
 COPT = -O0 -g
else
 DESTPREFIX = build/release
 COPT = -O2 -g
endif
MAKEDIR = $(PWD)
DESTDIR = $(PWD)/$(DESTPREFIX)
CFLAGS = $(CPPFLAGS) $(COPT) $(RUBY_CFLAGS) $(WX_CPPFLAGS) $(CPP_EXTRA_FLAGS)
LDFLAGS = $(WX_LDFLAGS) $(RUBY_LDFLAGS) $(LD_EXTRA_FLAGS)
export CFLAGS
export LDFLAGS
export DESTDIR
export CC
export CPP
export AR
export TARGET_PLATFORM
export TARGET_ARCH
export RANLIB
export PWD

release: all

debug: all

all: $(DESTPREFIX) $(DESTPREFIX)/$(PRODUCT)

$(DESTPREFIX) :
	mkdir -p $(DESTPREFIX)

amber11 : ../amber11/src/antechamber/*.[ch] ../amber11/src/sqm/*.f ../amber11/src/config.h
	make -f ../Makefile_amber11

ortep3/ortep3$(EXE_SUFFIX) :
	make -f ../Makefile_ortep3

ifeq ($(TARGET_PLATFORM),MSW)
EXTRA_OBJECTS = OpenGL_extensions.o
RESOURCE = molby_rc.o
$(DESTPREFIX)/$(RESOURCE) : molby.rc
	$(TOOL_PREFIX)windres -i molby.rc -o $(DESTPREFIX)/$(RESOURCE) -I$(WX_DIR)/include
endif

depend: cleandep $(DESTPREFIX) $(OBJECTS:%.o=$(DESTPREFIX)/%.d) $(EXTRA_OBJECTS:%.o=$(DESTPREFIX)/%.d)
	cat $(DESTPREFIX)/*.d > build/Makefile.depend

cleandep:
	rm -f build/Makefile.depend

-include build/Makefile.depend

#  For unknown reasons, OpenGL_extensions.c needs to be compiled with -O0.
$(DESTPREFIX)/OpenGL_extensions.o : ../wxSources/OpenGL_extensions.c
	$(CC) -c $< -o $@ $(CPPFLAGS) -O0 -g $(CPP_EXTRA_FLAGS) $(RUBY_CFLAGS) $(WX_CPPFLAGS)

$(DESTPREFIX)/%.d : ../wxSources/%.cpp
	$(CPP) -MM $< >$@ $(subst -arch ppc,,$(CFLAGS))

$(DESTPREFIX)/%.d : ../wxSources/%.c
	$(CC) -MM $< >$@ $(subst -arch ppc,,$(CFLAGS))

$(DESTPREFIX)/%.o : ../wxSources/%.cpp
	$(CPP) -c $< -o $@ $(CFLAGS)

$(DESTPREFIX)/%.o : ../wxSources/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(DESTPREFIX)/MolLib.a : ../MolLib/*.[ch] ../MolLib/MD/*.[ch]
	mkdir -p $(DESTPREFIX)/MolLib/MD; cd ../MolLib; $(MAKE)

$(DESTPREFIX)/Ruby_bind.a : ../MolLib/Ruby_bind/*.[ch]
	mkdir -p $(DESTPREFIX)/MolLib/Ruby_bind; cd ../MolLib/Ruby_bind; $(MAKE)

ALL_OBJECTS = $(OBJECTS) $(EXTRA_OBJECTS) $(LIBS) $(RESOURCE)
DESTOBJECTS = $(addprefix $(DESTPREFIX)/,$(ALL_OBJECTS))
$(DESTPREFIX)/$(EXECUTABLE) : $(DESTOBJECTS) ../revisionInfo.txt
ifeq ($(TARGET_PLATFORM),MSW)
	sh ../record_build_date.sh --with-svn-status >buildInfo.c
endif
	$(CC) -c buildInfo.c -o $(DESTPREFIX)/buildInfo.o $(CFLAGS)
	$(CPP) -o $@ $(DESTOBJECTS) $(DESTPREFIX)/buildInfo.o $(CFLAGS) $(LDFLAGS)

$(DESTPREFIX)/$(PRODUCT) : $(DESTPREFIX)/$(EXECUTABLE) ../Scripts/*.rb ../bitmaps/bitmaps amber11 ortep3/ortep3$(EXE_SUFFIX)
ifeq ($(TARGET_PLATFORM),MSW)
	echo PWD = $(PWD)
	rm -rf $(DESTPREFIX)/$(PRODUCT_DIR)
	mkdir -p $(DESTPREFIX)/$(PRODUCT_DIR)
	cp $(DESTPREFIX)/$(EXECUTABLE) "$(DESTPREFIX)/$(PRODUCT_DIR)/$(FINAL_EXECUTABLE)"
	cp -r ../Scripts $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r ../bitmaps/bitmaps $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r amber11 $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r ortep3 $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r ../JANPA $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r ../docs/MolbyDoc $(DESTPREFIX)/$(PRODUCT_DIR)
	mkdir -p $(DESTPREFIX)/$(PRODUCT_DIR)/Scripts/lib
	for i in $(RUBY_EXTLIB); do cp $(RUBY_DIR)/lib/$$i $(DESTPREFIX)/$(PRODUCT_DIR)/Scripts/lib; done
endif

ifeq ($(TARGET_PLATFORM),MSW)
ifneq ($(DEBUG),1)
install: setup

setup: $(DESTPREFIX) $(DESTPREFIX)/$(PRODUCT_DIR)/$(FINAL_EXECUTABLE)
ifneq ($(DEBUG),1)
	for i in $(DESTPREFIX)/$(PRODUCT_DIR)/*.exe $(DESTPREFIX)/$(PRODUCT_DIR)/amber11/bin/*.exe $(DESTPREFIX)/$(PRODUCT_DIR)/ortep3/*.exe; do $(TOOL_PREFIX)strip "$$i"; done
endif
	mkdir -p ../latest_binaries
ifneq ($(WINE_PATH),)
	($(WINE_PATH)/wine ../../Inno\ Setup\ 5/ISCC.exe molby$(PRODUCT_SUFFIX).iss || exit 1)
else
	(/c/Program\ Files\ \(x86\)/Inno\ Setup\ 5/iscc molby$(PRODUCT_SUFFIX).iss || exit 1)
endif
	mv Output/$(SETUP_NAME).exe ../latest_binaries
	(cd build/release/$(PRODUCT_DIR) && rm -rf $(MAKEDIR)/../latest_binaries/MolbyWin$(PRODUCT_SUFFIX).zip && zip -r $(MAKEDIR)/../latest_binaries/MolbyWin$(PRODUCT_SUFFIX).zip * -x \*.DS_Store \*.svn*)
endif
endif

clean:
	rm -f $(DESTPREFIX)/*.o $(DESTPREFIX)/*.a $(DESTPREFIX)/$(EXECUTABLE)
	rm -rf $(DESTPREFIX)/$(PRODUCT_DIR)
	rm -f $(DESTPREFIX)/MolLib/*.o
	rm -f $(DESTPREFIX)/MolLib/MD/*.o
	rm -f $(DESTPREFIX)/MolLib/Ruby_bind/*.o

