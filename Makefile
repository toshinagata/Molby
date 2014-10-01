ifeq ($(TARGET_PLATFORM),MAC)
 WX_DIR = $(PWD)/../../wxWidgets-3.0.0
 WX_LIB_DIR = $(WX_DIR)/osx-build/lib
 WX_ARCH_DIR = $(WX_LIB_DIR)/wx/include/osx_cocoa-unicode-static-3.0
 WX_CPPFLAGS = -isystem $(WX_ARCH_DIR) -isystem $(WX_DIR)/include -D_FILE_OFFSET_BITS=64 -D__WXMAC__ -D__WXOSX__ -D__WXOSX_COCOA__
 WX_LDFLAGS = -L$(WX_LIB_DIR)  -framework IOKit -framework Carbon -framework Cocoa -framework AudioToolbox -framework System -framework OpenGL -framework QuickTime -lwx_osx_cocoau-3.0 -lwx_osx_cocoau_gl-3.0 -framework WebKit -lwxregexu-3.0 -lwxtiff-3.0 -lwxjpeg-3.0 -lwxpng-3.0 -lz -lpthread -liconv 
 CPP_EXTRA_FLAGS = -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5 -arch ppc -arch i386 -DUSE_RUBY=1 -g -isystem $(PWD)/../../fftw-3.3.2/osx-build/include
 LD_EXTRA_FLAGS = -framework Accelerate -framework GLUT -L$(PWD)/../../fftw-3.3.2/osx-build/lib -lfftw3
 RUBY_DIR = $(PWD)/../../ruby-1.8.7-p160
 RUBY_CFLAGS = -isystem $(RUBY_DIR)/osx-build/include
 RUBY_LDFLAGS = -L$(RUBY_DIR)/osx-build/lib -lruby-static -lenc -ltrans
 EXECUTABLE = Molby
 EXE_SUFFIX =
endif

ifeq ($(TARGET_PLATFORM),MSW)
 WX_DIR = $(PWD)/../../wxWidgets-3.0.0
 WX_LIB_DIR = $(WX_DIR)/msw-build/lib
 WX_ARCH_DIR = $(WX_LIB_DIR)/wx/include/msw-unicode-static-3.0
 WX_CPPFLAGS = -isystem $(WX_ARCH_DIR) -isystem $(WX_DIR)/include -D_LARGEFIILE_SOURCE=unknown -D__WXMSW__
 WX_LDFLAGS = -L$(WX_LIB_DIR) -Wl,--subsystem,windows -mwindows -lwx_mswu_gl-3.0 -lopengl32 -lglu32 -lwx_mswu-3.0 -lwxregexu-3.0 -lwxexpat-3.0 -lwxtiff-3.0 -lwxjpeg-3.0 -lwxpng-3.0 -lwxzlib-3.0 -lrpcrt4 -loleaut32 -lole32 -luuid -lwinspool -lwinmm -lshell32 -lcomctl32 -lcomdlg32 -ladvapi32 -lwsock32 -lgdi32
 CPP_EXTRA_FLAGS = -isystem $(PWD)/../../CLAPACK-3.1.1.1-mingw/INCLUDE -isystem $(PWD)/../../fftw-3.3.2/msw-build/include -I$(PWD)/../MolLib
 LD_EXTRA_FLAGS = -L$(PWD)/../../CLAPACK-3.1.1.1-mingw/lib -L$(PWD)/../../fftw-3.3.2/msw-build/lib -llapackMinGW -lblasMinGW -lf2c_nomain -lfftw3 -static-libgcc
 RUBY_DIR = $(PWD)/../../ruby-2.0.0-p353
 RUBY_CFLAGS = -isystem $(RUBY_DIR)/msw-build/include/ruby-2.0.0 -I$(RUBY_DIR) -I$(RUBY_DIR)/msw-build/include/ruby-2.0.0/i386-mingw32
 RUBY_LDFLAGS = -L$(RUBY_DIR)/msw-build/lib -lmsvcrt-ruby200-static -lmsvcrt-ruby200 -lws2_32 -lshlwapi -limagehlp -lenc -ltrans
 EXECUTABLE = _Molby.exe_
 FINAL_EXECUTABLE = Molby.exe
 EXE_SUFFIX = .exe
endif

WXLIB_LIST = core,base,gl,adv
OBJECTS = ConsoleFrame.o GlobalParameterFrame.o GlobalParameterFilesFrame.o MoleculeView.o MyApp.o MyCommand.o MyDocument.o MyGLCanvas.o MySlider.o MyClipboardData.o ProgressFrame.o MyListCtrl.o MyDocManager.o wxKillAddition.o RubyDialogFrame.o MyIPCSupport.o MyVersion.o MyThread.o MyProgressIndicator.o modalwindow.o MyTextCtrl.o
LIBS = MolLib.a Ruby_bind.a
RUBY_EXTLIB = scanf.rb

ifeq ($(TARGET_PLATFORM),MAC)
PRODUCT = Molby.app
else
PRODUCT_DIR = Molby
PRODUCT = $(PRODUCT_DIR)/$(EXECUTABLE)
endif

CPP = g++
CC = gcc

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
CFLAGS = $(CPPFLAGS) $(COPT) $(CPP_EXTRA_FLAGS) $(RUBY_CFLAGS) $(WX_CPPFLAGS)
LDFLAGS = $(WX_LDFLAGS) $(LD_EXTRA_FLAGS) $(RUBY_LDFLAGS)
export CFLAGS
export LDFLAGS
export DESTDIR
export CC
export TARGET_PLATFORM

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
EXTRA_OBJECTS = listctrl.o window_msw.o textctrl_msw.o OpenGL_extensions.o
RESOURCE = molby_rc.o
#  The following HOMETEMP kludges are to work around a bug where '#include "..."' 
#  does not work when the include path is on the C: drive whereas the source is 
#  on the Z: drive. 2009.7.24. Toshi Nagata
HOMETEMP = $(HOME)/__molby_temp_build__
$(DESTPREFIX)/$(RESOURCE) : molby.rc
	mkdir -p $(HOMETEMP)/msw_build $(HOMETEMP)/bitmaps
	cp molby.rc $(HOMETEMP)/msw_build
	cp ../bitmaps/*.ico $(HOMETEMP)/bitmaps
	(cd $(HOMETEMP)/msw_build; windres -i molby.rc -o molby_rc.o -I$(WX_DIR)/include)
	cp $(HOMETEMP)/msw_build/molby_rc.o $@
	rm -rf $(HOMETEMP)
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
$(DESTPREFIX)/$(EXECUTABLE) : $(DESTOBJECTS)
ifeq ($(TARGET_PLATFORM),MAC)
	sh ../record_build_date.sh --with-svn-status
endif
ifeq ($(TARGET_PLATFORM),MSW)
	sh ../record_build_date.sh
endif
	$(CC) -c buildInfo.c -o $(DESTPREFIX)/buildInfo.o $(CFLAGS)
	$(CPP) -o $@ $(DESTOBJECTS) $(DESTPREFIX)/buildInfo.o $(CFLAGS) $(LDFLAGS)

$(DESTPREFIX)/$(PRODUCT) : $(DESTPREFIX)/$(EXECUTABLE) ../Scripts/*.rb amber11 ortep3/ortep3$(EXE_SUFFIX)
ifeq ($(TARGET_PLATFORM),MAC)
	rm -rf $(DESTPREFIX)/$(PRODUCT)
	mkdir -p $(DESTPREFIX)/$(PRODUCT)/Contents/MacOS
	mkdir -p $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	cp -f Info.plist $(DESTPREFIX)/$(PRODUCT)/Contents
	echo -n "APPL????" > $(DESTPREFIX)/$(PRODUCT)/Contents/PkgInfo
	cp -r ../Scripts $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	cp -r amber11 $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	cp -r ortep3 $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	mkdir -p $(DESTPREFIX)/$(PRODUCT)/Contents/Resources/Scripts/lib
	for i in $(RUBY_EXTLIB); do cp $(RUBY_DIR)/lib/$$i $(DESTPREFIX)/$(PRODUCT)/Contents/Resources/Scripts/lib; done
	cp $(DESTPREFIX)/$(EXECUTABLE) $(DESTPREFIX)/$(PRODUCT)/Contents/MacOS
endif
ifeq ($(TARGET_PLATFORM),MSW)
	rm -rf $(DESTPREFIX)/$(PRODUCT_DIR)
	mkdir -p $(DESTPREFIX)/$(PRODUCT_DIR)
	cp $(DESTPREFIX)/$(EXECUTABLE) $(DESTPREFIX)/$(PRODUCT_DIR)/$(FINAL_EXECUTABLE)
	cp `which mingwm10.dll` $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r ../Scripts $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r amber11 $(DESTPREFIX)/$(PRODUCT_DIR)
	cp -r ortep3 $(DESTPREFIX)/$(PRODUCT_DIR)
	mkdir -p $(DESTPREFIX)/$(PRODUCT_DIR)/Scripts/lib
	for i in $(RUBY_EXTLIB); do cp $(RUBY_DIR)/lib/$$i $(DESTPREFIX)/$(PRODUCT_DIR)/Scripts/lib; done
endif

ifeq ($(TARGET_PLATFORM),MSW)
setup: build/release/$(PRODUCT_DIR)/$(FINAL_EXECUTABLE)
	mkdir -p ../Products
	(/c/Program\ Files\ \(x86\)/Inno\ Setup\ 5/iscc molby.iss && mv Output/SetupMolbyWin.exe ../Products)
	(cd build/release/$(PRODUCT_DIR) && rm -rf $(MAKEDIR)/../Products/MolbyWin.zip && zip -r $(MAKEDIR)/../Products/MolbyWin.zip * -x \*.DS_Store \*.svn*)
endif

clean:
	rm -f $(DESTPREFIX)/*.o $(DESTPREFIX)/*.a $(DESTPREFIX)/$(EXECUTABLE)
	rm -rf $(DESTPREFIX)/$(PRODUCT_DIR)
	rm -f $(DESTPREFIX)/MolLib/*.o
	rm -f $(DESTPREFIX)/MolLib/MD/*.o
	rm -f $(DESTPREFIX)/MolLib/Ruby_bind/*.o

#	rm -f $(EXECUTABLE) $(OBJECTS)
#	rm -rf $(PRODUCT)
#	cd ../MolLib; $(MAKE) clean
#	cd ../MolLib/Ruby_bind; $(MAKE) clean

