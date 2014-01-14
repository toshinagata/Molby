ifeq ($(TARGET_PLATFORM),MAC)
 WXCONFIG_PREFIX = $(HOME)/Development/wxMac/osx-build/
 CPP_EXTRA_FLAGS = -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4 -arch ppc -arch i386 -DUSE_RUBY=1 -g -I$(PWD)/../../fftw-3.3.2/osx-build/include
 LD_EXTRA_FLAGS = -framework Accelerate -framework GLUT -L$(PWD)/../../fftw-3.3.2/osx-build/lib -lfftw3
 RUBY_DIR = $(HOME)/Development/ruby-1.8.7-static
 RUBY_CFLAGS = -I$(RUBY_DIR)
 RUBY_LDFLAGS = -L$(RUBY_DIR) -lruby-static
 EXECUTABLE = Molby
 EXE_SUFFIX =
endif

ifeq ($(TARGET_PLATFORM),MSW)
 WX_DIR = $(PWD)/../../wxWidgets-3.0.0
 WX_LIB_DIR = $(WX_DIR)/msw-build-3/lib
 WX_ARCH_DIR = $(WX_LIB_DIR)/wx/include/msw-unicode-static-3.0
 WX_CPPFLAGS = -I$(WX_ARCH_DIR) -I$(WX_DIR)/include -D_LARGEFIILE_SOURCE=unknown -D__WXMSW__
 WX_LDFLAGS = -L$(WX_LIB_DIR) -Wl,--subsystem,windows -mwindows -lwx_mswu_gl-3.0 -lopengl32 -lglu32 -lwx_mswu-3.0 -lwxregexu-3.0 -lwxexpat-3.0 -lwxtiff-3.0 -lwxjpeg-3.0 -lwxpng-3.0 -lwxzlib-3.0 -lrpcrt4 -loleaut32 -lole32 -luuid -lwinspool -lwinmm -lshell32 -lcomctl32 -lcomdlg32 -ladvapi32 -lwsock32 -lgdi32
 CPP_EXTRA_FLAGS = -O2 -I$(PWD)/../../CLAPACK-3.1.1.1-mingw/INCLUDE -I$(PWD)/../../fftw-3.3.2/msw-build-3/include
 LD_EXTRA_FLAGS = -L$(PWD)/../../CLAPACK-3.1.1.1-mingw/lib -L$(PWD)/../../fftw-3.3.2/msw-build-3/lib -llapackMinGW -lblasMinGW -lf2c_nomain -lfftw3
 RUBY_DIR = $(PWD)/../../ruby-1.8.7-p160
 RUBY_CFLAGS = -I$(RUBY_DIR)/msw-build-3/include
 RUBY_LDFLAGS = -L$(RUBY_DIR)/msw-build-3/lib -lmsvcrt-ruby18-static -lws2_32
 EXECUTABLE = _Molby.exe_
 FINAL_EXECUTABLE = Molby.exe
 EXE_SUFFIX = .exe
endif

WXLIB_LIST = core,base,gl,adv
OBJECTS = ConsoleFrame.o GlobalParameterFrame.o GlobalParameterFilesFrame.o MoleculeView.o MyApp.o MyCommand.o MyDocument.o MyGLCanvas.o MySlider.o MyClipboardData.o ProgressFrame.o MyListCtrl.o MyDocManager.o wxKillAddition.o RubyDialogFrame.o MyIPCSupport.o MyVersion.o MyThread.o MyProgressIndicator.o
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
CFLAGS = $(CPPFLAGS) $(CPP_EXTRA_FLAGS) $(RUBY_CFLAGS) $(WX_CPPFLAGS)
LDFLAGS = $(WX_LDFLAGS) $(LD_EXTRA_FLAGS) $(RUBY_LDFLAGS)
DESTPREFIX = build
DESTDIR = $(PWD)/$(DESTPREFIX)
export CFLAGS
export LDFLAGS
export DESTDIR
export CC

all: $(DESTPREFIX) $(DESTPREFIX)/$(PRODUCT)

$(DESTPREFIX) :
	mkdir -p $(DESTPREFIX)

amber11 : ../amber11/src/antechamber/*.[ch] ../amber11/src/sqm/*.f ../amber11/src/config.h
	make -f ../Makefile_amber11

$(DESTPREFIX)/mopac-build/mopac/mopac606$(EXE_SUFFIX) : 
	make -f ../Makefile_mopac606_nbo

ifeq ($(TARGET_PLATFORM),MSW)
EXTRA_OBJECTS = listctrl.o window_msw.o textctrl_msw.o
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
	cat $(DESTPREFIX)/*.d > $(DESTPREFIX)/Makefile.depend

cleandep:
	rm -f $(DESTPREFIX)/Makefile.depend

-include $(DESTPREFIX)/Makefile.depend

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
	$(CC) -c $(DESTPREFIX)/buildInfo.c -o $(DESTPREFIX)/buildInfo.o
	$(CPP) -o $@ $(DESTOBJECTS) $(DESTPREFIX)/buildInfo.o $(CFLAGS) $(LDFLAGS)

$(DESTPREFIX)/$(PRODUCT) : $(DESTPREFIX)/$(EXECUTABLE) ../Scripts/*.rb amber11 $(DESTPREFIX)/mopac-build/mopac/mopac606$(EXE_SUFFIX)
ifeq ($(TARGET_PLATFORM),MAC)
	rm -rf $(DESTPREFIX)/$(PRODUCT)
	mkdir -p $(DESTPREFIX)/$(PRODUCT)/Contents/MacOS
	mkdir -p $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	cp -f Info.plist $(DESTPREFIX)/$(PRODUCT)/Contents
	echo -n "APPL????" > $(DESTPREFIX)/$(PRODUCT)/Contents/PkgInfo
	cp -r ../Scripts $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	cp -r amber11 $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
	cp -r $(DESTPREFIX)/mopac-build/mopac $(DESTPREFIX)/$(PRODUCT)/Contents/Resources
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
	cp -r $(DESTPREFIX)/mopac-build/mopac $(DESTPREFIX)/$(PRODUCT_DIR)
	mkdir -p $(DESTPREFIX)/$(PRODUCT_DIR)/Scripts/lib
	for i in $(RUBY_EXTLIB); do cp $(RUBY_DIR)/lib/$$i $(DESTPREFIX)/$(PRODUCT_DIR)/Scripts/lib; done
endif

ifeq ($(TARGET_PLATFORM),MSW)
setup: $(DESTPREFIX)/$(PRODUCT_DIR)/$(FINAL_EXECUTABLE)
	/c/Program\ Files\ \(x86\)/Inno\ Setup\ 5/iscc molby.iss
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

