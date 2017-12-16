ifdef DESTDIR
OBJDIR = $(DESTDIR)/MolLib
else
DESTDIR = .
OBJDIR = .
endif

OBJECTS = IntGroup.o MainView.o MainViewCommon.o MolAction.o Molecule.o \
  Object.o Dcd.o Parameter.o Trackball.o Types.o Missing.o
MD_OBJECTS = MDCore.o MDEwald.o MDForce.o MDGraphite.o MDPressure.o MDSurface.o

ALL_OBJECTS = $(addprefix $(OBJDIR)/,$(OBJECTS)) $(addprefix $(OBJDIR)/MD/,$(MD_OBJECTS))

$(OBJDIR)/%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(OBJDIR)/MD/%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(DESTDIR)/MolLib.a : $(ALL_OBJECTS)
	$(AR) $(ARFLAGS) $(DESTDIR)/MolLib.a $(ALL_OBJECTS)
	$(RANLIB) $(DESTDIR)/MolLib.a

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJDIR)/MD/*.o
	rm -f $(DESTDIR)/*.a
