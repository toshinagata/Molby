ifdef DESTDIR
OBJDIR = $(DESTDIR)/MolLib
else
DESTDIR = .
OBJDIR = .
endif

OBJECTS = ruby_types.o ruby_bind.o ruby_dialog.o ruby_md.o
ALL_OBJECTS = $(addprefix $(OBJDIR)/,$(OBJECTS))

$(OBJDIR)/%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(DESTDIR)/Ruby_bind.a : $(ALL_OBJECTS)
	$(AR) $(ARFLAGS) $(DESTDIR)/Ruby_bind.a $(ALL_OBJECTS)
	ranlib $(DESTDIR)/Ruby_bind.a

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(DESTDIR)/*.a
