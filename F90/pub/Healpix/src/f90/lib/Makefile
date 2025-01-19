# -*- Makefile -*-

fac	= lib
inc     =
VPATH	= $(BINDIR)
#exec	= $(BINDIR)/$(fac)
#src	= $(fac).f90
modobj	= gifmod.o
gifobj	= gifout.o gifstr.o
libobj	= gd.o gdfontt.o gdfonts.o gdfontmb.o gdfontl.o gdfontg.o
obj 	= $(libobj) $(gifobj) $(modobj)

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
LIB:=$(LIBDIR)/libhpxgif$(LIBSUFFIX)

# add internal name of shared library
ifneq ($(strip $(FLAGNAMELIB)),)
NAMELIB = $(FLAGNAMELIB)$(notdir $(LIB))
else
NAMELIB =
endif

# library creation
.DEFAULT_GOAL := all

all: $(OD)_mkdir $(LIB)

$(LIB): $(OBJ)
	$(AR) $@ $(NAMELIB) $(OBJ)
	@cp -p $@ $(LIBDIR)/libgif$(LIBSUFFIX)
	@for p in $(OD)/*.$(MOD); do \
	cp $$p $(INCDIR)/; \
	done

# create build sub-directory
%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi


# objects creation
$(OD)/%.o: $(SD)/%.f90 $(inc)
	$(FC) $(FFLAGS) -c $< -o $@ $(MODDIR)$(OD)

$(OD)/%.o: $(SD)/%.c $(inc)
	$(CC) $(CFLAGS) -c $< -o $@



