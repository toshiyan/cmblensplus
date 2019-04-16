# -*- Makefile -*-

fac     = sharp
inc     = 
obj	= libsharp_all.o

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
###LIB:=$(OD)/libsharp_healpix_f$(LIBSUFFIX)
LIB:=$(LIBDIR)/libsharp_healpix_f$(LIBSUFFIX)

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

# create build sub-directory
%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi 

# object creation
$(OD)/%.o: $(SD)/%.c $(inc)
	$(CC) $(CFLAGS) -I c_utils -I libfftpack -I libsharp -I libsharp_healpix_f -c $< -o $@

