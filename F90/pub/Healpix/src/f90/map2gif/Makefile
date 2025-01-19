# -*- Makefile -*-

fac	= map2gif
inc     = 
VPATH	= $(BINDIR)
exec	= $(BINDIR)/$(fac)
#src	= euler_matrix_new.f90 $(fac).f90
obj 	= euler_matrix_new.o   $(fac).o

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
LIB:=$(LIBDIR)/libhealpix$(LIBSUFFIX)


%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi 


$(OD)/%.o: $(SD)/%.f90 $(inc) $(OD)_mkdir
	$(FC) $(FFLAGS) -c $< -o $@ $(MODDIR)$(OD)

default:  $(exec) 

$(exec): $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $(exec) $(OBJ) $(LDFLAGS)

