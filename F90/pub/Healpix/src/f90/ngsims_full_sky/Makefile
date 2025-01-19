# -*- Makefile -*-

fac1    = sky_ng_sim
fac2    = sky_ng_sim_bin
inc     =
VPATH	= $(BINDIR)
exec1	= $(BINDIR)/$(fac1)
exec2	= $(BINDIR)/$(fac2)
subs    = sky_sub.o  sub_ngpdf_powergauss.o  sub_ngpdf_sho.o
obj1 	= $(subs)  sky_ng_sim.o
obj2 	= $(subs)  sky_ng_sim_bin.o

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac1)
SD :=$(shell pwd)
OBJ1:=$(obj1:%=$(OD)/%)
OBJ2:=$(obj2:%=$(OD)/%)
SUBS:=$(subs:%=$(OD)/%)
LIB:=$(LIBDIR)/libhealpix$(LIBSUFFIX)


%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi


$(OD)/%.o: $(SD)/%.f90 $(inc) $(OD)_mkdir
	$(FC) $(FFLAGS) -c $< -o $@ $(MODDIR)$(OD)

$(OD)/%.o: $(SD)/%.F90 $(inc) $(OD)_mkdir
	@echo compiling $(SD)/$*.F90
	@case $(OS) in \
	Darwin ) \
	ln -s -f $< $*.junk.F90; \
	$(FC) $(FFLAGS) $(PPFLAGS) -c $*.junk.F90 -o $@ $(MODDIR)$(OD) ; \
	rm $*.junk.F90;; \
	OSF ) \
	ln -s -f $< $*.f90; \
	$(FC) $(FFLAGS) -cpp       -c $(basename $<).f90 -o $@ $(MODDIR)$(OD) ;; \
	* ) \
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@  $(MODDIR)$(OD);; \
	esac

default: $(exec1) $(exec2)

$(exec1): $(OBJ1) $(LIB)
	$(FC)  $(FFLAGS) -o $(exec1) $(OBJ1) $(PGLIBS) $(LDFLAGS)

$(exec2): $(OBJ2) $(LIB)
	$(FC)  $(FFLAGS) -o $(exec2) $(OBJ2) $(PGLIBS) $(LDFLAGS)

# spell out dependencies to allow parallel compilation (make -j)
$(OD)/sub_ngpdf_powergauss.o: $(OD)/sky_sub.o

$(OD)/sub_ngpdf_sho.o:        $(OD)/sky_sub.o

$(OD)/sky_ng_sim.o:           $(SUBS)

$(OD)/sky_ng_sim_bin.o:       $(SUBS)

