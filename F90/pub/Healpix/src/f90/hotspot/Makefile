# -*- Makefile -*-

fac	= hotspot
inc     =
VPATH	= $(BINDIR)
exec	= $(BINDIR)/$(fac)
#src	= maxima_tools.f90 HotSpots.F90
obj 	= maxima_tools.o HotSpots.o

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
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

default:  $(exec)

$(exec): $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $(exec) $(OBJ) $(LDFLAGS)

# spell out dependencies to allow parallel compilation (make -j)
$(OD)/HotSpots.o: $(OD)/maxima_tools.o



