# -*- Makefile -*-


incalm  = alm_map_ss_inc.F90 alm_map_dd_inc.F90

incpix = convert_inplace_1d_inc.f90 \
	convert_inplace_nd_inc.f90 \
	convert_nest2ring_nd_inc.f90 \
	convert_nest2ring_1d_inc.f90 \
	convert_ring2nest_nd_inc.f90 \
	convert_ring2nest_1d_inc.f90 \
	remove_dipole_inc.f90 \
	apply_mask_inc.f90 \
	pixel_routines.F90

incfits = fits_s_inc.F90 fits_d_inc.F90

incudgrade = udgrade_s_inc.f90 udgrade_d_inc.f90

incnumrec = source_gen_sort_inc.f90

obj0  = healpix_types.o indmed.o utilities.o
obj1  = bit_manipulation.o extension.o long_intrinsic.o obsolete.o ran_tools_dist.o rngmod.o
obj2  = healpix_fft.o misc_utils.o
obj3  = head_fits.o num_rec.o statistics.o paramfile_io.o
#obj3  = head_fits.o num_rec.o statistics_stack.o paramfile_io.o
obj4  = pix_tools.o
obj5  = coord_v_convert.o fitstools.o udgrade_nr.o
obj6  = mask_tools.o alm_tools.o
obj7  = healpix_modules.o

obj   = sharp_healpix_f.o healpix_sharp_f90.o $(obj0) $(obj1) $(obj2) $(obj3) $(obj4) $(obj5) $(obj6) $(obj7)

fac = mod
MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
LIB:=$(LIBDIR)/libhealpix$(LIBSUFFIX)

# add internal name of shared library
ifneq ($(strip $(FLAGNAMELIB)),)
NAMELIB = $(FLAGNAMELIB)$(notdir $(LIB))
else
NAMELIB =
endif

# library creation
.DEFAULT_GOAL := all

all: $(OD)_mkdir $(LIB)

$(LIB): $(OBJ) $(incpix) $(incalm) $(incfits) $(incudgrade)
	$(AR) $@ $(NAMELIB) $(OBJ)
	@for p in $(OD)/*.$(MOD); do \
	cp $$p $(INCDIR)/; \
	done

# create build sub-directory
%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi

# objects creation
$(obj1:%=$(OD)/%): $(obj0:%=$(OD)/%)
$(obj2:%=$(OD)/%): $(obj1:%=$(OD)/%)
$(obj3:%=$(OD)/%): $(obj2:%=$(OD)/%)
$(obj4:%=$(OD)/%): $(obj3:%=$(OD)/%)
$(obj5:%=$(OD)/%): $(obj4:%=$(OD)/%)
$(obj6:%=$(OD)/%): $(obj5:%=$(OD)/%)
$(obj7:%=$(OD)/%): $(obj6:%=$(OD)/%)

$(OD)/%.o: $(SD)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OD)/%.o: $(SD)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(MODDIR)$(OD)

$(OD)/%.o: $(SD)/%.F90
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

$(OD)/long_intrinsic.o: $(SD)/long_intrinsic.F90
	@echo compiling $(SD)/long_intrinsic.F90
	@case $(OS) in \
	Darwin ) \
	ln -s -f $< long_intrinsic.junk.F90; \
	$(FC) $(FFLAGS) $(PPFLAGS) $(FI8FLAG) -c long_intrinsic.junk.F90 -o $@ $(MODDIR)$(OD); \
	rm long_intrinsic.junk.F90;; \
	OSF ) \
	ln -s -f $< long_intrinsic.f90; \
	$(FC) $(FFLAGS) -cpp $(FI8FLAG) -c $(basename $<).f90 -o $@ $(MODDIR)$(OD) ;; \
	* ) \
	$(FC) $(FFLAGS) $(PPFLAGS) $(FI8FLAG) -c $< -o $@ $(MODDIR)$(OD) ;; \
	esac

$(OD)/pix_tools.o:pix_tools.F90  $(incpix)

$(OD)/alm_tools.o:alm_tools.F90  $(incalm)

$(OD)/fitstools.o: fitstools.F90 $(incfits)

$(OD)/udgrade_nr.o: udgrade_nr.F90 $(incudgrade)

$(OD)/num_rec.o: num_rec.F90 $(incnumrec)

$(OD)/mask_tools.o: mask_tools.F90


# source generation
alm_map_ss_inc.F90: alm_map_template.F90 gen_alm_code
	@./gen_alm_code

alm_map_dd_inc.F90: alm_map_template.F90 gen_alm_code alm_map_ss_inc.F90
	@./gen_alm_code

fits_s_inc.F90: fits_template.F90 gen_fits_code
	@./gen_fits_code

fits_d_inc.F90: fits_template.F90 gen_fits_code fits_s_inc.F90
	@./gen_fits_code

udgrade_s_inc.f90: udgrade_template.f90 gen_udgrade_code
	@./gen_udgrade_code

udgrade_d_inc.f90: udgrade_template.f90 gen_udgrade_code udgrade_s_inc.f90
	@./gen_udgrade_code



