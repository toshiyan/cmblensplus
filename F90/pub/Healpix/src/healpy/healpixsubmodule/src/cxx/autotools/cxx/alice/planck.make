PKG:=alice

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=
CXXBIN:=alice3
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)

ALLOBJ:=alice3.o
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)


ODEP:=$(HDR_$(PKG)) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
BDEP:=$(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir

$(BINDIR)/alice3: $(OD)/alice3.o $(BDEP)

all_cxxbin+=$(CXXBIN)
