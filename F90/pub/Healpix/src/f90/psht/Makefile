# -*- Makefile -*-

libsrc	= libsharp_all.c
libobj	= libsharp_all.o

allobj	= $(libobj)

%.o: %.c
	$(CC) $(CFLAGS) -I c_utils -I libfftpack -I libsharp -I libsharp_healpix_f -c $< -o $@

default: libsharp_healpix_f.a

libsharp_healpix_f.a: $(libobj)
	$(AR) libsharp_healpix_f.a $(libobj)
	@cp -p libsharp_healpix_f.a $(LIBDIR)/
