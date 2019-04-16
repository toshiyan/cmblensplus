#!/bin/sh

BINPATH=../$HEALPIX_TARGET/bin

time $BINPATH/syn_alm_cxx fwhm_arcmin=60 infile=cl.fits rand_seed=1234 nlmax=512 outfile=\!test.alm polarisation=true && \
time $BINPATH/alm2map_cxx nlmax=512 infile=test.alm outfile=\!test.fits nside=256 polarisation=true && \
time $BINPATH/map2tga test.fits test.tga -bar -title "Synthesized Map" && \
time $BINPATH/anafast_cxx nlmax=512 infile=test.fits outfile_alms=\!test2.alm outfile=\!test_cl.fits polarisation=true iter_order=3 && \
time $BINPATH/alm2map_cxx nlmax=512 infile=test2.alm outfile=\!test2.fits nside=256 polarisation=true && \
time $BINPATH/map2tga test2.fits test2.tga -bar -title "Reconstructed Map" && \
time $BINPATH/udgrade_cxx infile=test2.fits outfile=\!test3.fits polarisation=false nside=8 && \
time $BINPATH/map2tga test3.fits test3.tga -bar -title "Downgraded Map" && \
time $BINPATH/map2tga test3.fits test4.tga -bar -interpol -title "Downgraded, Interpolated Map" && \
time $BINPATH/alm2map_cxx nlmax=512 infile=test.alm outfile=!test4.fits nside=317 polarisation=true && \
time $BINPATH/map2tga test4.fits test5.tga -bar -title "Synthesized Map (Nside=317)" && \
time $BINPATH/median_filter_cxx test.fits '!test5.fits' 60 && \
time $BINPATH/map2tga test5.fits test6.tga -bar -title "Median-filtered map (1 degree)" && \
time $BINPATH/smoothing_cxx nlmax=512 infile=test.fits outfile=\!test7.fits polarisation=true fwhm_arcmin=300 && \
time $BINPATH/alice3 ell=200 kernel_steps=200 steps=400 in=test7.fits nside=256 out=\!test_alice3 && \
time $BINPATH/map2tga test_alice3_background.fits test_alice1.tga -pal 0 -bar -interpol -title "Alice background texture" && \
time $BINPATH/map2tga test_alice3_texture.fits test_alice2.tga -pal 0 -bar -interpol -title "Alice convolved texture" && \
time $BINPATH/map2tga test_alice3_mod_texture.fits test_alice3.tga -pal 0 -bar -interpol -title "Alice modulated convolved texture"

display test.tga test2.tga test3.tga test4.tga test5.tga test6.tga test_alice?.tga
