static void X(hp_map2alm) (int nside, int lmax, int mmax, FLT *map,
  FLT *alm, double *zbounds, double *wgt)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info *ainfo;

  CHECK_STACK_ALIGN(8);
  sharp_make_healpix_geom_info_2 (nside, wgt, zbounds[0], zbounds[1], &ginfo);
  sharp_make_rectangular_alm_info (lmax,mmax,1,&ainfo);

  sharp_execute(SHARP_MAP2ALM,0,&alm,&map,ginfo,ainfo,1,FLAG,NULL,NULL);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

static void X(hp_map2alm_pol) (int nside, int lmax, int mmax, FLT *map,
  FLT *alm, double *zbounds, double *wgt)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info *ainfo;
  ptrdiff_t npix = 12*(ptrdiff_t)nside*nside;
  void *mapptr[2], *almptr[2];

  CHECK_STACK_ALIGN(8);
  sharp_make_healpix_geom_info_2 (nside, wgt, zbounds[0], zbounds[1], &ginfo);
  sharp_make_rectangular_alm_info (lmax,mmax,3,&ainfo);

  sharp_execute(SHARP_MAP2ALM,0,&alm,&map,ginfo,ainfo,1,FLAG,0,0);
  mapptr[0]=map+npix; mapptr[1]=map+2*npix;
  almptr[0]=alm+2; almptr[1]=alm+4;
  sharp_execute(SHARP_MAP2ALM,2,&almptr[0],&mapptr[0],ginfo,ainfo,1,FLAG,0,0);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

static void X(hp_map2alm_spin) (int nside, int lmax, int mmax, int spin,
  FLT *map, FLT *alm, double *zbounds, double *wgt)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info *ainfo;
  ptrdiff_t npix = 12*(ptrdiff_t)nside*nside;
  void *mapptr[2], *almptr[2];

  CHECK_STACK_ALIGN(8);
  sharp_make_healpix_geom_info_2 (nside, wgt, zbounds[0], zbounds[1], &ginfo);
  sharp_make_rectangular_alm_info (lmax,mmax,2,&ainfo);

  mapptr[0]=map; mapptr[1]=map+npix;
  almptr[0]=alm; almptr[1]=alm+2;
  sharp_execute(SHARP_MAP2ALM,spin,&almptr[0],&mapptr[0],ginfo,ainfo,1,FLAG,0,0);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

static void X(hp_alm2map) (int nside, int lmax, int mmax, FLT *alm,
  FLT *map)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info *ainfo;

  CHECK_STACK_ALIGN(8);
  sharp_make_healpix_geom_info (nside, 1, &ginfo);
  sharp_make_rectangular_alm_info (lmax,mmax,1,&ainfo);

  sharp_execute(SHARP_ALM2MAP,0,&alm,&map,ginfo,ainfo,1,FLAG,0,0);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

static void X(hp_alm2map_pol) (int nside, int lmax, int mmax, FLT *alm,
  FLT *map)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info *ainfo;
  ptrdiff_t npix = 12*(ptrdiff_t)nside*nside;
  void *mapptr[2], *almptr[2];

  CHECK_STACK_ALIGN(8);
  sharp_make_healpix_geom_info (nside, 1, &ginfo);
  sharp_make_rectangular_alm_info (lmax,mmax,3,&ainfo);

  sharp_execute(SHARP_ALM2MAP,0,&alm,&map,ginfo,ainfo,1,FLAG,0,0);
  mapptr[0]=map+npix; mapptr[1]=map+2*npix;
  almptr[0]=alm+2; almptr[1]=alm+4;
  sharp_execute(SHARP_ALM2MAP,2,&almptr[0],&mapptr[0],ginfo,ainfo,1,FLAG,0,0);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

static void X(hp_alm2map_spin) (int nside, int lmax, int mmax, int spin,
  FLT *alm, FLT *map)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info *ainfo;
  ptrdiff_t npix = 12*(ptrdiff_t)nside*nside;
  void *mapptr[2], *almptr[2];

  CHECK_STACK_ALIGN(8);
  sharp_make_healpix_geom_info (nside, 1, &ginfo);
  sharp_make_rectangular_alm_info (lmax,mmax,2,&ainfo);

  mapptr[0]=map; mapptr[1]=map+npix;
  almptr[0]=alm; almptr[1]=alm+2;
  sharp_execute(SHARP_ALM2MAP,2,&almptr[0],&mapptr[0],ginfo,ainfo,1,FLAG,0,0);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

FCALLSCSUB7(X(hp_map2alm),Y(SHARP_HP_MAP2ALM_X),Z(sharp_hp_map2alm_x),INT,INT,INT,PFLT,PFLT,PDOUBLE,PDOUBLE)
FCALLSCSUB7(X(hp_map2alm_pol),Y(SHARP_HP_MAP2ALM_POL_X),Z(sharp_hp_map2alm_pol_x),INT,INT,INT,PFLT,PFLT,PDOUBLE,PDOUBLE)
FCALLSCSUB8(X(hp_map2alm_spin),Y(SHARP_HP_MAP2ALM_SPIN_X),Z(sharp_hp_map2alm_spin_x),INT,INT,INT,INT,PFLT,PFLT,PDOUBLE,PDOUBLE)
FCALLSCSUB5(X(hp_alm2map),Y(SHARP_HP_ALM2MAP_X),Z(sharp_hp_alm2map_x),INT,INT,INT,PFLT,PFLT)
FCALLSCSUB5(X(hp_alm2map_pol),Y(SHARP_HP_ALM2MAP_POL_X),Z(sharp_hp_alm2map_pol_x),INT,INT,INT,PFLT,PFLT)
FCALLSCSUB6(X(hp_alm2map_spin),Y(SHARP_HP_ALM2MAP_SPIN_X),Z(sharp_hp_alm2map_spin_x),INT,INT,INT,INT,PFLT,PFLT)
