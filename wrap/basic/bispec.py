import libbasic

def cl_flat(cpmodel,z,dz,zs,lmax,k,pk0,zn=0,kn=0,pktype='T12',cltype='kk',dNdz=None):
  """
  Compute lensing bispectrum analytically
 
  Args:
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :z[*zn*] (*double*): redshift points for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :k[*kn*] (*double*): k for the matter power spectrum in unit of [*h/Mpc*]
    :pk0[*kn*] (*double*): the linear matter power spectrum at z=0 in unit of [*Mpc^3/h^3*]

  Args(optional):
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :cltype (*str*): kk, gk, or gg
    :dNdz[*zn*] (*double*): redshift distribution of galaxy, only used when btype includes g

  Returns:
    :cl[*l*] (*double*): power spectrum from LSS contributions at [*lmin,lmax*]

  Usage:
    :cl = basic.bispec.cl_flat(cpmodel,z,dz,zs,lmax,k,pk0,zn,kn,pktype,cltype,dNdz):
  """
  if zn == 0: zn=len(z)
  if kn == 0: kn=len(k)
  if dNdz is None: dNdz = z*0.
  if len(dNdz) != zn: print('size of dNdz is strange')
  return libbasic.bispec.cl_flat(cpmodel,z,dz,zs,lmax,k,pk0,zn,kn,pktype,cltype,dNdz)

def bispeclens(shap,cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,zn=0,kn=0,lan=0.,kan=0.,pktype='T12',ltype='',btype='kkk',dNdz=None):
  """
  Compute lensing bispectrum analytically
 
  Args:
    :shap (*str*): shape of the bispectrum (equi, fold, sque, or isos)
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :k[*kn*] (*double*): k for the matter power spectrum in unit of [*h/Mpc*]
    :pk0[*kn*] (*double*): the linear matter power spectrum at z=0 in unit of [*Mpc^3/h^3*]

  Args(optional):
    :lan, kan (*double*): parameters for the modified gravity extension, default to lan=kan=1 (GR)
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :ltype (*str*): fullsky correction (full) or not
    :btype (*str*): bispectrum type, i.e., kkk (lens-lens-lens), gkk (density-lens-lens), ggk (density-density-lens), or ggg (density-density-density)
    :dNdz[*zn*] (*double*): redshift distribution of galaxy, only used when btype includes g

  Returns:
    :bl0[*l*] (*double*): lensing bispectrum from LSS contributions at [*lmin,lmax*]
    :bl1[*l*] (*double*): lensing bispectrum from post-Born contributions at [*lmin,lmax*]

  Usage:
    :bl0,bl1 = basic.bispec.bispeclens(shap,cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,lan,kan,zn,kn,pktype,ltype,btype,dNdz):
  """
  if zn == 0: zn=len(z)
  if kn == 0: kn=len(k)
  if dNdz is None: dNdz = z*0.
  if len(dNdz) != zn: print('size of dNdz is strange')
  return libbasic.bispec.bispeclens(shap,cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,lan,kan,zn,kn,pktype,ltype,btype,dNdz)

def bispeclens_bin(shap,cpmodel,model,z,dz,zs,lmin,lmax,bn,k,pk0,zn=0,kn=0,lan=0.,kan=0.,pktype='T12',btype='kkk',dNdz=None):
  """
  Compute binned lensing bispectrum analytically
 
  Args:
    :shap (*str*): shape of the bispectrum (equi, fold, sque, or isos)
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :bn (*int*): number of multipole bins
    :k[*kn*] (*double*): k for the matter power spectrum in unit of [*h/Mpc*]
    :pk0[*kn*] (*double*): the linear matter power spectrum at z=0 in unit of [*Mpc^3/h^3*]

  Args(optional):
    :lan, kan (*double*): parameters for the modified gravity extension, default to lan=kan=1 (GR)
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :btype (*str*): bispectrum type, i.e., kkk (lens-lens-lens), gkk (density-lens-lens), ggk (density-density-lens), or ggg (density-density-density)
    :dNdz[*zn*] (*double*): redshift distribution of galaxy, only used when btype includes g

  Returns:
    :bc[*bn*] (*double*): multipole bin centers
    :bl0[*bn*] (*double*): binned lensing bispectrum from LSS contributions
    :bl1[*bn*] (*double*): binned lensing bispectrum from post-Born contributions

  Usage:
    :bc,bl0,bl1 = basic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype,btype,dNdz):
  """
  if zn == 0: zn=len(z)
  if kn == 0: kn=len(k)
  if dNdz is None: dNdz = z*0.
  if len(dNdz) != zn: print('size of dNdz is strange')
  return libbasic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype,btype,dNdz)

def bispeclens_snr(cpmodel,model,z,dz,zs,lmin,lmax,cl,k,pk0,zn=0,kn=0,pktype='T12',btype='kkk',dNdz=None,cgg=None):
  """
  Compute SNR of lensing bispectrum analytically
 
  Args:
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :cl[*l*] (*int*): observed lensing spectrum at 0<=l<=lmax
    :k[*kn*] (*double*): k for the matter power spectrum in unit of [*h/Mpc*]
    :pk0[*kn*] (*double*): the linear matter power spectrum at z=0 in unit of [*Mpc^3/h^3*]

  Args(optional):
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :btype (*str*): bispectrum type, i.e., kkk (lens-lens-lens), gkk (density-lens-lens), ggk (density-density-lens), or ggg (density-density-density)
    :dNdz[*zn*] (*double*): redshift distribution of galaxy, only used when btype includes g
    :cgg[*l*] (*double*): observed galaxy spectrum

  Returns:
    :snr (*double*): total SNR

  Usage:
    :snr = basic.bispec.bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,pktype,btype,dNdz,cgg):
  """
  if zn == 0: zn=len(z)
  if kn == 0: kn=len(k)
  if dNdz is None: dNdz = z*0.
  if len(dNdz) != zn: print('size of dNdz is strange')
  if cgg is None: cgg = cl*0.
  if len(cgg) != lmax+1: print('size of cgg is strange')
  return libbasic.bispec.bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,pktype,btype,dNdz,cgg)

def bispeclens_gauss_bin(shap,bn,lmin,lmax,cl):
  """
  Compute binned bispectrum analytically for the quadratic gaussian model

  Args:
    :shap (*str*): shape of the bispectrum (equi, fold, sque, or isos)
    :bn (*int*): number of multipole bins
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :cl[*l*] (*double*): the power spectrum at [*0:lmax+1*]

  Returns:
    :bc[*bn*] (*double*): multipole bin centers
    :bl[*bn*] (*double*): binned bispectrum

  Usage:
    :bc,bl = basic.bispec.bispeclens_gauss_bin(shap,bn,lmin,lmax,cl):
  """
  return libbasic.bispec.bispeclens_gauss_bin(shap,bn,lmin,lmax,cl)

def zpoints(zmin,zmax,zn,zspace=1):
  """
  Precomputing interpolation points for z

  Args: 
    :zmin/zmax (*double*): minimum/maximum redshifts
    :zn (*int*): number of redshifts

  Args(optional):
    :zspace (*int*): type of spacing

  Returns:
    :z[*zn*] (*double*): redshifts
    :dz[*zn*] (*double*): redshift intervals

  Usage:
    :z,dz = basic.bispec.zpoints(zmin,zmax,zn,zspace):
  """
  return libbasic.bispec.zpoints(zmin,zmax,zn,zspace)

def skewspeclens(cpmodel,model,zmin,zmax,zn,zs,ols,lmin,lmax,k,pk0,pktype='T12',pb=True,theta=0.0,Om=0.3,H0=70.,w0=-1.,wa=0.,mnu=0.06,ns=0.965,bn=0,kn=0,verbose=True):
  """
 Compute skew spectrum using a matter bispectrum fitting formula

  Args:
    :cpmodel (*str*): cosmological parameter model (model0, modelw, modelp, or input)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :zmin/zmax (*double*): minimum/maximum z for z-integral
    :zn (*int*): number of redshifts for the z-integral
    :zs[*2*] (*double*): source redshifts where zs[*2*] is used for the squared map
    :lmin/lmax (*int*): minimum/maximum multipoles of alms included in the skew spectrum
    :ols[*bn*] (*int*): output multipoles to be computed for skew spectrum
    :k[*kn*] (*double*): k for the matter power spectrum [*h/Mpc*]
    :pk0 (*double*): the linear matter power spectrum at z=0 [*Mpc^3/h^3*]

  Args(optional):
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :theta (*double*): kappa map resolution in arcmin
    :pb (*bool*): with post-Born correction or not (default=True)
    :verbose (*bool*): output messages

  Returns:
    :skew (*double*): skew-spectrum

  Usage:
    :skew = basic.bispec.skewspeclens(cpmodel,model,zmin,zmax,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pktype,pb,Om,H0,w0,wa,mnu,ns,verbose):
  """
  if bn==0: bn=len(ols)
  if kn == 0: kn=len(k)
  return libbasic.bispec.skewspeclens(cpmodel,model,zmin,zmax,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pktype,pb,Om,H0,w0,wa,mnu,ns,verbose)

