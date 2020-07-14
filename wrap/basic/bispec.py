import libbasic

def bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan=0.,kan=0.,pktype='T12',ltype=''):
  """
  Compute lensing bispectrum analytically
 
  Args:
    :shap (*str*): shape of the bispectrum (equi, fold, sque, or isos)
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :zn (*int*): number of redshifts for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :k[*kn*] (*double*): k for the matter power spectrum
    :pk0 (*double*): the linear matter power spectrum at z=0
    :kn (*int*): size of k

  Args(optional):
    :lan, kan (*double*): parameters for the modified gravity extension, default to lan=kan=1 (GR)
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :ltype (*str*): fullsky correction (full) or not

  Returns:
    :bl0[*l*] (*double*): lensing bispectrum from LSS contributions at [*lmin,lmax*]
    :bl1[*l*] (*double*): lensing bispectrum from post-Born contributions at [*lmin,lmax*]

  Usage:
    :bl0,bl1 = basic.bispec.bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,pktype,ltype):
  """
  return libbasic.bispec.bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,pktype,ltype)

def bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan=0.,kan=0.,pktype='T12'):
  """
  Compute binned lensing bispectrum analytically
 
  Args:
    :shap (*str*): shape of the bispectrum (equi, fold, sque, or isos)
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :zn (*int*): number of redshifts for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :bn (*int*): number of multipole bins
    :k[*kn*] (*double*): k for the matter power spectrum
    :pk0 (*double*): the linear matter power spectrum at z=0
    :kn (*int*): size of k

  Args(optional):
    :lan, kan (*double*): parameters for the modified gravity extension, default to lan=kan=1 (GR)
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)

  Returns:
    :bc[*bn*] (*double*): multipole bin centers
    :bl0[*bn*] (*double*): binned lensing bispectrum from LSS contributions
    :bl1[*bn*] (*double*): binned lensing bispectrum from post-Born contributions

  Usage:
    :bc,bl0,bl1 = basic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype):
  """
  return libbasic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype)

def bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,pktype='T12'):
  """
  Compute SNR of lensing bispectrum analytically
 
  Args:
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :zn (*int*): number of redshifts for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zs[*3*] (*double*): source redshifts
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :cl[*l*] (*int*): observed lensing spectrum at 0<=l<=lmax
    :k[*kn*] (*double*): k for the matter power spectrum
    :pk0 (*double*): the linear matter power spectrum at z=0
    :kn (*int*): size of k

  Args(optional):
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)

  Returns:
    :snr (*double*): total SNR

  Usage:
    :snr = basic.bispec.bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,pktype):
  """
  return libbasic.bispec.bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,pktype)

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

def skewspeclens(cpmodel,model,z,dz,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta=0.0,pktype='T12',pb=True):
  """
 Compute skew spectrum using a matter bispectrum fitting formula

  Args:
    :cpmodel (*str*): cosmological parameter model (model0, modelw, or modelp)
    :model (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
    :z[*zn*] (*double*): redshift points for the z-integral
    :dz[*zn*] (*double*): interval of z
    :zn (*int*): number of redshifts for the z-integral
    :zs[*2*] (*double*): source redshifts where zs[*2*] is used for the squared map
    :lmin/lmax (*int*): minimum/maximum multipoles of the bispectrum
    :bn (*int*): number of multipoles for skew spectrum
    :ols[*bn*] (*int*): multipoles to be computed for skew spectrum
    :k[*kn*] (*double*): k for the matter power spectrum
    :pk0 (*double*): the linear matter power spectrum at z=0
    :kn (*int*): size of k

  Args(optional):
    :pktype (*str*): fitting formula for the matter power spectrum (Lin, S02 or T12)
    :theta (*double*): kappa map resolution in arcmin
    :pb (*bool*): with post-Born correction or not (default=True)

  Returns:
    :skew (*double*): skew-spectrum

  Usage:
    :skew = basic.bispec.skewspeclens(cpmodel,model,z,dz,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pktype,pb):
  """
  return libbasic.bispec.skewspeclens(cpmodel,model,z,dz,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pktype,pb)

