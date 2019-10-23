import libbasic

def bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,pktype= 'T12',ltype= ''):
  """
  Usage:
    :bl0,bl1 = basic.bispec.bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,pktype,ltype):
  """
  return libbasic.bispec.bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,pktype,ltype)

def bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype= 'T12',ltype= ''):
  """
  Usage:
    :bc,bl0,bl1 = basic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype,ltype):
  """
  return libbasic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,pktype,ltype)

def bispeclens_gauss_bin(shap,bn,lmin,lmax,cl):
  """
  Usage:
    :bl = basic.bispec.bispeclens_gauss_bin(shap,bn,lmin,lmax,cl):
  """
  return libbasic.bispec.bispeclens_gauss_bin(shap,bn,lmin,lmax,cl)

def zpoints(zmin,zmax,zn,zspace= 1):
  """
  Usage:
    :z,dz = basic.bispec.zpoints(zmin,zmax,zn,zspace):
  """
  return libbasic.bispec.zpoints(zmin,zmax,zn,zspace)

def skewspeclens(cpmodel,model,z,dz,zn,zs,olmin,olmax,lmin,lmax,k,pk0,kn,sigma,W,pktype= 'T12'):
  """
  Usage:
    :skew = basic.bispec.skewspeclens(cpmodel,model,z,dz,zn,zs,olmin,olmax,lmin,lmax,k,pk0,kn,sigma,W,pktype):
  """
  return libbasic.bispec.skewspeclens(cpmodel,model,z,dz,zn,zs,olmin,olmax,lmin,lmax,k,pk0,kn,sigma,W,pktype)

