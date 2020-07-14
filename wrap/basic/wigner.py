import libbasic

def wigner_3j(l2,l3,m2,m3):
  """
  Compute wigner3j for all possible l1 where w3j = (j1,j2,j3/m1,m2,m3)
 
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

  Returns:
    :bl0[*l*] (*double*): lensing bispectrum from LSS contributions at [*lmin,lmax*]
    :bl1[*l*] (*double*): lensing bispectrum from post-Born contributions at [*lmin,lmax*]

  Usage:
    :w3j = basic.wigner.wigner_3j(l2,l3,m2,m3):
  """
  return libbasic.wigner.wigner_3j(l2,l3,m2,m3)

