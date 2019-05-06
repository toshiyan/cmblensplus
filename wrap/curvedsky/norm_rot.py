import libcurvedsky

def qeb(lmax,rlmin,rlmax,fCEE,OCEE,OCBB):
  """
  Usage:
    :Al = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,fCEE,OCEE,OCBB):
  """
  return libcurvedsky.norm_rot.qeb(lmax,rlmin,rlmax,fCEE,OCEE,OCBB)

