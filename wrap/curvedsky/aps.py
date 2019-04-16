import devcurvedsky

def calccl_spc(lmax,alm1,alm2=0,norm=0):
  """
  Usage:
    - e.g., cl = devcurvedsky.aps.calccl_spc(lmax,alm1,alm2,norm)
  """
  if alm2==0: alm2= 0
  if norm==0: norm= 1
  return devcurvedsky.aps.calccl_spc(lmax,alm1,alm2,norm)

def calccl_dpc(lmax,alm1,alm2=0,norm=0):
  """
  Usage:
    - e.g., cl = devcurvedsky.aps.calccl_dpc(lmax,alm1,alm2,norm)
  """
  if alm2==0: alm2= 0
  if norm==0: norm= 1
  return devcurvedsky.aps.calccl_dpc(lmax,alm1,alm2,norm)

def alm2bcl(bn,lmax,alm1,alm2=0,oL=0,norm=0):
  """
  Usage:
    - e.g., cb = devcurvedsky.aps.alm2bcl(bn,lmax,alm1,alm2,oL,norm)
  """
  if alm2==0: alm2= 0
  if oL==0: oL= 0
  if norm==0: norm= 1
  return devcurvedsky.aps.alm2bcl(bn,lmax,alm1,alm2,oL,norm)

