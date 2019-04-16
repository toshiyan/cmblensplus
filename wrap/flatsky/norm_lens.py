import libflatsky

def qtt(nx,ny,D,rL,OC,CT,eL):
  """
  Usage:
    - e.g., Ag,Ac = flatsky.norm_lens.qtt(nx,ny,D,rL,OC,CT,eL)
  """
  return flatsky.norm_lens.qtt(nx,ny,D,rL,OC,CT,eL)

def qte(nx,ny,D,rL,OT,OE,TE,eL):
  """
  Usage:
    - e.g., Ag,Ac = flatsky.norm_lens.qte(nx,ny,D,rL,OT,OE,TE,eL)
  """
  return flatsky.norm_lens.qte(nx,ny,D,rL,OT,OE,TE,eL)

def qtb(nx,ny,D,OT,OB,TE,rL,eL):
  """
  Usage:
    - e.g., Ag,Ac = flatsky.norm_lens.qtb(nx,ny,D,OT,OB,TE,rL,eL)
  """
  return flatsky.norm_lens.qtb(nx,ny,D,OT,OB,TE,rL,eL)

def qee(nx,ny,D,OE,EE,rL,eL):
  """
  Usage:
    - e.g., Ag,Ac = flatsky.norm_lens.qee(nx,ny,D,OE,EE,rL,eL)
  """
  return flatsky.norm_lens.qee(nx,ny,D,OE,EE,rL,eL)

def qeb(nx,ny,D,OE,OB,EE,rL,eL):
  """
  Usage:
    - e.g., Ag,Ac = flatsky.norm_lens.qeb(nx,ny,D,OE,OB,EE,rL,eL)
  """
  return flatsky.norm_lens.qeb(nx,ny,D,OE,OB,EE,rL,eL)

