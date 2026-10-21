# cmblensplus/flatsky/__init__.py
import importlib as _importlib

_SUBMODULES = ('utils','bispec','ffttools','rec_lens','norm_lens','rec_rot','norm_rot','rec_tau','norm_tau','rec_src','norm_src','norm_kxt','norm_kxs')
__all__ = list(_SUBMODULES)

def __getattr__(name):
    if name in _SUBMODULES:
        mod = _importlib.import_module(f"{__name__}.{name}")
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(set(globals()) | set(_SUBMODULES))
    
#from . import utils
#from . import bispec
#from . import ffttools
#from . import rec_lens
#from . import norm_lens
#from . import rec_rot
#from . import norm_rot
#from . import rec_tau
#from . import norm_tau
#from . import rec_src
#from . import norm_src
#from . import norm_kxt
#from . import norm_kxs

#__all__ = ["utils","bispec","ffttools","rec_lens","norm_lens","rec_rot","norm_rot","rec_tau","norm_tau","rec_src","norm_src","norm_kxt","norm_kxs",]
