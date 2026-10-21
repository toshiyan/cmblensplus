
import importlib as _importlib

_SUBMODULES = ('delens', 'bispec', 'flat', 'galaxy', 'cosmofuncs', 'wigner_funcs')
__all__ = list(_SUBMODULES)

def __getattr__(name):
    if name in _SUBMODULES:
        mod = _importlib.import_module(f"{__name__}.{name}")
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(set(globals()) | set(_SUBMODULES))
    
#from . import delens
#from . import bispec
#from . import flat
#from . import galaxy
#from . import cosmofuncs
#from . import wigner_funcs

#__all__ = ['delens','bispec','flat','galaxy','cosmofuncs','wigner_funcs']

