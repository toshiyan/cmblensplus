
env -u PYTHONPATH MPLBACKEND=Agg pytest -q test/smoke.py

env -u PYTHONPATH python - <<'PY'
import cmblensplus
print("package:", cmblensplus.__file__)

import cmblensplus.basic as basic
import cmblensplus.flatsky as flatsky

print("basic:", basic.__file__)
print("flatsky:", flatsky.__file__)
PY

env -u PYTHONPATH python - <<'PY' 
import cmblensplus.basic as basic
import cmblensplus.flatsky as flatsky

print("basic public names:")
for x in dir(basic):
    if not x.startswith("_"):
        print(" ", x)
                             
print("\nflatsky public names:")    
for x in dir(flatsky):
    if not x.startswith("_"):
        print(" ", x)
PY

env -u PYTHONPATH python - <<'PY'
import cmblensplus.basic as basic
import cmblensplus.flatsky as flatsky

assert hasattr(basic, "wigner_funcs")
assert hasattr(basic, "bispec")
assert hasattr(basic, "flat")
assert hasattr(flatsky, "utils")
assert hasattr(flatsky, "rec_lens")
assert hasattr(flatsky, "norm_lens")

print("smoke test passed")
PY
