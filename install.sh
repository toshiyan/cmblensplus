
#conda create -n py_v3_13 -c conda-forge python=3.13 numpy scipy matplotlib astropy healpy ducc0 meson ninja meson-python pkg-config compilers ipykernel pip tqdm

rm -rf build
make -C fortran_internal/src_dft clean
make -C fortran_internal/src_utils clean
env -u PYTHONPATH python -m pip install --no-build-isolation --no-cache-dir -e ".[dev]"

