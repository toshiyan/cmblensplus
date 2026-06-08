rm -rf build
make -C fortran_internal/src_dft clean
make -C fortran_internal/src_utils clean
env -u PYTHONPATH python -m pip install --no-build-isolation --no-cache-dir -e ".[dev]"

