
#!/bin/sh

# setup directory

cwd=$(pwd)

fftw=${cwd}/FFTW
cfitsio=${cwd}/cfitsio
healpix=${cwd}/Healpix
lapack=${cwd}/LAPACK95


# Helper function for downloading, extracting, and navigating back
download_and_extract() {

    local url="$1"
    local output_dir="$2"
    local strip_components="${3:-0}"

    mkdir -p "${output_dir}"
    cd "${output_dir}"
    wget -q "${url}" -O package.tar.gz
    tar -xzf package.tar.gz --strip-components="${strip_components}"
    rm -f package.tar.gz
    cd "${cwd}"

}

# Remove and create a clean directory
reset_dir() {
    local dir="$1"
    [ -d "${dir}" ] && rm -rf "${dir}"
    mkdir -p "${dir}"
}


for args in "$@"; do
    # FFTW
    case "${args}" in
    FFTW | all)
        echo '---- Install FFTW ----'
        reset_dir "${fftw}"
        download_and_extract "https://www.fftw.org/fftw-3.3.10.tar.gz" "${fftw}/build"
        cd "${fftw}/build/fftw-3.3.10"
        ./configure --prefix=${fftw}/ --enable-openmp --enable-shared --enable-static
        make && make install
        cd ${cwd}
        ;;
    # cfitsio
    cfitsio | all)
        echo '---- Installing CFITSIO ----'
        reset_dir "${cfitsio}"
        download_and_extract "https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.5.0.tar.gz" "${cfitsio}/build"
        cd "${cfitsio}/build/cfitsio-4.5.0"
        ./configure --prefix="${cfitsio}"
        make && make install
        cd "${cwd}"
        ;;
    # Healpix
    healpix | all)
        echo '---- Installing Healpix ----'
        reset_dir "${healpix}"
        download_and_extract "https://sourceforge.net/projects/healpix/files/Healpix_3.80/Healpix_3.80_2021Jun22.tar.gz" "${healpix}" 1
        export LD_LIBRARY_PATH="${healpix}/lib:${cfitsio}/lib:$LD_LIBRARY_PATH"
        export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | awk '!seen[$0]++' | tr '\n' ':' | sed 's/:$//')
        cd "${healpix}"
        FC=ifort CC=gcc FITSDIR="${cfitsio}/lib" ./configure --auto=f90
        make && make test
        cd "${cwd}"
        ;;
    lapack | all)
        echo '---- Installing LAPACK and LAPACK95 ----'
        reset_dir "${lapack}"
        mkdir -p "${lapack}/lib" "${lapack}/mod" "${lapack}/tmp0" "${lapack}/tmp1"
        wget -q "http://www.netlib.org/lapack/lapack.tgz" -O "${lapack}/tmp0/lapack.tgz"
        wget -q "http://www.netlib.org/lapack95/lapack95.tgz" -O "${lapack}/tmp1/lapack95.tgz"
        tar -xzf "${lapack}/tmp0/lapack.tgz" -C "${lapack}/tmp0" --strip-components=1
        tar -xzf "${lapack}/tmp1/lapack95.tgz" -C "${lapack}/tmp1" --strip-components=1
        cd "${lapack}/tmp0"
        cp INSTALL/make.inc.ifort make.inc
        sed -i "s/FFLAGS =/FFLAGS = -fPIC/g" make.inc
        sed -i "s/FFLAGS_NOOPT =/FFLAGS_NOOPT = -fPIC/g" make.inc
        make
        mv *.a "${lapack}/lib/"
        cd "${lapack}/tmp1"
        sed -i "s/ = f95 -free/ = ifort -fPIC/g" make.inc
        sed -i "s/ = f95 -fixed/ = ifort -fPIC/g" make.inc
        cd "${lapack}/tmp1/SRC"
        make single_double_complex_dcomplex
        mv ${lapack}/tmp1/lapack95.a ${lapack}/lib/liblapack95.a
        mv ${lapack}/tmp1/lapack95_modules/* ${lapack}/mod/
        rm -rf "${lapack}/tmp0" "${lapack}/tmp1"
        cd "${cwd}"
        ;;
    *)
        echo "Unknown argument: ${args}"
        ;;
    esac
done

