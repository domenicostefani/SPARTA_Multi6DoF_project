set ROOT_DIR=..\

pushd %ROOT_DIR%
    cmake -B build -G "Visual Studio 17 2022" -DSAF_ENABLE_SOFA_READER_MODULE=1 -DSAF_ENABLE_NETCDF=1  -DFFTW3F_LIBRARY=D:/develop-farina-proj/libs/win/fftw-3.3.5-dll64/libfftw3f-3.lib -DFFTW3_INCLUDE_DIR=D:/develop-farina-proj/libs/win/fftw-3.3.5-dll64
popd