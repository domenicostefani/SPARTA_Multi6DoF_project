@echo off
setlocal enabledelayedexpansion
set ROOT_DIR=..\


pushd %ROOT_DIR%
    cmake -B build -G "Visual Studio 17 2022" -DSAF_ENABLE_SOFA_READER_MODULE=1 -DSAF_ENABLE_NETCDF=1
    cd build
    "C:\Program Files\Microsoft Visual Studio\2022\Community\Msbuild\Current\Bin\MSBuild.exe" audio_plugins\_SPARTA_Multi6DoFconv_\sparta_Multi6DoFconv.sln /p:Configuration=Release /m /p:PostBuildEventUseInBuild=false
    if not exist "VST" mkdir -p VST

    copy audio_plugins\_SPARTA_Multi6DoFconv_\sparta_Multi6DoFconv_artefacts\Release\VST\sparta_multi6DoFconv.dll VST
popd
