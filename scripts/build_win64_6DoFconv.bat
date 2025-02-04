@echo off
@REM NExt command is a shortcut to exit the script if any command fails. Put %e% at the end of each command to enable this feature.
set "e=|| exit /b"

setlocal enabledelayedexpansion
set ROOT_DIR=..\


pushd %ROOT_DIR%
    cmake -B build -G "Visual Studio 17 2022" -DSAF_ENABLE_SOFA_READER_MODULE=1 -DSAF_ENABLE_NETCDF=1  %e%
    cd build %e%
    "C:\Program Files\Microsoft Visual Studio\2022\Community\Msbuild\Current\Bin\MSBuild.exe" audio_plugins\_SPARTA_6DoFconv_\sparta_6DoFconv.sln /p:Configuration=Release /m /p:PostBuildEventUseInBuild=false %e%
    if not exist "VST" mkdir -p VST %e%

    copy audio_plugins\_SPARTA_6DoFconv_\sparta_6DoFconv_artefacts\Release\VST\sparta_6DoFconv.dll VST %e%
popd
