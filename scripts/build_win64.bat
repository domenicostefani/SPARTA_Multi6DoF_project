echo off
setlocal enabledelayedexpansion
set ROOT_DIR=..\


pushd %ROOT_DIR%
    cmake -B build -G "Visual Studio 17 2022" -DSAF_ENABLE_SOFA_READER_MODULE=1 -DSAF_ENABLE_NETCDF=1
    cd build
    @REM "C:\Program Files\Microsoft Visual Studio\2022\Community\Msbuild\Current\Bin\MSBuild.exe" ALL_BUILD.vcxproj /p:Configuration=Release /m
    "C:\Program Files\Microsoft Visual Studio\2022\Community\Msbuild\Current\Bin\MSBuild.exe" sparta_head.sln /p:Configuration=Release /m /p:PostBuildEventUseInBuild=false
    if not exist "VST" mkdir -p VST
    @REM Cp all mathcing audio_plugins/_SPARTA*/*artefacts/Release/VST/*.dll to VST; do it only with CMD commands


    for /d %%i in (.\audio_plugins\_SPARTA*) do (
        for /d %%j in (%%i\*artefacts) do (
            @REM Put j/Release/VST/ inside a tml variable
            set "tml=%%j\Release\VST"
            @REM Echo
            @REM echo Copying from !tml!
            @REM Copy all .dll files from tml to VST
            for %%k in (!tml!\*.dll) do (
                @REM Create variable with filename only (basename)
                echo Copying %%k to .\VST
                copy "%%k" .\VST 1>NUL
            )
        )
    )

popd