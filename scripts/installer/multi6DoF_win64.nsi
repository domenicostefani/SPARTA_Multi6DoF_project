; Compile script using Nullsoft Scriptable Install System (NSIS) on windows

;--------------------------------
!include x64.nsh
!include "MUI2.nsh"

; load the version from file
; !define /file VERSION "../VERSION"
!define VERSION "1.0.0"

; The name of the installer
Name "multi6DoF_v${VERSION}_win64"

; The file to write
!system 'mkdir "../_WIN_RELEASE" 2> NUL'
!define OUTFILE "../_WIN_RELEASE/multi6DoF_v${VERSION}_win64.exe"
OutFile ${OUTFILE}

; Build Unicode installer
Unicode True

; The default installation directory
InstallDir $PROGRAMFILES64\Steinberg\VSTPlugins\multi6DoF_v${VERSION}_win64\"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------
; Pages
!insertmacro MUI_PAGE_WELCOME

!define MUI_TEXT_WELCOME_INFO_TITLE "Multi6DoF v${VERSION}"

!insertmacro MUI_PAGE_LICENSE "../../audio_plugins\_SPARTA_Multi6DoFconv_\README.md"
!insertmacro MUI_LANGUAGE "English"
Page components SelectComponents
Page directory
Page instfiles

Function SelectComponents 
  ; SetBrandingImage ..\..\sparta_screenshot.png
FunctionEnd

;--------------------------------

; The stuff to install
; Installer sections



; multi6DoF_VST2
Section "Multi6DoF VST2" multi6DoF_VST2
    SetOutPath "$INSTDIR"
    File /r "..\..\build\audio_plugins\_SPARTA_Multi6DoFconv_\sparta_multi6DoFconv_artefacts\Release\VST\sparta_multi6DoFconv.dll"
    
    ${DisableX64FSRedirection}
SectionEnd


; libfftw3f-3.dll
; Section "FFTW3f library" libfftw3f-3.dll
;     SetOutPath "$SYSDIR"
;     File "..\win-libs\x64\libfftw3f-3.dll"
; SectionEnd

; saf_ipp_custom.dll
Section "Spatial Audio Framework - IPP custom" saf_ipp_custom.dll
    SetOutPath "$SYSDIR"
;    File "C:\Windows\System32\saf_ipp_custom.dll"
    File "saf_ipp_custom.dll"
SectionEnd

; saf_mkl_custom_lp64.dll
Section "Spatial Audio Framework - MKL custom" saf_mkl_custom_lp64.dll
    SetOutPath "$SYSDIR"
;    File "C:\Windows\System32\saf_mkl_custom_lp64.dll"
    File "saf_mkl_custom_lp64.dll"
SectionEnd