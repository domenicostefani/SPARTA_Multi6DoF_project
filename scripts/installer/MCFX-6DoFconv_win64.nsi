; this NSIS plugin provides a modern user interface
!include x64.nsh
!include "MUI2.nsh"

; this is the name of your software
Name "MCFX-6DoFconv"

; these variables will hold the vst installation paths the user picked
Var Vst3InstDir
Var Vst2InstDir


; Request application privileges for Windows Vista
RequestExecutionLevel admin

!define APP_NAME "MCFX-6DoFconv"
!define APP_VERSION "1.0.0"
!define APP_DIR "$PROGRAMFILES64\SPARTAMod\${APP_NAME}"
;!define APP_ICON "exampleplugin.ico"

; this is a customized finish page text
!define MUI_FINISHPAGE_TEXT "Thank you for installing ${APP_NAME}."

; icon for the installer and uninstaller (needs to be in the same folder where the makensis command is called)
;!define MUI_ICON "exampleplugin.ico"
;!define MUI_UNICON "exampleplugin.ico"

; your optional end user agreement page texts, not that later we will include a text file with the full license
!define MUI_LICENSE_PAGE_TITLE "End User License Agreement"
!define MUI_LICENSE_PAGE_TEXT "Please read the following license agreement carefully."
!define MUI_LICENSE_PAGE_CUSTOMTEXT "This license applies to both the binary and the VST3 plugin software."
!define MUI_TEXT_LICENSE_TITLE "End User License Agreement"

; defines the output name of the installer
Outfile "${APP_NAME}-${APP_VERSION}-Install.exe"
; defines the main installation folder (user can customize it)
InstallDir "${APP_DIR}"

; this ensures the installer has admin rights
RequestExecutionLevel admin

; this section is the main installing procedure
Section "${APP_NAME} Uninstaller" SEC01
    ${DisableX64FSRedirection}
    ; this part covers the app installation directory
    SetOutPath "$INSTDIR"
    ;File "..\..\build\audio_plugins\_SPARTA_MCFX-6DoFconv_\sparta_MCFX-6DoFconv_artefacts\Release\VST\sparta_MCFX-6DoFconv.dll"
    ;File "exampleplugin.ico"


    ; this will create the startmenu shortcut
    ;CreateShortCut "$SMSTARTUP\ExampleApp.lnk" "$INSTDIR\ExampleApp.exe" "" "$INSTDIR\ExampleApp.exe,0"
    ;CreateShortCut "$SMPROGRAMS\MyCompany\ExampleApp.lnk" "$INSTDIR\ExampleApp.exe" "" "$INSTDIR\exampleplugin.ico" 0

    ; this will write an uninstaller to the install location
    WriteUninstaller "$INSTDIR\${APP_NAME}-${APP_VERSION}-Uninstall.exe"
    CreateShortCut "$SMPROGRAMS\SPARTAMod\${APP_NAME}-${APP_VERSION}-Uninstall.lnk" "$INSTDIR\${APP_NAME}-${APP_VERSION}-Uninstall.exe" "" "" 0
SectionEnd

; MCFX-6DoFconv_VST3 section
Section "MCFX-6DoFconv VST3" MCFX-6DoFconv_VST3
    ; and this part will take care of copying the vst3 plugin to the location picked by user
    SetOutPath "$Vst3InstDir\sparta_MCFX-6DoFconv.vst3"
    ; the /r is for recursive (as it's a folder)
    File /nonfatal /a /r "..\..\build\audio_plugins\_SPARTA_MCFX-6DoFconv_\sparta_MCFX-6DoFconv_artefacts\Release\VST3\sparta_MCFX-6DoFconv.vst3\" #note back slash at the end
SectionEnd

; MCFX-6DoFconv_VST2 section
Section "MCFX-6DoFconv VST2" MCFX-6DoFconv_VST2
    ; and this part will take care of copying the vst2 plugin to the location picked by user
    SetOutPath "$Vst2InstDir\${APP_NAME}"
    File /nonfatal /a "..\..\build\audio_plugins\_SPARTA_MCFX-6DoFconv_\sparta_MCFX-6DoFconv_artefacts\Release\VST\sparta_MCFX-6DoFconv.dll"
SectionEnd


; libfftw3f-3.dll
Section "FFTW3f library" libfftw3f-3.dll
    SetOutPath "$SYSDIR"
    File libfftw3f-3.dll
SectionEnd

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


; here we describe all the files to be removed during uninstallation
Section "Uninstall"
   ;Delete "$INSTDIR\ExampleApp.exe"
   ;Delete "$INSTDIR\exampleplugin.ico"
   ;Delete "$SMPROGRAMS\MyCompany\ExampleApp.lnk"
   ;Delete "$SMSTARTUP\ExampleApp.lnk"
   Delete "$INSTDIR\${APP_NAME}-${APP_VERSION}-Uninstall.exe"
   RMDir /r "$Vst3InstDir\sparta_MCFX-6DoFconv.vst3"
   RMDir /r "$Vst2InstDir\${APP_NAME}"
   RMDir "$INSTDIR"
SectionEnd

; this instantiates some variables and ensure folder exists, as well as setting 64bits paths
Function .onInit
   SetRegView 64
   CreateDirectory "$SMPROGRAMS\MyCompany"
   StrCpy $Vst3InstDir "$PROGRAMFILES64\Common Files\VST3"
   StrCpy $Vst2InstDir "$PROGRAMFILES64\Steinberg\VSTPlugins"
FunctionEnd

; this instantiates the same kind of things but for the uninstaller
Function un.onInit
   SetRegView 64
   StrCpy $Vst3InstDir "$PROGRAMFILES64\Common Files\VST3"
   StrCpy $Vst2InstDir "$PROGRAMFILES64\Steinberg\VSTPlugins"
FunctionEnd

; then we proceed to describe all the page with MUI plugin macros
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "../../audio_plugins\_SPARTA_MCFX-6DoFconv_\README.md"

Page components SelectComponents

Function SelectComponents 
  ; SetBrandingImage ..\..\sparta_screenshot.png
FunctionEnd

!insertmacro MUI_PAGE_DIRECTORY  # this open a first time the browser page to set INSTDIR (default variable)

; we then proceed to set variables that will apply to the next MUI_PAGE_DIRECTORY page.
; Note that we also change MUI_DIRECTORYPAGE_VARIABLE, which will cause the directory page
; to change this new variable instead of INSTDIR
!define MUI_PAGE_HEADER_TEXT "Steinberg VST3 plugin location"
!define MUI_PAGE_HEADER_SUBTEXT "Choose the folder in which to install the ${APP_NAME} VST3 plugin."
!define MUI_DIRECTORYPAGE_TEXT_TOP "The installer will install the VST3 plugin for ${APP_NAME} in the following folder. To install in a differenct folder, click Browse and select another folder. Click Next to continue."
!define MUI_DIRECTORYPAGE_VARIABLE $Vst3InstDir
!insertmacro MUI_PAGE_DIRECTORY # and here we open the directory page with our custom settings

; Now MUI_DIRECTORYPAGE_VARIABLE for the VST2 plugin
!define MUI_PAGE_HEADER_TEXT "VST2 plugin location"
!define MUI_PAGE_HEADER_SUBTEXT "Choose the folder in which to install the ${APP_NAME} VST2 plugin."
!define MUI_DIRECTORYPAGE_TEXT_TOP "The installer will install the VST2 plugin for ${APP_NAME} in the following folder. To install in a differenct folder, click Browse and select another folder. Click Next to continue."
!define MUI_DIRECTORYPAGE_VARIABLE $Vst2InstDir
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

; the same process will apply to the uninstaller.
; Note that we could very much have registered the VST3 install location in the Windows registry
; instead of asking for it again.
!insertmacro MUI_UNPAGE_WELCOME
!define MUI_PAGE_HEADER_TEXT "Steinberg VST3 plugin location"
!define MUI_PAGE_HEADER_SUBTEXT "Choose the folder the ${APP_NAME} VST3 plugin was installed to."
!define MUI_DIRECTORYPAGE_TEXT_TOP "The uninstaller will remove the VST3 plugin from the following folder. Click Next to continue."
!define MUI_DIRECTORYPAGE_VARIABLE $Vst3InstDir
!insertmacro MUI_UNPAGE_DIRECTORY

!define MUI_PAGE_HEADER_TEXT "VST2 plugin location"
!define MUI_PAGE_HEADER_SUBTEXT "Choose the folder the ${APP_NAME} VST2 plugin was installed to."
!define MUI_DIRECTORYPAGE_TEXT_TOP "The uninstaller will remove the VST2 plugin from the following folder. Click Next to continue."
!define MUI_DIRECTORYPAGE_VARIABLE $Vst2InstDir
!insertmacro MUI_UNPAGE_DIRECTORY


!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

; finally we tell the installer to default to english
!insertmacro MUI_LANGUAGE "English"
