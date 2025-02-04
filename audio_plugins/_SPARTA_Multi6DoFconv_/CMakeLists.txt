
cmake_minimum_required(VERSION 3.15)

project(sparta_multi6DoFconv VERSION 1.0.0)
message(STATUS "  ${PROJECT_NAME}")



IF(UNIX AND NOT APPLE AND NOT ANDROID)
	SET(WITH_ZITA_CONVOLVER FALSE CACHE BOOL "build with zita-convolver (better performance under linux)")
ENDIF(UNIX AND NOT APPLE AND NOT ANDROID)

IF(NOT APPLE)
	SET (WITH_FFTW3 TRUE)
ENDIF(NOT APPLE)

################################################
# additional cmake find packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

IF (WITH_ZITA_CONVOLVER OR WIN32)
	find_package(FFTW3 COMPONENTS single)
ENDIF (WITH_ZITA_CONVOLVER OR WIN32)

IF(UNIX AND NOT APPLE AND NOT ANDROID)
	find_package(FFTW3 COMPONENTS single threads)
ENDIF(UNIX AND NOT APPLE AND NOT ANDROID)

include_directories(${FFTW3_INCLUDE_DIR})







# juce helper function to create a library, which is configured as a JUCE plugin
juce_add_plugin(${PROJECT_NAME}
    VERSION                     "${PROJECT_VERSION}"
    COMPANY_NAME                "AALTO"                         
    IS_SYNTH                    FALSE                               
    NEEDS_MIDI_INPUT            FALSE                      
    NEEDS_MIDI_OUTPUT           FALSE                  
    IS_MIDI_EFFECT              FALSE                       
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE           
    COPY_PLUGIN_AFTER_BUILD     TRUE     # Where does this go in Windows + Linux?
    PLUGIN_MANUFACTURER_CODE    RAPD                
    PLUGIN_CODE                 M6DF                     
    FORMATS                     VST      # ( VST AAX Unity VST3 AU AUv3 Standalone )
    PRODUCT_NAME                "${PROJECT_NAME}"
)   

# Generates "JuceHeader.h", and adds its path to the include paths of the target
juce_generate_juce_header(${PROJECT_NAME})  

# Add plugin source files
target_sources(${PROJECT_NAME} 
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PluginEditor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PluginProcessor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sceneView.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PluginEditor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PluginProcessor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sceneView.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MCFX_ConvolverData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MCFX_ConvolverData.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MCFX_MtxConv.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MCFX_MtxConv.h
)

# Add any extra JUCE-specific pre-processor definitions
target_compile_definitions(${PROJECT_NAME}
PUBLIC 
    JUCE_DISPLAY_SPLASH_SCREEN=0
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
)
 
# Link with saf example, and the required juce modules
target_link_libraries(${PROJECT_NAME} 
PRIVATE 
    saf_example_tvconv
    saf_example_rotator
    juce::juce_recommended_config_flags
    juce::juce_recommended_lto_flags
    juce::juce_audio_plugin_client
    juce::juce_audio_utils
    juce::juce_osc
)

IF(WITH_FFTW3)
	target_link_libraries (${PROJECT_NAME} PRIVATE
		${FFTW3F_LIBRARY}
		${FFTW3F_THREADS_LIBRARY}
	)
ENDIF(WITH_FFTW3)