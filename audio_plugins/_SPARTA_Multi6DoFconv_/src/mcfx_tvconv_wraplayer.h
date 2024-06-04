/*
 * Copyright 2020 Leo McCormack
 * Copyright 2024 Domenico Stefani
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Leo McCormack is responsible for developing ./SDKs/Spatial_Audio_Framework/examples/include/tvconv.h
 * Domenico Stefani is responsible for adapting the code for the integration of the MCFX Convolver (by Matthias Kronlachner)
 */

/**
 * @file mcfx_tvconv_wraplayer.h
 * @brief Wrapper/compatibility layer between the MCFX convolver and the tvconv API.
 * @author Domenico Stefani
 * @date 25 May 2024
 */

#ifndef __MCFX_TVCONV_WRAPLAYER_H_INCLUDED__
#define __MCFX_TVCONV_WRAPLAYER_H_INCLUDED__

#define PER_POS_CONVOLVER_MODE 1
#define SINGLE_CONVOLVER_MODE 2
#define CROSSFADED_CONVOLVERS_MODE 3

#define MCFX_CONVOLVER_MODE CROSSFADED_CONVOLVERS_MODE
#define PRIMING                    // Used to fill the convolver buffers with coherend data before switching to it
#define PRIMING_BUF_SIZE_BLOCKS 8  // This does something bad only if too small, but there is some other kind of discontinuity even when large
// TODO: determine priming buffer size with length of the longest IR PARETITION

// Next defines are here for debugging purposes
// #define DISABLE_FILTER_REPLACEMENT // [FOR DEBUG]
// #define DISABLE_ENTIRELY_CONV_CHANGE_BUT_RESET // [FOR DEBUG]
// #define DISABLE_CONV_RESET // [FOR DEBUG]
// #define ZEROS_AT_CHANGEBLOCK_END // [FOR DEBUG] Sets 4 samples at 0 at the end of priming blocks

#define CROSSFADE_RELEASE 1
#define CROSSFADE_DEBUG_MUTE 2     // [FOR DEBUG]
#define CROSSFADE_DEBUG_DISABLE 3  // [FOR DEBUG]
//---
#define CROSSFADE CROSSFADE_RELEASE
// #define CROSSFADE CROSSFADE_DEBUG_MUTE // [FOR DEBUG]
// #define CROSSFADE CROSSFADE_DEBUG_DISABLE // [FOR DEBUG]

// #define TEST_PRIMING_BUFFER_STATE  // [FOR DEBUG] This overwrites entire buffer with the priming buffer, to check if the priming buffer is correctly filled

#ifdef USE_ZITA_CONVOLVER
    #error "Not implemented yet (This was in MCFX but has not been correctly ported here)"
#endif

#include "_common.h"
#include "md_malloc.h"      //TODO: see whether to better integrate this with SAF
#include "saf_utilities.h"  //TODO: see whether to better integrate this with SAF
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    #define SAF_SOFA_READER_MODULE  // TODO: see whether to better integrate this with SAF// Or at least copy the licence comment from saf.h
    #include "saf_sofa_reader.h"    //TODO: see whether to better integrate this with SAF// Or at least copy the licence comment from saf.h
#endif                              /* SAF_ENABLE_SOFA_READER_MODULE */
#include <JuceHeader.h>
#include <stddef.h>  //TODO: see wether to keep or move (this is for the definition of NULL)

#include "MCFX_ConvolverData.h"
#include "MCFX_MtxConv.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** SOFA loader error codes */
typedef enum {
    /** Not initialized */
    SAF_TVCONV_NOT_INIT,
    /** Loading file */
    SAF_TVCONV_SOFA_LOADING,
    /** None of the error checks failed */
    SAF_TVCONV_SOFA_OK,
    /** Not a SOFA file, or no such file was found in the specified location */
    SAF_TVCONV_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH,
    /** Dimensions of the SOFA data were not as expected */
    SAF_TVCONV_SOFA_ERROR_DIMENSIONS_UNEXPECTED,
    /** The data-type of the SOFA data was not as expected */
    SAF_TVCONV_SOFA_ERROR_FORMAT_UNEXPECTED,
    /** NetCDF is not thread safe! */
    SAF_TVCONV_SOFA_ERROR_NETCDF_IN_USE

} SAF_TVCONV_ERROR_CODES;

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of McfxConv
 *
 * @param[in] phMcfxConv (&) address of McfxConv handle
 */
void tvconv_create(void** const phMcfxConv);

/**
 * Destroys an instance of McfxConv
 *
 * @param[in] phMcfxConv (&) address of McfxConv handle
 */
void tvconv_destroy(void** const phMcfxConv);

/**
 * Initialises an instance of McfxConv with default settings
 *
 * @param[in] hMcfxConv     McfxConv handle
 * @param[in] samplerate    Host samplerate.
 * @param[in] hostBlockSize Host frame/block size
 */
void tvconv_init(void* const hMcfxConv,
                 int samplerate,
                 int hostBlockSize);

void mcfxConv_process(void* const hMcfxConv,
                      juce::AudioBuffer<float>& buffer,
                      int nInputs,
                      int nOutputs,
                      bool isNonRealtime = true);

/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1. Re-initialising all settings/variables,
 * as McfxConv is currently configured, at next available opportunity.
 */
void tvconv_refreshParams(void* const hMcfxConv);

/**
 * Checks whether things have to be reinitialised, and does so if it is needed
 */
void tvconv_checkReInit(void* const hMcfxConv);

/** Reads IRs and positions from the current sofa file path. */
void tvconv_setFiltersAndPositions(void* const hMcfxConv);

/** Sets current sofa file path. */
void tvconv_setSofaFilePath(void* const hMcfxConv, const char* path);

/** Asks the MCFX engine to reinitialise the convolver, required after  */
void tvconv_reinitConvolver(void* const hMcfxConv);

/**
 *  Sets the target listener position.
 *
 *  @param[in] dim      dimension of the coordinate to be set (0 is x, 1 is y,
 *  *                   and 2 is z).
 *  @param[in] position new position to be set.
 */
void tvconv_setTargetPosition(void* const hMcfxConv, float position, int dim);

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int tvconv_getFrameSize(void);

/** Returns the number input channels */
int tvconv_getNumInputChannels(void* const hMcfxConv);

/** Returns the number of output channels (the same as the number of channels in the loaded sofa file) */
int tvconv_getNumOutputChannels(void* const hMcfxConv);

/** Returns the currect host block size */
int tvconv_getHostBlockSize(void* const hMcfxConv);

/** Returns the number of IR channels in the loaded sofa file */
int tvconv_getNumIRs(void* const hMcfxConv);

/** Returns the number of listener positions in the loaded sofa file */
int tvconv_getNumListenerPositions(void* const hMcfxConv);

/** Returns the current coordinate of dimension dim  (0 ... NUM_DIMENSIONS-1) */
float tvconv_getListenerPosition(void* const hMcfxConv, int index, int dim);

/** Returns the index of the current IR position */
int tvconv_getListenerPositionIdx(void* const hMcfxConv);

/** Returns the current coordinate of dimension dim  (0 ... NUM_DIMENSIONS-1) */
float tvconv_getTargetPosition(void* const hMcfxConv, int dim);

/** Returns the source coordinate of dimension dim  (0 ... NUM_DIMENSIONS-1) */
float tvconv_getSourcePosition(void* const hMcfxConv, int dim);

/** Returns minimum cooridinate of dimension dim (0 ... NUM_DIMENSIONS-1) */
float tvconv_getMinDimension(void* const hMcfxConv, int dim);

/** Returns minimum cooridinate of dimension dim  (0 ... NUM_DIMENSIONS-1) */
float tvconv_getMaxDimension(void* const hMcfxConv, int dim);

/** Returns the current filter length, in samples */
int tvconv_getIRLength(void* const hMcfxConv);

/** Returns the samplerate of the loaded filters  */
int tvconv_getIRFs(void* const hMcfxConv);

/* Returns true if the MCFX engine had to resample IRs to match the host samplerate */
int tvconv_getResampledIR(void* const hMcfxConv);

/** Returns the samperate of the host */
int tvconv_getHostFs(void* const hMcfxConv);

/** Returns the processing delay in samples (may be used for delay compensation features) */
int tvconv_getProcessingDelay(void* const hMcfxConv);

/** Returns the current Sofa file path */
char* tvconv_getSofaFilePath(void* const hMcfxConv);

/** Returns the current Sofa file error state */
SAF_TVCONV_ERROR_CODES tvconv_getSofaErrorState(void* const hMcfxConv);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS tvconv_getCodecStatus(void* const hMcfxConv);

/** Returns the number of skipped processing cycles */
int mcfxConv_getSkippedCyclesCount(void* const hMcfxConv);

/** Return the minimum number of input channels required by the MCFX convolver engine for the current SOFA configuration */
int mcfxConv_getMinInCh(void* const hMcfxConv);

/** Return the minimum number of output channels required by the MCFX convolver engine for the current SOFA configuration */
int mcfxConv_getMinOutCh(void* const hMcfxConv);

/** Return the number of convolutions computed by the MCFX convolver engine for the current SOFA configuration  */
int mcfxConv_getNumConv(void* const hMcfxConv);

/** Return the maximum partition size of the MCFX convolver engine for the current SOFA configuration */
unsigned int mcfxConv_getMaxPartitionSize(void* const hMcfxConv);

/** Return the buffer size of the MCFX convolver engine for the current SOFA configuration */
unsigned int mcfxConv_getBufferSize(void* const hMcfxConv);

/** Return the buffer size of the MCFX convolver engine for the current SOFA configuration */
unsigned int mcfxConv_getConvBufferSize(void* const hMcfxConv);

/** Set the buffer size of the MCFX convolver engine, which corresponds to the size of the first conv. partition in the non-uniform Conv. partitioning scheme*/
void mcfxConv_setConvBufferSize(void* const hMcfxConv, unsigned int convBufferSize);

/** Set the Maximum partition size (non uniform conv partitioning) of the MCFX convolver engine*/
void mcfxConv_setMaxPartitionSize(void* const hMcfxConv, unsigned int maxPartitionSize);

#if (MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE) && defined(MCFX_CONVOLVER_MODE) && defined(CROSSFADED_CONVOLVERS_MODE)  //[ds 2024]
float mcfxConv_getMaxCrossfadeTimeS(void* const hMcfxConv, bool* minReached = NULL, bool* maxReached = NULL);
void mcfxConv_DoubleCrossfadeTime(void* const hMcfxConv, bool* maxReached);
void mcfxConv_HalveCrossfadeTime(void* const hMcfxConv, bool* minReached);

void mcfxConv_setCrossfadeTime_samples(void* const hMcfxConv, unsigned int crossfadeTime_samples);
void mcfxConv_setCrossfadeTime_s(void* const hMcfxConv, float crossfadeTime_s);
void mcfxConv_setCrossfadeTime_ms(void* const hMcfxConv, float crossfadeTime_ms);
float mcfxConv_getCrossfadeTime_s(void* const hMcfxConv);
float mcfxConv_getCrossfadeTime_ms(void* const hMcfxConv);
unsigned int mcfxConv_getCrossfadeTime_samples(void* const hMcfxConv);
#endif

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MCFX_TVCONV_WRAPLAYER_H_INCLUDED__ */