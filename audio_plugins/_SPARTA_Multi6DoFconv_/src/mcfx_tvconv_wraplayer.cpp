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
 * Leo McCormack is responsible for developing SDKs\Spatial_Audio_Framework\examples\src\tvconv\tvconv.c , tvconv_internal.h,  tvconv_internal.c
 * Domenico Stefani is responsible for adapting the code for the integration of the MCFX Convolver (by Matthias Kronlachner)
 */

/**
 * @file mcfx_tvconv_wraplayer.cpp
 * @brief Wrapper/compatibility layer between the MCFX convolver and the tvconv API.
 * @author Domenico Stefani
 * @date 25 May 2024
 */

#include "mcfx_tvconv_wraplayer.h"


/* ZITA defines*/
#ifdef USE_ZITA_CONVOLVER //TODO: fix ZITA with new multi convolver paradigm
    #define CONVPROC_SCHEDULER_PRIORITY 0
    #define CONVPROC_SCHEDULER_CLASS SCHED_FIFO
    #define ZITA_THREAD_SYNC_MODE true
#endif

#define MAX_PART_SIZE 8192
#define MAX_NCONVOLVERS 21
#define DEF_NCONVOLVERS 5

// ------ Internal functions ------ //
void internal_setCodecStatus(void* const hTVCnv, CODEC_STATUS newStatus);
void internal_findNearestNeigbour(void* const hTVCnv);
void internal_setMinMaxDimensions(void* const hTVCnv);

struct CrossfadedConvInfo {
    int pos_idx; // The position index for which the convolver is still active
    bool isProcessing;
    bool isFadingIn, isFadingOut;
    int age; // how many samples the current convolver has been active for. Used for convolver "Stealing" when the oldest convolver is replaced
    float gain = 1.0f;
};

typedef struct _internal_MCFXConv_struct {
// int hopSize, fftSize, nBins;
// int length_h, nIRs, nCHout;
// int numFilterBlocks;
// void* hFFT;
// float* x_pad, *hx_n,
//         *z_n, *z_n_last, *z_n_last2,
//         *y_n_overlap, *y_n_overlap_last,
//         *out1, *out2,
//         *fadeIn, *fadeOut,
//         *outFadeIn, *outFadeOut;
// float_complex* X_n, *HX_n;
// float_complex*** Hpart_f;
// int posIdx_last, posIdx_last2;
#ifdef USE_ZITA_CONVOLVER //TODO: fix ZITA with new multi convolver paradigm
    Convproc* zita_conv; /* zita-convolver engine class instances */
#else

    #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
        juce::OwnedArray<MtxConvMaster>* mtxconv_s; // One convolver per position, no crossfade  (You can hear clicks when changing positions, plus very inefficient in terms of CPU & memory usage)
    #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
        MtxConvMaster* mtxconv; // Single convolver with no crossfade between positions (You can hear clicks when changing positions)
    #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
        juce::OwnedArray<MtxConvMaster>* mtxconv_xfd_s; // Multiple convolvers with crossfade between positions
        juce::OwnedArray<CrossfadedConvInfo>* mtxconv_xfd_info_s;
    #else
        #error "Invalid MCFX_CONVOLVER_MODE"
    #endif

#endif
    // Array with n Convolver data structures, each containing the filter data for a listener position
    juce::OwnedArray<ConvolverData>* convdata_s;
// ConvolverData conv_data;
    unsigned int _MaxPartSize;     // maximum size of the partition
    bool safemode_;  // this will add some latency for hosts that might send partial blocks, done automatically based on host type
    bool _configLoaded;  // is a configuration successfully loaded?
    Atomic<int> _skippedCycles;  // the number of skipped cycles do to unfinished partitions
    
    int _min_in_ch;
    int _min_out_ch;
    int _num_conv;

    int latencySamples = 0;


    int posIdx_last;

    #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
    #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
    #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
        int current_active_convolver_idx = -1;
    #else
        #error "Invalid MCFX_CONVOLVER_MODE"
    #endif
} Internal_Conv_struct;

// --- End of internal functions -- //


/* ====================================== */
/*  MCFX Create, Destroy, Apply functions */
/* ====================================== */

void mcfxConv_create(
    void** const phMCFXc,
    int _BufferSize,
    unsigned int& _ConvBufferSize_ref,
    float** H, /* nIRs x FLAT(nCHout x length_h) */
    int length_h,
    int nIRs,
    int nCHout,
    int initIdx,
    double sampleRate,
    int ir_fs,
    unsigned int& _MaxPartSize_ref
   #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
   #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
   #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
    ,int num_xfd_convolvers
   #else
    #error "Invalid MCFX_CONVOLVER_MODE"
   #endif
    ) {

    *phMCFXc = malloc1d(sizeof(Internal_Conv_struct));
    Internal_Conv_struct* h = (Internal_Conv_struct*)(*phMCFXc);

    /* INIT Mcfx params */
    h->_configLoaded = false;
    h->_skippedCycles = 0;   
    h->_min_in_ch = 0;
    h->_min_out_ch = 0;
    h->_num_conv = 0;
    /* End of Mcfx param init */

    h->posIdx_last = -1;
    
    // Ensure that the first partition size (or conv buffer size) is at least as large as the block size 
    if (_ConvBufferSize_ref < _BufferSize) _ConvBufferSize_ref = _BufferSize;
    _ConvBufferSize_ref = nextPowerOfTwo(_ConvBufferSize_ref);

    int nCHin = 1;  // TODO: pass this as argument as soon as multisource suppoert is enabled


    #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
            h->mtxconv_s = new juce::OwnedArray<MtxConvMaster>();
            for (int listenerPosition = 0; listenerPosition < nIRs; listenerPosition++) {
                h->mtxconv_s->add(new MtxConvMaster);
            }
    #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
        h->mtxconv = new MtxConvMaster;
    #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
        h->mtxconv_xfd_s = new juce::OwnedArray<MtxConvMaster>();
        h->mtxconv_xfd_info_s = new juce::OwnedArray<CrossfadedConvInfo>();

        for (int conv_idx = 0; conv_idx < num_xfd_convolvers; conv_idx++) {
            h->mtxconv_xfd_s->add(new MtxConvMaster);
            h->mtxconv_xfd_info_s->add(new CrossfadedConvInfo);
        }
    #else
        #error "Invalid MCFX_CONVOLVER_MODE"
    #endif

    h->convdata_s = new juce::OwnedArray<ConvolverData>();

    // Load the IR matrices in the way that the MtxConvMaster class expects
    // through the ConvolverData class
    for (int listenerPosition = 0; listenerPosition < nIRs; listenerPosition++) {
        h->convdata_s->add(new ConvolverData);
        ConvolverData* curConvData = h->convdata_s->getUnchecked(listenerPosition);
        curConvData->setSampleRate(sampleRate);

        float* curPositionH = H[listenerPosition];
        for (int out_ch = 0; out_ch < nCHout; out_ch++) {
            float* curOutchanH = curPositionH + out_ch * length_h;
            for (int in_ch = 0; in_ch < nCHin; in_ch++) {
                float* curSingleIR = curOutchanH;  // TODO: change this as soon as multipleSource support is added; use in_ch to address the correct IR
                std::cout << "Loading IR for -> posIdx: " << listenerPosition << " out_ch: " << out_ch << " in_ch: " << in_ch << std::endl;

                float* const* curSingleIRPtr = &curSingleIR;
                AudioBuffer<float> TempAudioBuffer(curSingleIRPtr, 1, length_h);  // This replaced the following line from MCFX: AudioBuffer<float> TempAudioBuffer(1,256); AND the loading of the ir from wav file

                float gain = 1.f;  // Feature of MCFX conf, not used in this implementation
                int delay = 0;     // Feature of MCFX conf, not used in this implementation

                {
                    double src_samplerate = (double)ir_fs;
                    // Differently from mcfx convolver, in_ch and out_ch are zero-based indexes so the first channel is 0, second 1 etc

                    curConvData->addIR(in_ch, out_ch, 0, delay, 0, TempAudioBuffer, 0, src_samplerate);

                    // String debug; //TODO: fix MCFX debug print
                    // debug << "conv # " << conv_data.getNumIRs() << " " << IrFilename.getFullPathName() << " (" << TempAudioBuffer.getNumSamples() << " smpls) loaded";
                    // if (src_samplerate != _SampleRate)
                    // debug << " (resampled to " << conv_data.getLength(conv_data.getNumIRs()-1) <<  " smpls)";
                    // DebugPrint(debug);
                    // this->_presetLoadState = PresetLoadState::Loaded;

                }


            } // end check channel assignment
            
        }
#ifdef USE_ZITA_CONVOLVER  // TODO: fix ZITA with new multi convolver paradigm
        int err = 0;

        unsigned int options = 0;

        options |= Convproc::OPT_FFTW_MEASURE;
        options |= Convproc::OPT_VECTOR_MODE;

        zita_conv.set_options(options);
        zita_conv.set_density(0.5);

        printf("max length: %lli \n", conv_data.getMaxLength());

        err = zita_conv.configure(conv_data.getNumInputChannels(), conv_data.getNumOutputChannels(), (unsigned int)conv_data.getMaxLength(), _BufferSize, _ConvBufferSize_ref, Convproc::MAXPART);

        for (int i = 0; i < conv_data.getNumIRs(); i++) {
            err = zita_conv.impdata_create(conv_data.getInCh(i), conv_data.getOutCh(i), 1, (float*)conv_data.getIR(i)->getReadPointer(0), 0, (unsigned int)conv_data.getLength(i));
        }

        zita_conv.print();
        zita_conv.start_process(CONVPROC_SCHEDULER_PRIORITY, CONVPROC_SCHEDULER_CLASS);

#else

        // _min_in_ch = conv_data.getNumInputChannels();
        // _min_out_ch = conv_data.getNumOutputChannels();
        // _num_conv = conv_data.getNumIRs();
        // Since we are here initializing multiple convolvers, we replace the previous three lines with setting the various values as a MAX of their current value and the new value
        // These should anyway be all equal in case of a properly created SOFA
        if (h->_min_in_ch!=0 && h->_min_in_ch != curConvData->getNumInputChannels())
            throw std::invalid_argument("Number of input channels in the SOFA file is not consistent across all listener positions");//TODO: Replace this with a non-fatal error
        h->_min_in_ch = jmax(h->_min_in_ch, curConvData->getNumInputChannels());
        if (h->_min_out_ch!=0 && h->_min_out_ch != curConvData->getNumOutputChannels())
            throw std::invalid_argument("Number of output channels in the SOFA file is not consistent across all listener positions");//TODO: Replace this with a non-fatal error
        h->_min_out_ch = jmax(h->_min_out_ch, curConvData->getNumOutputChannels());
        if (h->_num_conv!=0 && h->_num_conv != curConvData->getNumIRs())
            throw std::invalid_argument("Number of IRs in the SOFA file is not consistent across all listener positions");//TODO: Replace this with a non-fatal error
        h->_num_conv = jmax(h->_num_conv, curConvData->getNumIRs());
    }



    _MaxPartSize_ref = jmin(MAX_PART_SIZE, nextPowerOfTwo(_MaxPartSize_ref)); // Min, in case _MaxPartSize has been reduced by the user and reinit is performed as a consequence

    // try autodetecting host and deciding whether we need safemode (to avoid having to add another user parameter - let's see how this works for testers)
    PluginHostType me;
    h->safemode_ = me.isAdobeAudition() || me.isPremiere() || me.isSteinberg(); // probably an incomplete list



    #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
        int convolversToInit = nIRs;
    #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
        int convolversToInit = 1;
    #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
        int convolversToInit = num_xfd_convolvers;
    #else
     #error "Invalid MCFX_CONVOLVER_MODE"
    #endif



    // Here we initialize ALL the convolvers used.
    // if SINGLECONVOLVER is defined, only one convolver is initialized
    // if CROSSFADED_MULTIPLE_CONVOLVERS is defined, a set number of convolvers (< or > than the number of positions) is initialized (each conv with data from a different position)
    // if neither is defined, a convolver is initialized for each position, with the data from that position
    for (int conv_idx = 0; conv_idx < convolversToInit; conv_idx++) {
        int data_idx = conv_idx%h->convdata_s->size(); // in the case of CROSSFADED_MULTIPLE_CONVOLVERS and more convolvers than positions, we cycle through the positions
        ConvolverData* curConvData = h->convdata_s->getUnchecked(data_idx);


    #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
        MtxConvMaster* curMtxConv = h->mtxconv_s->getUnchecked(conv_idx);
    #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
        MtxConvMaster* curMtxConv = h->mtxconv;
    #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
        MtxConvMaster* curMtxConv = h->mtxconv_xfd_s->getUnchecked(conv_idx);

        CrossfadedConvInfo* curXfdInfo = h->mtxconv_xfd_info_s->getUnchecked(conv_idx);
        // CREATE CrossfadedConvInfo
        curXfdInfo->pos_idx = data_idx;
        curXfdInfo->isProcessing = false;
        curXfdInfo->isFadingIn = false;
        curXfdInfo->isFadingOut = false;
        curXfdInfo->age = 0;
        curXfdInfo->gain = 1.0f;
        // ---
    #else
        #error "Invalid MCFX_CONVOLVER_MODE"
    #endif

        curMtxConv->Configure(
            curConvData->getNumInputChannels(),
            curConvData->getNumOutputChannels(),
            _BufferSize, 
            curConvData->getMaxLength(),
            _ConvBufferSize_ref,
            _MaxPartSize_ref,
            h->safemode_);

        for (int i=0; i < curConvData->getNumIRs(); i++)
        {
            // if (threadShouldExit()) return; //TODO: remove thread stuff
            curMtxConv->AddFilter(curConvData->getInCh(i), curConvData->getOutCh(i), *curConvData->getIR(i));
            // no delay and length yet!
        }
        curMtxConv->StartProc(); //TODO: figure out when to enable
    }
#endif

    h->_configLoaded = true;
    h->latencySamples = h->safemode_ ? _ConvBufferSize_ref : _ConvBufferSize_ref-_BufferSize;
    h->_skippedCycles.set(0);
}

void mcfxConv_destroy(void** const phMCFXc) {
    Internal_Conv_struct* h = (Internal_Conv_struct*)(*phMCFXc);

    if (h != NULL) {
#ifdef USE_ZITA_CONVOLVER  // TODO: fix ZITA with new multi convolver paradigm
        h->zita_conv->stop_process();
        h->zita_conv->cleanup();
        delete h->zita_conv;
#else  // Use MtxConvolver

    #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
        // Stop processing with the convolver for each listener position (all convolvers)
        for (int listenerPosition = 0; listenerPosition < h->mtxconv_s->size(); listenerPosition++) {
            MtxConvMaster* curMtxConv = h->mtxconv_s->getUnchecked(listenerPosition);
            curMtxConv->StopProc();
            curMtxConv->Cleanup();
        }
        h->mtxconv_s->clear();  // Clear the array (which internally deletes the objects)
        delete h->mtxconv_s;    // Delete the array
    #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
        h->mtxconv->StopProc();
        h->mtxconv->Cleanup();
        delete h->mtxconv;
    #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
        // Stop processing for all convolvers
        for (int conv_idx = 0; conv_idx < h->mtxconv_xfd_s->size(); conv_idx++) {
            MtxConvMaster* curMtxConv = h->mtxconv_xfd_s->getUnchecked(conv_idx);
            curMtxConv->StopProc();
            curMtxConv->Cleanup();
        }
        h->mtxconv_xfd_s->clear();  // Clear convolver and info arrays (whoch internally deletes the objects)
        h->mtxconv_xfd_info_s->clear();
        delete h->mtxconv_xfd_s;  // Delete the arrays
        delete h->mtxconv_xfd_info_s;
    #else
        #error "Invalid MCFX_CONVOLVER_MODE"
    #endif

        h->convdata_s->clear();
        delete h->convdata_s;
#endif  // End ZITA/MtxConvolver ifdef
    }
    free(h);
    h = NULL;
    *phMCFXc = NULL;
}

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MIN_FRAME_SIZE ( 512 )
#define MAX_FRAME_SIZE ( 8192 )
#define NUM_DIMENSIONS ( 3 )

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/** Structure for a vector */
typedef float vectorND[NUM_DIMENSIONS];

/** Main structure for MCFX convolver  */
typedef struct _mcfxconv {
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];
    float outFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];

    /* Internal buffers */
    float** inputFrameTD;
    float** outputFrameTD;

    /* internal */
    void* hMCFXConv; /**< INTERNAL MCFX convolver handle */
    // [[deprecated("Replaced by _BufferSize from the MCFX convolver")]]
    // int hostBlockSize;          /**< current host block size */
    // [[deprecated("No longer used since MCXF convolver does not clamp the block size")]]
    // int hostBlockSize_clamped;  /**< Clamped between MIN and #MAX_FRAME_SIZE */
    // [[deprecated("Replaced by _SampleRate from the MCFX convolver")]]
    // int host_fs;                /**< current samplerate of the host */
    int reInitFilters;   /**< FLAG: 0: do not reinit, 1: reinit, 2: reinit in progress */
    int nOutputChannels; /**< number of output channels (same as the number of channels in the loaded wav) */

    int ir_fs;
    float** irs;     /**< npositionsx x (FLAT: nfilters x filter_length) */
    int nIrChannels; /**< number of filters per position */
    int ir_length;

    /* positions */
    vectorND* listenerPositions; /**< The listener positions; nListenerPositions x 3 */
    int nListenerPositions;
    vectorND minDimensions; /**< Minimum values across all dimensions */
    vectorND maxDimensions; /**< Maximum values across all dimensions */
    int position_idx;
    vectorND sourcePosition;

    /* flags/status */
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    PROC_STATUS procStatus;

    /* user parameters */
    int nInputChannels; /**< number of input channels */
    vectorND targetPosition;
    char* sofa_filepath;
    SAF_TVCONV_ERROR_CODES sofa_file_error;

    /*MCFX paramters*/
    double _SampleRate;
    unsigned int _BufferSize;      // size of the processing Block -> This replaces tvconv's hostBlockSize
    unsigned int _ConvBufferSize;  // size of the head convolution block (possibility to make it larger in order to reduce CPU load)

    bool _isProcessing;
    unsigned int _MaxPartSize;
    bool _paramReload;  // vst parameter to allow triggering reload of configuration

    Array<int> _conv_in;   // list with input routing
    Array<int> _conv_out;  // list with output routing

    /* Added when porting MCFX to 6DoF*/
    bool resampled_ir = false;

#if MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
    /* For Crossfade we use a number of convolvers that handle multiple convolutions while switching between matrices */
    int num_xfd_convolvers;
#endif
} McfxConvData;

void tvconv_create(void** const phMcfxConv) {
    McfxConvData* pData = (McfxConvData*)malloc1d(sizeof(McfxConvData));
    *phMcfxConv = (void*)pData;

    /* Default user parameters */
    pData->nInputChannels = 1;

    /* internal values */
    pData->_BufferSize = -1; /* force initialisation */
    pData->inputFrameTD = NULL;
    pData->outputFrameTD = NULL;
    pData->hMCFXConv = NULL;
    pData->irs = NULL;
    pData->reInitFilters = 1;
    pData->nIrChannels = 0;
    pData->ir_length = 0;
    pData->ir_fs = 0;
    pData->nOutputChannels = 0;
    pData->sofa_filepath = NULL;
    pData->sofa_file_error = SAF_TVCONV_NOT_INIT;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));

    /* positions */
    pData->listenerPositions = NULL;
    pData->nListenerPositions = 0;
    pData->position_idx = 0;
    for (int d = 0; d < NUM_DIMENSIONS; d++){
        pData->sourcePosition[d] = 0.0f;
        pData->targetPosition[d] = 0.0f;
        pData->minDimensions[d] = 0.0f;
        pData->maxDimensions[d] = 0.0f;
    }

    /* flags/status */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = (char *)malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;

    /* MCFX parameters */ 
    pData->_isProcessing = false;
    pData->_MaxPartSize = MAX_PART_SIZE;
    pData->_paramReload = false;
    
#if MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE //[ds 2024]
    pData->num_xfd_convolvers = DEF_NCONVOLVERS; // Default number of convolvers for crossfade
#endif
}

void tvconv_destroy(void** const phMcfxConv) {
    McfxConvData* pData = (McfxConvData*)(*phMcfxConv);

    if (pData != NULL){
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            if (pData->codecStatus == CODEC_STATUS_INITIALISING)
                std::cout << "pData->codecStatus == CODEC_STATUS_INITIALISING" << std::endl
                          << std::flush;
            if (pData->procStatus == PROC_STATUS_ONGOING)
                std::cout << "pData->procStatus == PROC_STATUS_ONGOING" << std::endl
                          << std::flush;
            SAF_SLEEP(10);
        }

        if (pData->inputFrameTD != NULL) free(pData->inputFrameTD);
        if (pData->outputFrameTD != NULL) free(pData->outputFrameTD);
        if (pData->irs != NULL) free(pData->irs);
        if (pData->listenerPositions != NULL) free(pData->listenerPositions);
        if (pData->hMCFXConv != NULL) mcfxConv_destroy(&(pData->hMCFXConv));
        free(pData);
        pData = NULL;
    }
}

void tvconv_init(
    void* const hConvData,
    int sampleRate,
    int hostBlockSize) {
    McfxConvData *pData = (McfxConvData*)(hConvData);

    pData->_SampleRate = sampleRate; // pData->_SampleRate replaced saf's original host_fs

    if(pData->_BufferSize != hostBlockSize){
        pData->_BufferSize = hostBlockSize;
        pData->_ConvBufferSize = hostBlockSize; //TODO: check how this works with SPARTAS way of handling blocksize in Init function calls 
        // pData->hostBlockSize_clamped = SAF_CLAMP(pData->_BufferSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);mcfxConv_replaceConvData
        pData->reInitFilters = 1;
        internal_setCodecStatus(hConvData, CODEC_STATUS_NOT_INITIALISED); 
    }
    tvconv_checkReInit(hConvData);
}

#if MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
/**
 * @brief Replace the current convolver data with the one at the new listener position, without allocating new memory
 * 
 * @param innerConvStruct 
 * @param newListenerPos 
 */
void mcfxConv_replaceConvData(void* const hMcfxConv, int newListenerPos) {
    if (hMcfxConv == NULL) return;
    McfxConvData* pData = (McfxConvData*)(hMcfxConv);
    if (pData->hMCFXConv == NULL) return;
    Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv);

    #if 0 //TODO: remove
    curMtxConv->Configure(
        curConvData->getNumInputChannels(),
        curConvData->getNumOutputChannels(),
        _BufferSize,
        curConvData->getMaxLength(),
        _ConvBufferSize_ref,
        _MaxPartSize_ref,
        h->safemode_);
    #endif
    if (newListenerPos >= h->convdata_s->size())
        return;

    ConvolverData* newPosConvData = h->convdata_s->getUnchecked(newListenerPos);
    MtxConvMaster* pMtxConv = h->mtxconv;

    for (int i = 0; i < newPosConvData->getNumIRs(); i++) {
        pMtxConv->ReplaceFilter(newPosConvData->getInCh(i), newPosConvData->getOutCh(i), *newPosConvData->getIR(i));
    }
}
#endif

void mcfxConv_process(void* const hMcfxConv,
                      juce::AudioBuffer<float>& buffer,
                      int nInputs,
                      int nOutputs) {


    McfxConvData* pData = (McfxConvData*)(hMcfxConv);

    int nSamples = buffer.getNumSamples();
    tvconv_checkReInit(hMcfxConv);

    // Process Block if filters are loaded and saf_matrixConv_apply is ready for it
    if (pData->reInitFilters == 0 && pData->codecStatus == CODEC_STATUS_INITIALISED) {
        pData->procStatus = PROC_STATUS_ONGOING; // TODO: check if this is the place for this line
        std::cout << "Currently Processing" << std::endl
                  << std::flush << std::fflush;
        Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv);

        //<-- From PluginProcessor.cpp of MCFXConvolver
        if (h != NULL && h->_configLoaded) { // TODO: make _configLoaded become true


        #if MCFX_CONVOLVER_MODE == PER_POS_CONVOLVER_MODE
            size_t convolverVectorSize = h->convdata_s->size();

            int current_pos_idx = pData->position_idx;
            if (current_pos_idx >= convolverVectorSize) {
                pData->procStatus = PROC_STATUS_NOT_ONGOING;
                return;  // TODO: check when this is the case
            }
            MtxConvMaster* curMtxConv = h->mtxconv_s->getUnchecked(current_pos_idx);
        #elif MCFX_CONVOLVER_MODE == SINGLE_CONVOLVER_MODE
            //// If a single convolver is in use, filters need to be replaced whenever the position changes // Cross-fading is not implemented yet

            MtxConvMaster* curMtxConv = h->mtxconv;
            int current_pos_idx = pData->position_idx;
            // If the current_pos_idx is different than the previous position, replace the filters and update previous position
            if (current_pos_idx != h->posIdx_last) {
                mcfxConv_replaceConvData((void*)pData, current_pos_idx);
                h->posIdx_last = current_pos_idx;
            }
        #elif MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE
                //TODO: Continue here with the crossfade implementation                
                MtxConvMaster* curMtxConv = h->mtxconv_xfd_s->getUnchecked(h->current_active_convolver_idx);
        #else
        #error "Invalid MCFX_CONVOLVER_MODE"
        #endif

            pData->_isProcessing = true;
            int curNumSamples = buffer.getNumSamples();

#ifdef USE_ZITA_CONVOLVER  // TODO: fix ZITA with new multi convolver paradigm
            for (int i = 0; i < jmin(h->conv_data.getNumInputChannels(), getTotalNumInputChannels()); i++) {
                float* indata = h->zita_conv->inpdata(i);
                memcpy(indata, buffer.getReadPointer(i), curNumSamples * sizeof(float));
            }

            h->zita_conv->process(ZITA_THREAD_SYNC_MODE);

            for (int i = 0; i < jmin(h->conv_data.getNumOutputChannels(), getNumOutputChannels()); i++) {
                float* outdata = h->zita_conv->outdata(i) + _ConvBufferPos;
                memcpy(buffer.getWritePointer(i), outdata, curNumSamples * sizeof(float));
            }
#else
            // curMtxConv->processBlock(buffer, buffer, isNonRealtime()); // if isNotRealtime always set to true!
            curMtxConv->processBlock(buffer, buffer, curNumSamples, true);  // try to always wait except - add a special flag to deactivate waiting...

            h->_skippedCycles.set(curMtxConv->getSkipCount());
#endif

            pData->_isProcessing = false;

        } else {  // config loaded

            // clear output in case no config is loaded!
            buffer.clear();
        }

        //<-- From PluginProcessor.cpp of MCFXConvolver
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

void tvconv_process(
    void* const hMcfxConv,
    float** const inputs,
    float** const outputs,
    int nInputs,
    int nOutputs,
    int nSamples) {

    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    int s, ch, i;
    int numInputChannels, numOutputChannels;

    tvconv_checkReInit(hMcfxConv);
    pData->procStatus = PROC_STATUS_ONGOING;

    numInputChannels = pData->nInputChannels;
    numOutputChannels = pData->nOutputChannels;

    // for(s=0; s<nSamples; s++){
    //     /* Load input signals into inFIFO buffer */
    //     for(ch=0; ch<SAF_MIN(SAF_MIN(nInputs,numInputChannels),MAX_NUM_CHANNELS); ch++)
    //         pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
    //     for(; ch<numInputChannels; ch++) /* Zero any channels that were not given */
    //         pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

    //     /* Pull output signals from outFIFO buffer */
    //     for(ch=0; ch<SAF_MIN(SAF_MIN(nOutputs, numOutputChannels),MAX_NUM_CHANNELS); ch++)
    //         outputs[ch][s] = pData->outFIFO[ch][pData->FIFO_idx];
    //     for(; ch<nOutputs; ch++) /* Zero any extra channels */
    //         outputs[ch][s] = 0.0f;

    //     /* Increment buffer index */
    //     pData->FIFO_idx++;

    //     /* Process frame if inFIFO is full and filters are loaded and saf_matrixConv_apply is ready for it */
    //     if (pData->FIFO_idx >= pData->hostBlockSize_clamped && pData->reInitFilters == 0 &&
    //         pData->codecStatus == CODEC_STATUS_INITIALISED) {
    //         pData->FIFO_idx = 0;

    //         /* Load time-domain data */
    //         for(i=0; i < numInputChannels; i++)
    //             utility_svvcopy(pData->inFIFO[i], pData->hostBlockSize_clamped, pData->inputFrameTD[i]);

    //         if(pData->hMCFXConv != NULL && pData->ir_length>0){
    //          mcfxConv_apply(pData->hMCFXConv,
    //                           FLATTEN2D(pData->inputFrameTD),
    //                           FLATTEN2D(pData->outputFrameTD),
    //                           pData->position_idx);
    //         }
    //         /* if the matrix convolver handle has not been initialised yet (i.e. no filters have been loaded) then zero the output */
    //         else
    //             memset(FLATTEN2D(pData->outputFrameTD), 0, MAX_NUM_CHANNELS * (pData->hostBlockSize_clamped)*sizeof(float));

    //         /* copy signals to output buffer */
    //         for (i = 0; i < SAF_MIN(numOutputChannels, MAX_NUM_CHANNELS); i++)
    //             utility_svvcopy(pData->outputFrameTD[i], pData->hostBlockSize_clamped, pData->outFIFO[i]);
    //     }
    //     else if(pData->FIFO_idx >= pData->hostBlockSize_clamped){
    //         /* clear outFIFO if codec was not ready */
    //         pData->FIFO_idx = 0;
    //         memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    //     }
    // }
    // pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/*sets*/

void tvconv_refreshParams(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    pData->reInitFilters = 1;
}

void tvconv_checkReInit(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);

    while (pData->procStatus == CODEC_STATUS_INITIALISING){
        SAF_SLEEP(10);
    }
    /* reinitialise if needed */
    if ((pData->reInitFilters == 1) && (pData->irs != NULL)) {
        pData->reInitFilters = 2;
        mcfxConv_destroy(&(pData->hMCFXConv));
        pData->hMCFXConv = NULL;
        /* if length of the loaded sofa file was not divisable by the specified number of inputs, then the handle remains NULL,
         * and no convolution is applied */
        // pData->hostBlockSize_clamped = SAF_CLAMP(pData->_BufferSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        if(pData->ir_length>0){
            bool checkAndNotifyResampling = false;
            if (pData->ir_fs != (int)pData->_SampleRate)
                checkAndNotifyResampling = true;
            mcfxConv_create(&(pData->hMCFXConv),
                              pData->_BufferSize, // This replaced pData->hostBlockSize_clamped,
                              pData->_ConvBufferSize, // PASSED AS REFERENCE
                              pData->irs,
                              pData->ir_length,
                              pData->nListenerPositions,
                              pData->nOutputChannels,
                              pData->position_idx,
                              pData->_SampleRate, // _sampleRate was passed as required by MCFX convdata, as double
                              pData->ir_fs,       // Original IR samplerate for resampling
                              pData->_MaxPartSize // PASSED AS REFERENCE
                             #if MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE //[ds 2024]
                              ,pData->num_xfd_convolvers
                             #endif
                              ); 
            if (checkAndNotifyResampling) {
                pData->resampled_ir = true;
            }
        
        }

        /* Resize buffers */
        pData->inputFrameTD  = (float**)realloc2d((void**)pData->inputFrameTD, 
                                                    MAX_NUM_CHANNELS,
                                                    pData->_BufferSize, /* This replaced pData->hostBlockSize_clamped,*/
                                                    sizeof(float));
        pData->outputFrameTD = (float**)realloc2d((void**)pData->outputFrameTD, 
                                                    MAX_NUM_CHANNELS,
                                                    pData->_BufferSize, /* This replaced pData->hostBlockSize_clamped,*/
                                                    sizeof(float));
        memset(FLATTEN2D(pData->inputFrameTD),0,MAX_NUM_CHANNELS*(pData->_BufferSize)*sizeof(float)); /*pData->_BufferSize replaced pData->hostBlockSize_clamped,pData->hostBlockSize_clamped*/

        /* reset FIFO buffers */
        pData->FIFO_idx = 0;
        memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
        memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));

        pData->reInitFilters = 0;
        pData->codecStatus = CODEC_STATUS_INITIALISED;
    }
}

void tvconv_setFiltersAndPositions(void* const hMcfxConv) {
        McfxConvData* pData = (McfxConvData*) hMcfxConv;
    #ifdef SAF_ENABLE_SOFA_READER_MODULE
        int i;
        vectorND tmp;
        SAF_SOFA_ERROR_CODES error;
        saf_sofa_container sofa;
    #endif

        if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
            return; /* re-init not required, or already happening */
        while (pData->procStatus == PROC_STATUS_ONGOING){
            /* re-init required, but we need to wait for the current processing loop to end */
            pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
            SAF_SLEEP(10);
        }

        /* for progress bar */
        pData->codecStatus = CODEC_STATUS_INITIALISING;
        strcpy(pData->progressBarText,"Initialising");
        pData->progressBar0_1 = 0.0f;

    #ifdef SAF_ENABLE_SOFA_READER_MODULE
        if(pData->sofa_filepath!=NULL){
            strcpy(pData->progressBarText,"Opening SOFA file");
            pData->progressBar0_1 = 0.2f;
            error = saf_sofa_open(&sofa, pData->sofa_filepath, SAF_SOFA_READER_OPTION_NETCDF);

            if(error==SAF_SOFA_OK){
                strcpy(pData->progressBarText,"Loading IRs");
                pData->progressBar0_1 = 0.5f;

                pData->ir_fs = (int)sofa.DataSamplingRate;
                pData->ir_length = sofa.DataLengthIR;
                pData->nIrChannels = sofa.nReceivers;
                pData->nListenerPositions = sofa.nListeners;

                /* copy only the first source position, because number of source positions might be incorrect in sofa */
                if(!strcmp(sofa.SourcePositionType, "spherical")){ // TODO: Change this to multiple sources
                    memcpy(tmp, sofa.SourcePosition, sizeof(vectorND));
                    unitSph2cart((float*)tmp, 1, 1, pData->sourcePosition);
                }
                else
                    memcpy(pData->sourcePosition, sofa.SourcePosition, sizeof(vectorND));

                pData->irs = (float**)realloc2d((void**)pData->irs, pData->nListenerPositions, pData->nIrChannels*pData->ir_length, sizeof(float));
                int tmp_length = pData->nIrChannels * pData->ir_length;
                for(i=0; i<pData->nListenerPositions; i++){
                    memcpy(pData->irs[i], &(sofa.DataIR[i*tmp_length]), tmp_length*sizeof(float));
                }

                strcpy(pData->progressBarText,"Loading positions");
                pData->progressBar0_1 = 0.8f;

                pData->listenerPositions = (vectorND*)realloc1d((void*)pData->listenerPositions, pData->nListenerPositions*sizeof(vectorND));
                memcpy(pData->listenerPositions, sofa.ListenerPosition, pData->nListenerPositions*sizeof(vectorND));
            }

            /* extra error handling */
            switch(error){
                case SAF_SOFA_OK: /** None of the error checks failed */
                    pData->sofa_file_error = SAF_TVCONV_SOFA_OK;
                    break;
                case SAF_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH:  /** Not a SOFA file, or no such file was found in the specified location */
                    pData->sofa_file_error = SAF_TVCONV_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH;
                    break;
                case SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED:      /** Dimensions of the SOFA data were not as expected */
                    pData->sofa_file_error = SAF_TVCONV_SOFA_ERROR_DIMENSIONS_UNEXPECTED;
                    break;
                case SAF_SOFA_ERROR_FORMAT_UNEXPECTED: /** The data-type of the SOFA data was not as expected */
                    pData->sofa_file_error = SAF_TVCONV_SOFA_ERROR_FORMAT_UNEXPECTED;
                    break;
                case SAF_SOFA_ERROR_NETCDF_IN_USE: /** NetCDF is not thread safe! */
                    pData->sofa_file_error = SAF_TVCONV_SOFA_ERROR_NETCDF_IN_USE;
                    break;
            }
        }
        saf_sofa_close(&sofa);
        pData->nOutputChannels = SAF_MIN(pData->nIrChannels, MAX_NUM_CHANNELS);
        internal_setMinMaxDimensions(hMcfxConv);
    #else
        pData->ir_length = 0;
        saf_print_warning("This example requires SAF_ENABLE_SOFA_READER_MODULE to do anything");
    #endif

        pData->position_idx = 0;
        pData->codecStatus = CODEC_STATUS_INITIALISED;
        pData->reInitFilters = 1;

        /* done! */
        strcpy(pData->progressBarText,"Done!");
        pData->progressBar0_1 = 1.0f;
}

void tvconv_reinitConvolver(void* const hMcfxConv) {
    if (hMcfxConv == NULL) return;
    ((McfxConvData*)(hMcfxConv))->reInitFilters = 1;
}

void tvconv_setSofaFilePath(void* const hMcfxConv, const char* path) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    pData->sofa_file_error = SAF_TVCONV_SOFA_LOADING;
    pData->sofa_filepath = (char*)malloc1d(strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    tvconv_setFiltersAndPositions(hMcfxConv);
}

void tvconv_setTargetPosition(void* const hMcfxConv, float position, int dim) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    pData->targetPosition[dim] = position;
    internal_findNearestNeigbour(hMcfxConv);
}

/*gets*/

int tvconv_getNumInputChannels(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->nInputChannels;
}

int tvconv_getNumOutputChannels(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->nOutputChannels;
}

int tvconv_getHostBlockSize(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->_BufferSize; // This replaced pData->hostBlockSize
}

int tvconv_getNumIRs(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->nIrChannels;
}

int tvconv_getNumListenerPositions(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->codecStatus==CODEC_STATUS_INITIALISED ? pData->nListenerPositions : 0;
}

float tvconv_getListenerPosition(void* const hMcfxConv, int index, int dim) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->codecStatus==CODEC_STATUS_INITIALISED ? pData->listenerPositions[index][dim] : 0.0f;
}

int tvconv_getListenerPositionIdx(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->position_idx;
}

float tvconv_getTargetPosition(void* const hMcfxConv, int dim) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->targetPosition[dim];
}

float tvconv_getSourcePosition(void* const hMcfxConv, int dim) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->sourcePosition[dim];
}

float tvconv_getMinDimension(void* const hMcfxConv, int dim) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->minDimensions[dim];
}

float tvconv_getMaxDimension(void* const hMcfxConv, int dim) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->maxDimensions[dim];
}

int tvconv_getIRLength(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->ir_length;
}

int tvconv_getIRFs(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->ir_fs;
}

int tvconv_getResampledIR(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->resampled_ir;
}

int tvconv_getHostFs(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->_SampleRate; // pData->_SampleRate replaced saf's original host_fs
}

int tvconv_getProcessingDelay(void* const hMcfxConv) {
    McfxConvData* pData = (McfxConvData*)(hMcfxConv);
    if (pData == NULL) return 0; // In the unlikely event that the handle is invalid, return 0 as the plugin latency in samples

    Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv);
    if (h == NULL)  return 0; // If no filters are loaded, return 0 as the plugin latency in samples

    return h->latencySamples;
}

char* tvconv_getSofaFilePath(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    if(pData->sofa_filepath!=NULL)
        return pData->sofa_filepath;
    else
        return "no_file";
}

SAF_TVCONV_ERROR_CODES tvconv_getSofaErrorState(void* const hMcfxConv)
{
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->sofa_file_error;
}

CODEC_STATUS tvconv_getCodecStatus(void* const hMcfxConv) {
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->codecStatus;
}


// --- internal functions --- //
// From tvconv_internal.h/.cpp

/** Sets codec status (see #CODEC_STATUS enum) */
void internal_setCodecStatus(void* const hTVCnv, CODEC_STATUS newStatus)
{
    McfxConvData *pData = (McfxConvData*)(hTVCnv);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

/** Finds the index holding the nearest neigbour to the selected position */
void internal_findNearestNeigbour(void* const hTVCnv)
{
    float dist = 0, minDist = 0;
    int i, d, min_idx = 0;
    McfxConvData *pData = (McfxConvData*)(hTVCnv);
    if (pData->nListenerPositions > 0 && pData->listenerPositions != NULL) {
        for(i = 0; i < pData->nListenerPositions; i++){
            for(d = 0; d < NUM_DIMENSIONS; d++)
                dist += (pData->targetPosition[d] - pData->listenerPositions[i][d]) *
                        (pData->targetPosition[d] - pData->listenerPositions[i][d]);
            
            if(dist < minDist || i == 0){
                minDist = dist;
                min_idx = i;
            }
            dist = 0;
        }
    }

    pData->position_idx = min_idx;
}

/**
 * Sets the smallest and the highest position of each dimension from the list of
 * positions
 */
void internal_setMinMaxDimensions(void* const hTVCnv)
{
    int i, d;
    McfxConvData *pData = (McfxConvData*)(hTVCnv);
    if(pData->listenerPositions != NULL){
        for(d = 0; d < NUM_DIMENSIONS; d++){
            pData->minDimensions[d] = pData->listenerPositions[0][d];
            pData->maxDimensions[d] = pData->listenerPositions[0][d];
            for(i = 1; i<pData->nListenerPositions; i++){
                if(pData->listenerPositions[i][d] < pData->minDimensions[d])
                    pData->minDimensions[d] = pData->listenerPositions[i][d];
                else if (pData->listenerPositions[i][d] > pData->maxDimensions[d])
                    pData->maxDimensions[d] = pData->listenerPositions[i][d];
            }
        }
        /* resetting current position to the start */
        for(d = 0; d < NUM_DIMENSIONS; d++)
            pData->targetPosition[d] = pData->minDimensions[d];
    }
}




/** From MCFX Convolver PluginProcessor.cpp */
void getIntFromLine(int& ret, String& line) {
    if (line.isEmpty()) return;
    ret = line.getIntValue();
    line = line.fromFirstOccurrenceOf(" ", false, true).trim();
}

/** From MCFX Convolver PluginProcessor.cpp */
void getFloatFromLine(float& ret, String& line) {
    if (line.isEmpty()) return;
    ret = line.getFloatValue();
    line = line.fromFirstOccurrenceOf(" ", false, true).trim();
}

//TODO: Remove this if it goes unused
//void unloadConfiguration(void* const hMcfxConv) {
//    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
//    Internal_Conv_struct *h = (Internal_Conv_struct*)(pData->hMCFXConv);
//    MtxConvMaster* curMtxConv = h->mtxconv_s.getUnchecked(pData->position_idx);
//    ConvolverData* currentConvData = h->convdata_s.getUnchecked(pData->position_idx);
//    // TODO: possibly need to check for h to be not NULL
//    
//    // delete configuration
//    h->_configLoaded = false;
//    pData->_conv_in.clear();
//    pData->_conv_out.clear();
//    pData->_min_in_ch = 0;
//    pData->_min_out_ch = 0;
//
//#ifdef USE_ZITA_CONVOLVER //TODO: fix ZITA with new multi convolver paradigm
//    h->zita_conv->stop_process();
//    h->zita_conv->cleanup();
//#else
//    curMtxConv->StopProc();
//    curMtxConv->Cleanup();
//#endif
//
//    currentConvData->clear();
//    std::cout << "Unloaded Convolution..." << std::endl;
//}

int mcfxConv_getSkippedCyclesCount(void* const hMcfxConv) //TODO: use in the plugin
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);

    if (pData->hMCFXConv == NULL) return 0; // If no filters are loaded, return 0 as the number of skipped cycles
    Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv); 

    return h->_skippedCycles.get();
}

int mcfxConv_getMinInCh(void* const hMcfxConv) //TODO: use in the plugin
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);

    if (pData->hMCFXConv == NULL) return 0; // If no filters are loaded, return 0 as the number of skipped cycles
    Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv); 

    return h->_min_in_ch;
}

int mcfxConv_getMinOutCh(void* const hMcfxConv) //TODO: use in the plugin
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);

    if (pData->hMCFXConv == NULL) return 0; // If no filters are loaded, return 0 as the number of skipped cycles
    Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv); 

    return h->_min_out_ch;
}

int mcfxConv_getNumConv(void* const hMcfxConv) //TODO: use in the plugin
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);

    if (pData->hMCFXConv == NULL) return 0; // If no filters are loaded, return 0 as the number of skipped cycles
    Internal_Conv_struct* h = (Internal_Conv_struct*)(pData->hMCFXConv); 

    return h->_num_conv;
}

unsigned int mcfxConv_getMaxPartitionSize(void* const hMcfxConv)
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->_MaxPartSize;
}



unsigned int mcfxConv_getBufferSize(void* const hMcfxConv)
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->_BufferSize;
}

unsigned int mcfxConv_getConvBufferSize(void* const hMcfxConv)
{
    if (hMcfxConv == NULL) return 0; // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);
    return pData->_ConvBufferSize;
}

void mcfxConv_setConvBufferSize(void* const hMcfxConv, unsigned int convBufferSize) {
    McfxConvData* pData = (McfxConvData*)(hMcfxConv);

    if (nextPowerOfTwo(convBufferSize) != pData->_ConvBufferSize) {
        pData->_ConvBufferSize = nextPowerOfTwo(convBufferSize);
        tvconv_reinitConvolver(hMcfxConv);
    }
}

void mcfxConv_setMaxPartitionSize(void* const hMcfxConv, unsigned int maxPartitionSize) {
    if (hMcfxConv == NULL) return;  // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData* pData = (McfxConvData*)(hMcfxConv);

    if (maxPartitionSize > MAX_PART_SIZE)
        return;

    if (nextPowerOfTwo(maxPartitionSize) != pData->_MaxPartSize) {
        pData->_MaxPartSize = nextPowerOfTwo(maxPartitionSize);
        tvconv_reinitConvolver(hMcfxConv);
    }
}

#if MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE //[ds 2024]
float mcfxConv_getMaxCrossfadeTimeS(void* const hMcfxConv, bool* minReached, bool* maxReached) {
    if (hMcfxConv == NULL) return 0; 
    McfxConvData *pData = (McfxConvData*)(hMcfxConv);

    //TODO: return 0 if mtxconv is NULL
    if (minReached)
        *minReached = (pData->num_xfd_convolvers <= 1 || pData->num_xfd_convolvers / 2 < 1);
    if (maxReached)
        *maxReached = (pData->num_xfd_convolvers >= MAX_NCONVOLVERS || pData->num_xfd_convolvers *2 > MAX_NCONVOLVERS);

    // Maximum crossfade time in seconds is the bufferSize * pData->numConvolvers-1 divided by the sample rate
    float maxCrossfadeTime = (float)(pData->_BufferSize * (pData->num_xfd_convolvers-1)) / (float)pData->_SampleRate;
    return maxCrossfadeTime;
}

void mcfxConv_DoubleCrossfadeTime(void* const hMcfxConv, bool* maxReached)
{
    if (hMcfxConv == NULL) return;  // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData* pData = (McfxConvData*)(hMcfxConv);

    if (pData->num_xfd_convolvers >= MAX_NCONVOLVERS || pData->num_xfd_convolvers *2 > MAX_NCONVOLVERS) {
        *maxReached = true;
        return; 
    }
    *maxReached = false;

    pData->num_xfd_convolvers *= 2;
    tvconv_reinitConvolver(hMcfxConv);
}

void mcfxConv_HalveCrossfadeTime(void* const hMcfxConv, bool* minReached)
{
    if (hMcfxConv == NULL) return;  // If the handle is invalid, return 0 as the number of skipped cycles
    McfxConvData* pData = (McfxConvData*)(hMcfxConv);

    if (pData->num_xfd_convolvers <= 1 || pData->num_xfd_convolvers / 2 < 1) {
        *minReached = true;
        return; 
    }
    *minReached = false;

    pData->num_xfd_convolvers /= 2;
    tvconv_reinitConvolver(hMcfxConv);
}
#endif