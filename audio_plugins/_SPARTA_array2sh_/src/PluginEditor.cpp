/*
  ==============================================================================

  This is an automatically generated GUI class created by the Projucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Projucer version: 5.3.2

  ------------------------------------------------------------------------------

  The Projucer is part of the JUCE library.
  Copyright (c) 2017 - ROLI Ltd.

  ==============================================================================
*/

//[Headers] You can add your own extra header files here...

//[/Headers]

#include "PluginEditor.h"


//[MiscUserDefs] You can add your own user definitions and misc code here...

//[/MiscUserDefs]

//==============================================================================
PluginEditor::PluginEditor (PluginProcessor* ownerFilter)
    : AudioProcessorEditor(ownerFilter)
{
    //[Constructor_pre] You can add your own custom stuff here..
    //[/Constructor_pre]

    presetCB.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (presetCB.get());
    presetCB->setEditableText (false);
    presetCB->setJustificationType (Justification::centredLeft);
    presetCB->setTextWhenNothingSelected (String());
    presetCB->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    presetCB->addItem (TRANS("Default"), 1);
    presetCB->addListener (this);

    presetCB->setBounds (88, 64, 120, 16);

    arrayTypeCB.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (arrayTypeCB.get());
    arrayTypeCB->setEditableText (false);
    arrayTypeCB->setJustificationType (Justification::centredLeft);
    arrayTypeCB->setTextWhenNothingSelected (TRANS("Spherical"));
    arrayTypeCB->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    arrayTypeCB->addItem (TRANS("Spherical"), 1);
    arrayTypeCB->addItem (TRANS("Cylindrical"), 2);
    arrayTypeCB->addListener (this);

    arrayTypeCB->setBounds (360, 344, 120, 16);

    QSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (QSlider.get());
    QSlider->setRange (4, 64, 1);
    QSlider->setSliderStyle (Slider::LinearHorizontal);
    QSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    QSlider->addListener (this);

    QSlider->setBounds (112, 96, 96, 16);

    rSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (rSlider.get());
    rSlider->setRange (0.01, 0.3, 0.001);
    rSlider->setSliderStyle (Slider::LinearHorizontal);
    rSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    rSlider->addListener (this);

    rSlider->setBounds (112, 128, 96, 16);

    RSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (RSlider.get());
    RSlider->setRange (0.01, 0.3, 0.001);
    RSlider->setSliderStyle (Slider::LinearHorizontal);
    RSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    RSlider->addListener (this);

    RSlider->setBounds (112, 160, 96, 16);

    cSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (cSlider.get());
    cSlider->setRange (200, 2000, 0.1);
    cSlider->setSliderStyle (Slider::LinearHorizontal);
    cSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    cSlider->addListener (this);

    cSlider->setBounds (360, 313, 120, 16);

    weightTypeCB.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (weightTypeCB.get());
    weightTypeCB->setEditableText (false);
    weightTypeCB->setJustificationType (Justification::centredLeft);
    weightTypeCB->setTextWhenNothingSelected (TRANS("Rigid"));
    weightTypeCB->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    weightTypeCB->addItem (TRANS("Rigid"), 1);
    weightTypeCB->addItem (TRANS("Open_omni"), 2);
    weightTypeCB->addItem (TRANS("Open_card"), 3);
    weightTypeCB->addItem (TRANS("Open_dipole"), 4);
    weightTypeCB->addListener (this);

    weightTypeCB->setBounds (360, 377, 120, 16);

    regTypeCB.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (regTypeCB.get());
    regTypeCB->setEditableText (false);
    regTypeCB->setJustificationType (Justification::centredLeft);
    regTypeCB->setTextWhenNothingSelected (TRANS("Soft Limiting"));
    regTypeCB->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    regTypeCB->addItem (TRANS("DaS"), 1);
    regTypeCB->addItem (TRANS("Soft Limiting"), 2);
    regTypeCB->addItem (TRANS("Tikhonov"), 3);
    regTypeCB->addListener (this);

    regTypeCB->setBounds (640, 280, 128, 16);

    regAmountSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (regAmountSlider.get());
    regAmountSlider->setRange (0, 80, 0.01);
    regAmountSlider->setSliderStyle (Slider::LinearHorizontal);
    regAmountSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    regAmountSlider->addListener (this);

    regAmountSlider->setBounds (640, 312, 128, 16);

    CHOrderingCB.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (CHOrderingCB.get());
    CHOrderingCB->setEditableText (false);
    CHOrderingCB->setJustificationType (Justification::centredLeft);
    CHOrderingCB->setTextWhenNothingSelected (TRANS("ACN"));
    CHOrderingCB->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    CHOrderingCB->addItem (TRANS("ACN"), 1);
    CHOrderingCB->addListener (this);

    CHOrderingCB->setBounds (640, 376, 128, 16);

    normalisationCB.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (normalisationCB.get());
    normalisationCB->setEditableText (false);
    normalisationCB->setJustificationType (Justification::centredLeft);
    normalisationCB->setTextWhenNothingSelected (TRANS("N3D"));
    normalisationCB->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    normalisationCB->addItem (TRANS("N3D"), 1);
    normalisationCB->addItem (TRANS("SN3D"), 2);
    normalisationCB->addListener (this);

    normalisationCB->setBounds (640, 408, 128, 16);

    maxFreqSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (maxFreqSlider.get());
    maxFreqSlider->setRange (5000, 24000, 1);
    maxFreqSlider->setSliderStyle (Slider::LinearHorizontal);
    maxFreqSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    maxFreqSlider->addListener (this);

    maxFreqSlider->setBounds (360, 408, 120, 16);

    gainSlider.reset (new Slider ("new slider"));
    addAndMakeVisible (gainSlider.get());
    gainSlider->setRange (-60, 60, 0.01);
    gainSlider->setSliderStyle (Slider::LinearHorizontal);
    gainSlider->setTextBoxStyle (Slider::TextBoxRight, false, 55, 20);
    gainSlider->addListener (this);

    gainSlider->setBounds (640, 344, 128, 16);

    degRadTB.reset (new ToggleButton ("new toggle button"));
    addAndMakeVisible (degRadTB.get());
    degRadTB->setButtonText (String());
    degRadTB->addListener (this);

    degRadTB->setBounds (186, 198, 23, 24);

    textButton.reset (new TextButton ("new button"));
    addAndMakeVisible (textButton.get());
    textButton->setButtonText (TRANS("Analyse"));
    textButton->addListener (this);
    textButton->setColour (TextButton::buttonColourId, Colour (0xff5c68a4));

    textButton->setBounds (712, 37, 72, 16);

    dispWindow.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (dispWindow.get());
    dispWindow->setEditableText (false);
    dispWindow->setJustificationType (Justification::centredLeft);
    dispWindow->setTextWhenNothingSelected (TRANS("EQ"));
    dispWindow->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    dispWindow->addItem (TRANS("EQ"), 1);
    dispWindow->addItem (TRANS("Coh"), 2);
    dispWindow->addItem (TRANS("L Diff"), 3);
    dispWindow->addListener (this);

    dispWindow->setBounds (642, 37, 63, 16);

    tb_loadJSON.reset (new TextButton ("new button"));
    addAndMakeVisible (tb_loadJSON.get());
    tb_loadJSON->setButtonText (TRANS("Import"));
    tb_loadJSON->setConnectedEdges (Button::ConnectedOnRight);
    tb_loadJSON->addListener (this);
    tb_loadJSON->setColour (TextButton::buttonColourId, Colour (0xff14889e));

    tb_loadJSON->setBounds (148, 39, 34, 14);

    tb_saveJSON.reset (new TextButton ("new button"));
    addAndMakeVisible (tb_saveJSON.get());
    tb_saveJSON->setButtonText (TRANS("Export"));
    tb_saveJSON->setConnectedEdges (Button::ConnectedOnLeft);
    tb_saveJSON->addListener (this);
    tb_saveJSON->setColour (TextButton::buttonColourId, Colour (0xff224d97));
    tb_saveJSON->setColour (TextButton::buttonOnColourId, Colour (0xff181f9a));

    tb_saveJSON->setBounds (182, 39, 34, 14);

    CBencodingOrder.reset (new ComboBox ("new combo box"));
    addAndMakeVisible (CBencodingOrder.get());
    CBencodingOrder->setEditableText (false);
    CBencodingOrder->setJustificationType (Justification::centredLeft);
    CBencodingOrder->setTextWhenNothingSelected (TRANS("Default"));
    CBencodingOrder->setTextWhenNoChoicesAvailable (TRANS("(no choices)"));
    CBencodingOrder->addListener (this);

    CBencodingOrder->setBounds (363, 276, 118, 20);


    //[UserPreSize]
    //[/UserPreSize]

    setSize (800, 450);


    //[Constructor] You can add your own custom stuff here..
	hVst = ownerFilter;
    openGLContext.setMultisamplingEnabled(true);
    openGLContext.attachTo(*this);

    /* create EQ window and analysis windows */
    addAndMakeVisible (eqviewIncluded = new eqview(556, 209, 100.0f, 20e3f, -40.0f, 100.0f, 48e3f )); /* TODO: switch to the more general "anaview"  */
    eqviewIncluded->setAlwaysOnTop(true);
    eqviewIncluded->setTopLeftPosition(228, 56);
    eqviewIncluded->setVisible(true);
    addAndMakeVisible (cohviewIncluded = new anaview(556, 209, 100.0f, 20e3f, 0.0f, 1.0f, TRANS("Spatial Corr. (T:I)"), 1, 48e3f ));
    cohviewIncluded->setAlwaysOnTop(true);
    cohviewIncluded->setTopLeftPosition(228, 56);
    cohviewIncluded->setVisible(false);
    addAndMakeVisible (ldiffviewIncluded = new anaview(556, 209, 100.0f, 20e3f, -30, 10, TRANS("Level Diff. (dB)"), 10.0f, 48e3f ));
    ldiffviewIncluded->setAlwaysOnTop(true);
    ldiffviewIncluded->setTopLeftPosition(228, 56);
    ldiffviewIncluded->setVisible(false);
    dispID = SHOW_EQ;
    needScreenRefreshFLAG = true;

    /* add master decoding order options */
    CBencodingOrder->addItem (TRANS("1st order"), ENCODING_ORDER_FIRST);
    CBencodingOrder->addItem (TRANS("2nd order"), ENCODING_ORDER_SECOND);
    CBencodingOrder->addItem (TRANS("3rd order"), ENCODING_ORDER_THIRD);
    CBencodingOrder->addItem (TRANS("4th order"), ENCODING_ORDER_FOURTH);
    CBencodingOrder->addItem (TRANS("5th order"), ENCODING_ORDER_FIFTH);
    CBencodingOrder->addItem (TRANS("6th order"), ENCODING_ORDER_SIXTH);
    CBencodingOrder->addItem (TRANS("7th order"), ENCODING_ORDER_SEVENTH);

    /* pass handles to data required for eq and analysis displays */
    int numFreqPoints, numCurves;
    float* freqVector = array2sh_getFreqVector(hVst->hA2sh, &numFreqPoints);
    float** solidCurves = array2sh_getbN_inv(hVst->hA2sh, &numCurves, &numFreqPoints);
    eqviewIncluded->setSolidCurves_Handle(freqVector, solidCurves, numFreqPoints, numCurves);
    float** faintCurves = array2sh_getbN_modal(hVst->hA2sh, &numCurves, &numFreqPoints);
    eqviewIncluded->setFaintCurves_Handle(freqVector, faintCurves, numFreqPoints, numCurves);
    float* dataHandle = array2sh_getSpatialCorrelation_Handle(hVst->hA2sh, &numCurves, &numFreqPoints);
    cohviewIncluded->setSolidCurves_Handle(freqVector, dataHandle, numFreqPoints, numCurves);
    dataHandle = array2sh_getLevelDifference_Handle(hVst->hA2sh, &numCurves, &numFreqPoints);
    ldiffviewIncluded->setSolidCurves_Handle(freqVector, dataHandle, numFreqPoints, numCurves);

    /* grab current parameter settings */
    CBencodingOrder->setSelectedId(array2sh_getEncodingOrder(hVst->hA2sh), dontSendNotification);
    arrayTypeCB->setSelectedId(array2sh_getArrayType(hVst->hA2sh), dontSendNotification);
    QSlider->setRange(array2sh_getMinNumSensors(hVst->hA2sh), array2sh_getMaxNumSensors(), 1);
    QSlider->setValue(array2sh_getNumSensors(hVst->hA2sh), dontSendNotification);
    rSlider->setValue(array2sh_getr(hVst->hA2sh), dontSendNotification);
    RSlider->setValue(array2sh_getR(hVst->hA2sh), dontSendNotification);
    cSlider->setValue(array2sh_getc(hVst->hA2sh), dontSendNotification);
    weightTypeCB->setSelectedId(array2sh_getWeightType(hVst->hA2sh), dontSendNotification);
    regTypeCB->setSelectedId(array2sh_getRegType(hVst->hA2sh), dontSendNotification);
    regAmountSlider->setValue(array2sh_getRegPar(hVst->hA2sh), dontSendNotification);
    CHOrderingCB->setSelectedId(array2sh_getChOrder(hVst->hA2sh), dontSendNotification);
    normalisationCB->setSelectedId(array2sh_getNormType(hVst->hA2sh), dontSendNotification);
    maxFreqSlider->setValue(array2sh_getMaxFreq(hVst->hA2sh), dontSendNotification);
    gainSlider->setValue(array2sh_getGain(hVst->hA2sh), dontSendNotification);
    showDegreesInstead = false;
    degRadTB->setToggleState(showDegreesInstead, dontSendNotification);

    /* Hide decoding orders that are unsuitable for number of sensors */
    for(int i=1; i<=7; i++)
        CBencodingOrder->setItemEnabled(i, (i+1)*(i+1) <= array2sh_getNumSensors(hVst->hA2sh) ? true : false);

    /* gain range should change with order */
    int curOrder = array2sh_getEncodingOrder(hVst->hA2sh);
    gainSlider->setRange (-10-(curOrder*curOrder), 5+curOrder, 0.01);

    /* sensor coord table */
    addAndMakeVisible (sensorCoordsVP = new Viewport ("new viewport"));
    sensorCoordsView_handle = new sensorCoordsView(ownerFilter, MAX_NUM_CHANNELS, array2sh_getNumSensors(hVst->hA2sh), showDegreesInstead);
    sensorCoordsVP->setViewedComponent (sensorCoordsView_handle);
    sensorCoordsVP->setScrollBarsShown (true, false);
    sensorCoordsVP->setAlwaysOnTop(true);
    sensorCoordsVP->setBounds(24, 224, 184, 200);
    sensorCoordsView_handle->setQ(array2sh_getNumSensors(hVst->hA2sh));

    /* disable unused parameters */
    bool shouldBeEnabled = array2sh_getWeightType(hVst->hA2sh) != 1 ? false : true; /* is it a rigid array? */
    RSlider->setEnabled(shouldBeEnabled);
    shouldBeEnabled = array2sh_getArrayType(hVst->hA2sh) != 1 ? false : true;  /* is it a cylindrical array? */
    weightTypeCB->setItemEnabled(3, shouldBeEnabled);
    weightTypeCB->setItemEnabled(4, shouldBeEnabled);
    if((array2sh_getWeightType(hVst->hA2sh) == 3) || (array2sh_getWeightType(hVst->hA2sh) == 4 ))
        weightTypeCB->setSelectedId(1, dontSendNotification);
    shouldBeEnabled = array2sh_getRegType(hVst->hA2sh) == 1 ? false : true;  /* is DaS regularisation selected? */
    regAmountSlider->setEnabled(shouldBeEnabled);
    regTypeCB->setItemEnabled(1, false); /* disable DaS */

    /* Presets */
    presetCB->setTextWhenNothingSelected (TRANS("Default"));
#ifdef ENABLE_SENNHEISER_AMBEO_PRESET
    presetCB->addItem (TRANS("Sennheiser Ambeo"), PRESET_SENNHEISER_AMBEO);
#endif
#ifdef ENABLE_CORE_SOUND_TETRAMIC_PRESET
    presetCB->addItem (TRANS("Core Sound TetraMic"), PRESET_CORE_SOUND_TETRAMIC);
#endif
#ifdef ENABLE_SOUND_FIELD_SPS200_PRESET
    presetCB->addItem (TRANS("Sound-field SPS200"), PRESET_SOUND_FIELD_SPS200);
#endif
#ifdef ENABLE_ZYLIA_1D_PRESET
    presetCB->addItem (TRANS("Zylia 1D"), PRESET_ZYLIA_1D);
#endif
#ifdef ENABLE_EIGENMIKE32_PRESET
    presetCB->addItem (TRANS("Eigenmike32"), PRESET_EIGENMIKE32);
#endif
#ifdef ENABLE_DTU_MIC_PRESET
    presetCB->addItem (TRANS("DTU mic"), PRESET_DTU_MIC );
#endif
#ifdef ENABLE_AALTO_HYDROPHONE_PRESET
    presetCB->addItem (TRANS("Aalto Hydro"), PRESET_AALTO_HYDROPHONE);
#endif
	/* Specify screen refresh rate */
    startTimer(80);//80); /*ms (40ms = 25 frames per second) */

    /* warnings */
    currentWarning = k_warning_none;

    //[/Constructor]
}

PluginEditor::~PluginEditor()
{
    //[Destructor_pre]. You can add your own custom destruction code here..
    //[/Destructor_pre]

    presetCB = nullptr;
    arrayTypeCB = nullptr;
    QSlider = nullptr;
    rSlider = nullptr;
    RSlider = nullptr;
    cSlider = nullptr;
    weightTypeCB = nullptr;
    regTypeCB = nullptr;
    regAmountSlider = nullptr;
    CHOrderingCB = nullptr;
    normalisationCB = nullptr;
    maxFreqSlider = nullptr;
    gainSlider = nullptr;
    degRadTB = nullptr;
    textButton = nullptr;
    dispWindow = nullptr;
    tb_loadJSON = nullptr;
    tb_saveJSON = nullptr;
    CBencodingOrder = nullptr;


    //[Destructor]. You can add your own custom destruction code here..
    eqviewIncluded = nullptr;
    sensorCoordsVP = nullptr;
    sensorCoordsView_handle = nullptr;
    //[/Destructor]
}

//==============================================================================
void PluginEditor::paint (Graphics& g)
{
    //[UserPrePaint] Add your own custom painting code here..
    //[/UserPrePaint]

    g.fillAll (Colours::white);

    {
        int x = 0, y = 30, width = 800, height = 420;
        Colour fillColour1 = Colour (0xff55636d), fillColour2 = Colour (0xff073642);
        Colour strokeColour = Colour (0xffa3a4a5);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setGradientFill (ColourGradient (fillColour1,
                                       408.0f - 0.0f + x,
                                       224.0f - 30.0f + y,
                                       fillColour2,
                                       688.0f - 0.0f + x,
                                       392.0f - 30.0f + y,
                                       true));
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 2);

    }

    {
        int x = 500, y = 264, width = 284, height = 172;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 12, y = 87, width = 204, height = 106;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 0, y = 0, width = 800, height = 32;
        Colour fillColour = Colour (0xff073642);
        Colour strokeColour = Colour (0xdcbdbdbd);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 3);

    }

    {
        int x = -19, y = 0, width = 195, height = 32;
        String text (TRANS("Array2SH"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (18.80f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centred, true);
    }

    {
        int x = 20, y = 57, width = 67, height = 30;
        String text (TRANS("Presets: "));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 20, y = 89, width = 92, height = 30;
        String text (TRANS("No. Sensors: "));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 20, y = 121, width = 116, height = 30;
        String text (TRANS("Array r (m):"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 20, y = 153, width = 124, height = 30;
        String text (TRANS("Sensor r (m):"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 12, y = 192, width = 204, height = 244;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 228, y = 264, width = 273, height = 172;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 244, y = 305, width = 172, height = 30;
        String text (TRANS("c (m/s):"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 12, y = 56, width = 204, height = 32;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 244, y = 337, width = 172, height = 30;
        String text (TRANS("Array Type:"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 244, y = 369, width = 172, height = 30;
        String text (TRANS("Weight Type:"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 524, y = 273, width = 172, height = 30;
        String text (TRANS("Reg. Type:"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 524, y = 305, width = 172, height = 30;
        String text (TRANS("Max Gain (dB):"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 523, y = 369, width = 172, height = 30;
        String text (TRANS("CH Order:"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 523, y = 401, width = 172, height = 30;
        String text (TRANS("Normalisation:"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 228, y = 56, width = 556, height = 209;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 523, y = 337, width = 172, height = 30;
        String text (TRANS("Post Gain (dB):"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 243, y = 401, width = 172, height = 30;
        String text (TRANS("Max Freq. (Hz):"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 75, y = 193, width = 85, height = 30;
        String text (TRANS("Azi   #   Elev"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 24, y = 224, width = 184, height = 200;
        Colour strokeColour = Colour (0x29b6b5b5);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 179, y = 193, width = 24, height = 23;
        String text (CharPointer_UTF8 ("\xc2\xb0"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 92, y = 30, width = 88, height = 30;
        String text (TRANS("Inputs"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 440, y = 30, width = 149, height = 30;
        String text (TRANS("Encoding Settings"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    {
        int x = 228, y = 264, width = 273, height = 41;
        Colour fillColour = Colour (0x13f4f4f4);
        Colour strokeColour = Colour (0x67a0a0a0);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.fillRect (x, y, width, height);
        g.setColour (strokeColour);
        g.drawRect (x, y, width, height, 1);

    }

    {
        int x = 243, y = 270, width = 172, height = 30;
        String text (TRANS("Encoding Order:"));
        Colour fillColour = Colours::white;
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setColour (fillColour);
        g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Bold"));
        g.drawText (text, x, y, width, height,
                    Justification::centredLeft, true);
    }

    //[UserPaint] Add your own custom painting code here..

	g.setColour(Colours::white);
	g.setFont(Font(11.00f, Font::plain));
    g.drawText(TRANS("Ver ") + JucePlugin_VersionString + BUILD_VER_SUFFIX + TRANS(", Build Date ") + __DATE__ + TRANS(" "),
		150, 16, 530, 11,
		Justification::centredLeft, true);

    /* display warning message */
    g.setColour(Colours::red);
    g.setFont(Font(11.00f, Font::plain));
    switch (currentWarning){
        case k_warning_none:
            break;
        case k_warning_frameSize:
            g.drawText(TRANS("Set frame size to multiple of ") + String(FRAME_SIZE),
                       getBounds().getWidth()-225, 16, 530, 11,
                       Justification::centredLeft, true);
            break;
        case k_warning_NinputCH:
            g.drawText(TRANS("Insufficient number of input channels (") + String(hVst->getTotalNumInputChannels()) +
                       TRANS("/") + String(array2sh_getNumSensors(hVst->hA2sh)) + TRANS(")"),
                       getBounds().getWidth()-225, 16, 530, 11,
                       Justification::centredLeft, true);
            break;
        case k_warning_NoutputCH:
            g.drawText(TRANS("Insufficient number of output channels (") + String(hVst->getTotalNumOutputChannels()) +
                       TRANS("/") + String(array2sh_getNSHrequired(hVst->hA2sh)) + TRANS(")"),
                       getBounds().getWidth()-225, 16, 530, 11,
                       Justification::centredLeft, true);
            break;
    }

    //[/UserPaint]
}

void PluginEditor::resized()
{
    //[UserPreResize] Add your own custom resize code here..
    //[/UserPreResize]

    //[UserResized] Add your own custom resize handling here..

    //[/UserResized]
}

void PluginEditor::comboBoxChanged (ComboBox* comboBoxThatHasChanged)
{
    //[UsercomboBoxChanged_Pre]

    //[/UsercomboBoxChanged_Pre]

    if (comboBoxThatHasChanged == presetCB.get())
    {
        //[UserComboBoxCode_presetCB] -- add your combo box handling code here..
        array2sh_setPreset(hVst->hA2sh, (int)presetCB->getSelectedId());

        /* grab current parameter settings */
        arrayTypeCB->setSelectedId(array2sh_getArrayType(hVst->hA2sh), dontSendNotification);
        CBencodingOrder->setSelectedId(array2sh_getEncodingOrder(hVst->hA2sh), dontSendNotification);
        int curOrder = CBencodingOrder->getSelectedId();
        QSlider->setRange((curOrder+1)*(curOrder+1), array2sh_getMaxNumSensors(), 1);
        QSlider->setValue(array2sh_getNumSensors(hVst->hA2sh), dontSendNotification);
        rSlider->setValue(array2sh_getr(hVst->hA2sh), dontSendNotification);
        RSlider->setValue(array2sh_getR(hVst->hA2sh), dontSendNotification);
        cSlider->setValue(array2sh_getc(hVst->hA2sh), dontSendNotification);
        weightTypeCB->setSelectedId(array2sh_getWeightType(hVst->hA2sh), dontSendNotification);
        regTypeCB->setSelectedId(array2sh_getRegType(hVst->hA2sh), dontSendNotification);
        regAmountSlider->setValue(array2sh_getRegPar(hVst->hA2sh), dontSendNotification);
        CHOrderingCB->setSelectedId(array2sh_getChOrder(hVst->hA2sh), dontSendNotification);
        normalisationCB->setSelectedId(array2sh_getNormType(hVst->hA2sh), dontSendNotification);
        maxFreqSlider->setValue(array2sh_getMaxFreq(hVst->hA2sh), dontSendNotification);
        gainSlider->setValue(array2sh_getGain(hVst->hA2sh), dontSendNotification);

        /* update view windows */
        sensorCoordsView_handle->setQ(array2sh_getNumSensors(hVst->hA2sh));

        needScreenRefreshFLAG = true;
        //[/UserComboBoxCode_presetCB]
    }
    else if (comboBoxThatHasChanged == arrayTypeCB.get())
    {
        //[UserComboBoxCode_arrayTypeCB] -- add your combo box handling code here..
        array2sh_setArrayType(hVst->hA2sh, arrayTypeCB->getSelectedId());
        needScreenRefreshFLAG = true;
        //[/UserComboBoxCode_arrayTypeCB]
    }
    else if (comboBoxThatHasChanged == weightTypeCB.get())
    {
        //[UserComboBoxCode_weightTypeCB] -- add your combo box handling code here..
        array2sh_setWeightType(hVst->hA2sh, weightTypeCB->getSelectedId());
        needScreenRefreshFLAG = true;
        //[/UserComboBoxCode_weightTypeCB]
    }
    else if (comboBoxThatHasChanged == regTypeCB.get())
    {
        //[UserComboBoxCode_regTypeCB] -- add your combo box handling code here..
        array2sh_setRegType(hVst->hA2sh, regTypeCB->getSelectedId());
        needScreenRefreshFLAG = true;
        //[/UserComboBoxCode_regTypeCB]
    }
    else if (comboBoxThatHasChanged == CHOrderingCB.get())
    {
        //[UserComboBoxCode_CHOrderingCB] -- add your combo box handling code here..
        array2sh_setChOrder(hVst->hA2sh, CHOrderingCB->getSelectedId());
        //[/UserComboBoxCode_CHOrderingCB]
    }
    else if (comboBoxThatHasChanged == normalisationCB.get())
    {
        //[UserComboBoxCode_normalisationCB] -- add your combo box handling code here..
        array2sh_setNormType(hVst->hA2sh, normalisationCB->getSelectedId());
        //[/UserComboBoxCode_normalisationCB]
    }
    else if (comboBoxThatHasChanged == dispWindow.get())
    {
        //[UserComboBoxCode_dispWindow] -- add your combo box handling code here..
        dispID = (DISP_WINDOW)dispWindow->getSelectedId();
        needScreenRefreshFLAG = 1;
        //[/UserComboBoxCode_dispWindow]
    }
    else if (comboBoxThatHasChanged == CBencodingOrder.get())
    {
        //[UserComboBoxCode_CBencodingOrder] -- add your combo box handling code here..
        int newOrder = CBencodingOrder->getSelectedId();
        array2sh_setEncodingOrder(hVst->hA2sh, newOrder);
        needScreenRefreshFLAG = true;
        //[/UserComboBoxCode_CBencodingOrder]
    }

    //[UsercomboBoxChanged_Post]
    bool shouldBeEnabled;
    shouldBeEnabled = array2sh_getWeightType(hVst->hA2sh) != 1 ? false : true; /* is it a rigid array? */
    RSlider->setEnabled(shouldBeEnabled);
    shouldBeEnabled = array2sh_getArrayType(hVst->hA2sh) != 1 ? false : true;  /* is it a cylindrical array? */
    weightTypeCB->setItemEnabled(3, shouldBeEnabled);
    weightTypeCB->setItemEnabled(4, shouldBeEnabled);
    if( (array2sh_getArrayType(hVst->hA2sh) == 2) && ((array2sh_getWeightType(hVst->hA2sh) == 3) || (array2sh_getWeightType(hVst->hA2sh) == 4)) )
        weightTypeCB->setSelectedId(1, sendNotification); /* force a suitable option */
    shouldBeEnabled = array2sh_getRegType(hVst->hA2sh) == 1 ? false : true;  /* is DaS regularisation selected? */
    regAmountSlider->setEnabled(shouldBeEnabled);
    //[/UsercomboBoxChanged_Post]
}

void PluginEditor::sliderValueChanged (Slider* sliderThatWasMoved)
{
    //[UsersliderValueChanged_Pre]
    //[/UsersliderValueChanged_Pre]

    if (sliderThatWasMoved == QSlider.get())
    {
        //[UserSliderCode_QSlider] -- add your slider handling code here..
        if(array2sh_getMinNumSensors(hVst->hA2sh)<= (int)QSlider->getValue())
            array2sh_setNumSensors(hVst->hA2sh, (int)QSlider->getValue());
        needScreenRefreshFLAG = true;
        //[/UserSliderCode_QSlider]
    }
    else if (sliderThatWasMoved == rSlider.get())
    {
        //[UserSliderCode_rSlider] -- add your slider handling code here..
        array2sh_setr(hVst->hA2sh, (float)rSlider->getValue());
        bool changeRToo = array2sh_getWeightType(hVst->hA2sh) != 1 ? true : false; /* is it a rigid array? */
        if(changeRToo){
            array2sh_setR(hVst->hA2sh, (float)rSlider->getValue());
            RSlider->setValue(rSlider->getValue(), dontSendNotification);
        }
        needScreenRefreshFLAG = true;
        //[/UserSliderCode_rSlider]
    }
    else if (sliderThatWasMoved == RSlider.get())
    {
        //[UserSliderCode_RSlider] -- add your slider handling code here..
        array2sh_setR(hVst->hA2sh, (float)RSlider->getValue());
        needScreenRefreshFLAG = true;
        //[/UserSliderCode_RSlider]
    }
    else if (sliderThatWasMoved == cSlider.get())
    {
        //[UserSliderCode_cSlider] -- add your slider handling code here..
        array2sh_setc(hVst->hA2sh, (float)cSlider->getValue());
        needScreenRefreshFLAG = true;
        //[/UserSliderCode_cSlider]
    }
    else if (sliderThatWasMoved == regAmountSlider.get())
    {
        //[UserSliderCode_regAmountSlider] -- add your slider handling code here..
        array2sh_setRegPar(hVst->hA2sh, (float)regAmountSlider->getValue());
        needScreenRefreshFLAG = true;
        //[/UserSliderCode_regAmountSlider]
    }
    else if (sliderThatWasMoved == maxFreqSlider.get())
    {
        //[UserSliderCode_maxFreqSlider] -- add your slider handling code here..
        array2sh_setMaxFreq(hVst->hA2sh, (float)maxFreqSlider->getValue());
        //[/UserSliderCode_maxFreqSlider]
    }
    else if (sliderThatWasMoved == gainSlider.get())
    {
        //[UserSliderCode_gainSlider] -- add your slider handling code here..
        array2sh_setGain(hVst->hA2sh, (float)gainSlider->getValue());
        //[/UserSliderCode_gainSlider]
    }

    //[UsersliderValueChanged_Post]
    //[/UsersliderValueChanged_Post]
}

void PluginEditor::buttonClicked (Button* buttonThatWasClicked)
{
    //[UserbuttonClicked_Pre]
    //[/UserbuttonClicked_Pre]

    if (buttonThatWasClicked == degRadTB.get())
    {
        //[UserButtonCode_degRadTB] -- add your button handler code here..
        showDegreesInstead = degRadTB->getToggleState();
        sensorCoordsView_handle->setUseDegreesInstead(showDegreesInstead);
        //[/UserButtonCode_degRadTB]
    }
    else if (buttonThatWasClicked == textButton.get())
    {
        //[UserButtonCode_textButton] -- add your button handler code here..
        array2sh_evaluateFilters(hVst->hA2sh);
        //[/UserButtonCode_textButton]
    }
    else if (buttonThatWasClicked == tb_loadJSON.get())
    {
        //[UserButtonCode_tb_loadJSON] -- add your button handler code here..
        FileChooser myChooser ("Load configuration...",
                               hVst->getLastDir().exists() ? hVst->getLastDir() : File::getSpecialLocation (File::userHomeDirectory),
                               "*.json");
        if (myChooser.browseForFileToOpen()) {
            File configFile (myChooser.getResult());
            hVst->setLastDir(configFile.getParentDirectory());
            hVst->loadConfiguration (configFile);
        }
        //[/UserButtonCode_tb_loadJSON]
    }
    else if (buttonThatWasClicked == tb_saveJSON.get())
    {
        //[UserButtonCode_tb_saveJSON] -- add your button handler code here..
        FileChooser myChooser ("Save configuration...",
                               hVst->getLastDir().exists() ? hVst->getLastDir() : File::getSpecialLocation (File::userHomeDirectory),
                               "*.json");
        if (myChooser.browseForFileToSave (true)) {
            File configFile (myChooser.getResult());
            hVst->setLastDir(configFile.getParentDirectory());
            hVst->saveConfigurationToFile (configFile);
        }
        //[/UserButtonCode_tb_saveJSON]
    }

    //[UserbuttonClicked_Post]
    //[/UserbuttonClicked_Post]
}



//[MiscUserCode] You can add your own definitions of your custom methods or any other code here...
void PluginEditor::timerCallback()
{
    /* these parameters can change internally */
    RSlider->setValue(array2sh_getR(hVst->hA2sh), dontSendNotification);
    int curOrder = CBencodingOrder->getSelectedId();
    QSlider->setRange((curOrder+1)*(curOrder+1), array2sh_getMaxNumSensors(), 1);
    QSlider->setValue(array2sh_getNumSensors(hVst->hA2sh), dontSendNotification);
    sensorCoordsView_handle->setQ(array2sh_getNumSensors(hVst->hA2sh));
    CBencodingOrder->setSelectedId(array2sh_getEncodingOrder(hVst->hA2sh), dontSendNotification);

	/* Some parameters shouldn't be enabled if playback is ongoing */
	if (hVst->getIsPlaying()) 
		textButton->setEnabled(false);
	else {
		textButton->setEnabled(true);
		array2sh_checkReInit(hVst->hA2sh);
	}

    /* draw magnitude/spatial-correlation/level-difference curves */
    if (needScreenRefreshFLAG || array2sh_getEvalReady(hVst->hA2sh)){
        switch(dispID){
            default:
            case SHOW_EQ:
                eqviewIncluded->setNumCurves(array2sh_getEncodingOrder(hVst->hA2sh)+1);
                eqviewIncluded->setVisible(true);
                cohviewIncluded->setVisible(false);
                ldiffviewIncluded->setVisible(false);
                eqviewIncluded->repaint();
                break;
            case SHOW_SPATIAL_COH:
                cohviewIncluded->setNumCurves(array2sh_getEncodingOrder(hVst->hA2sh)+1);
                eqviewIncluded->setVisible(false);
                cohviewIncluded->setVisible(true);
                ldiffviewIncluded->setVisible(false);
                cohviewIncluded->repaint();
                break;
            case SHOW_LEVEL_DIFF:
                ldiffviewIncluded->setNumCurves(array2sh_getEncodingOrder(hVst->hA2sh)+1);
                eqviewIncluded->setVisible(false);
                cohviewIncluded->setVisible(false);
                ldiffviewIncluded->setVisible(true);
                ldiffviewIncluded->repaint();
                break;
        }
        needScreenRefreshFLAG = false;
    }

    /* gain range should change with order */
    //int curOrder = array2sh_getEncodingOrder(hVst->hA2sh);
    gainSlider->setRange (-10-(curOrder*curOrder), 5+curOrder, 0.01);

    /* Hide decoding orders that are unsuitable for number of sensors */
    for(int i=1; i<=7; i++)
        CBencodingOrder->setItemEnabled(i, (i+1)*(i+1) <= array2sh_getNumSensors(hVst->hA2sh) ? true : false);

    /* display warning message, if needed */
    if ((hVst->getCurrentBlockSize() % FRAME_SIZE) != 0){
        currentWarning = k_warning_frameSize;
        repaint(0,0,getWidth(),32);
    }
    else if ((hVst->getCurrentNumInputs() < array2sh_getNumSensors(hVst->hA2sh))){
        currentWarning = k_warning_NinputCH;
        repaint(0,0,getWidth(),32);
    }
    else if ((hVst->getCurrentNumOutputs() < array2sh_getNSHrequired(hVst->hA2sh))){
        currentWarning = k_warning_NoutputCH;
        repaint(0,0,getWidth(),32);
    }
    else if(currentWarning){
        currentWarning = k_warning_none;
        repaint(0,0,getWidth(),32);
    }
}



//[/MiscUserCode]


//==============================================================================
#if 0
/*  -- Projucer information section --

    This is where the Projucer stores the metadata that describe this GUI layout, so
    make changes in here at your peril!

BEGIN_JUCER_METADATA

<JUCER_COMPONENT documentType="Component" className="PluginEditor" componentName=""
                 parentClasses="public AudioProcessorEditor, public Timer" constructorParams="PluginProcessor* ownerFilter"
                 variableInitialisers="AudioProcessorEditor(ownerFilter)" snapPixels="8"
                 snapActive="1" snapShown="1" overlayOpacity="0.330" fixedSize="1"
                 initialWidth="800" initialHeight="450">
  <BACKGROUND backgroundColour="ffffffff">
    <RECT pos="0 30 800 420" fill=" radial: 408 224, 688 392, 0=ff55636d, 1=ff073642"
          hasStroke="1" stroke="1.9, mitered, butt" strokeColour="solid: ffa3a4a5"/>
    <RECT pos="500 264 284 172" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <RECT pos="12 87 204 106" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <RECT pos="0 0 800 32" fill="solid: ff073642" hasStroke="1" stroke="2.7, mitered, butt"
          strokeColour="solid: dcbdbdbd"/>
    <TEXT pos="-19 0 195 32" fill="solid: ffffffff" hasStroke="0" text="Array2SH"
          fontname="Default font" fontsize="18.80000000000000071054" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="36" typefaceStyle="Bold"/>
    <TEXT pos="20 57 67 30" fill="solid: ffffffff" hasStroke="0" text="Presets: "
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="20 89 92 30" fill="solid: ffffffff" hasStroke="0" text="No. Sensors: "
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="20 121 116 30" fill="solid: ffffffff" hasStroke="0" text="Array r (m):"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="20 153 124 30" fill="solid: ffffffff" hasStroke="0" text="Sensor r (m):"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <RECT pos="12 192 204 244" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <RECT pos="228 264 273 172" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <TEXT pos="244 305 172 30" fill="solid: ffffffff" hasStroke="0" text="c (m/s):"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <RECT pos="12 56 204 32" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <TEXT pos="244 337 172 30" fill="solid: ffffffff" hasStroke="0" text="Array Type:"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="244 369 172 30" fill="solid: ffffffff" hasStroke="0" text="Weight Type:"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="524 273 172 30" fill="solid: ffffffff" hasStroke="0" text="Reg. Type:"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="524 305 172 30" fill="solid: ffffffff" hasStroke="0" text="Max Gain (dB):"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="523 369 172 30" fill="solid: ffffffff" hasStroke="0" text="CH Order:"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="523 401 172 30" fill="solid: ffffffff" hasStroke="0" text="Normalisation:"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <RECT pos="228 56 556 209" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <TEXT pos="523 337 172 30" fill="solid: ffffffff" hasStroke="0" text="Post Gain (dB):"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="243 401 172 30" fill="solid: ffffffff" hasStroke="0" text="Max Freq. (Hz):"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="75 193 85 30" fill="solid: ffffffff" hasStroke="0" text="Azi   #   Elev"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <RECT pos="24 224 184 200" fill="solid: 39a52a" hasStroke="1" stroke="1, mitered, butt"
          strokeColour="solid: 29b6b5b5"/>
    <TEXT pos="179 193 24 23" fill="solid: ffffffff" hasStroke="0" text="&#176;"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="92 30 88 30" fill="solid: ffffffff" hasStroke="0" text="Inputs"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <TEXT pos="440 30 149 30" fill="solid: ffffffff" hasStroke="0" text="Encoding Settings"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
    <RECT pos="228 264 273 41" fill="solid: 13f4f4f4" hasStroke="1" stroke="0.8, mitered, butt"
          strokeColour="solid: 67a0a0a0"/>
    <TEXT pos="243 270 172 30" fill="solid: ffffffff" hasStroke="0" text="Encoding Order:"
          fontname="Default font" fontsize="15.00000000000000000000" kerning="0.00000000000000000000"
          bold="1" italic="0" justification="33" typefaceStyle="Bold"/>
  </BACKGROUND>
  <COMBOBOX name="new combo box" id="abcd469891fabf2d" memberName="presetCB"
            virtualName="" explicitFocusOrder="0" pos="88 64 120 16" editable="0"
            layout="33" items="Default" textWhenNonSelected="" textWhenNoItems="(no choices)"/>
  <COMBOBOX name="new combo box" id="ee4eaf9d9b41f219" memberName="arrayTypeCB"
            virtualName="" explicitFocusOrder="0" pos="360 344 120 16" editable="0"
            layout="33" items="Spherical&#10;Cylindrical" textWhenNonSelected="Spherical"
            textWhenNoItems="(no choices)"/>
  <SLIDER name="new slider" id="93dd93c125dcb3b3" memberName="QSlider"
          virtualName="" explicitFocusOrder="0" pos="112 96 96 16" min="4.00000000000000000000"
          max="64.00000000000000000000" int="1.00000000000000000000" style="LinearHorizontal"
          textBoxPos="TextBoxRight" textBoxEditable="1" textBoxWidth="55"
          textBoxHeight="20" skewFactor="1.00000000000000000000" needsCallback="1"/>
  <SLIDER name="new slider" id="f6fe97a46cc051e2" memberName="rSlider"
          virtualName="" explicitFocusOrder="0" pos="112 128 96 16" min="0.01000000000000000021"
          max="0.29999999999999998890" int="0.00100000000000000002" style="LinearHorizontal"
          textBoxPos="TextBoxRight" textBoxEditable="1" textBoxWidth="55"
          textBoxHeight="20" skewFactor="1.00000000000000000000" needsCallback="1"/>
  <SLIDER name="new slider" id="41d575604e9dcd3a" memberName="RSlider"
          virtualName="" explicitFocusOrder="0" pos="112 160 96 16" min="0.01000000000000000021"
          max="0.29999999999999998890" int="0.00100000000000000002" style="LinearHorizontal"
          textBoxPos="TextBoxRight" textBoxEditable="1" textBoxWidth="55"
          textBoxHeight="20" skewFactor="1.00000000000000000000" needsCallback="1"/>
  <SLIDER name="new slider" id="c3b6d9fc71650ee4" memberName="cSlider"
          virtualName="" explicitFocusOrder="0" pos="360 313 120 16" min="200.00000000000000000000"
          max="2000.00000000000000000000" int="0.10000000000000000555"
          style="LinearHorizontal" textBoxPos="TextBoxRight" textBoxEditable="1"
          textBoxWidth="55" textBoxHeight="20" skewFactor="1.00000000000000000000"
          needsCallback="1"/>
  <COMBOBOX name="new combo box" id="ea26910fd5e03b81" memberName="weightTypeCB"
            virtualName="" explicitFocusOrder="0" pos="360 377 120 16" editable="0"
            layout="33" items="Rigid&#10;Open_omni&#10;Open_card&#10;Open_dipole"
            textWhenNonSelected="Rigid" textWhenNoItems="(no choices)"/>
  <COMBOBOX name="new combo box" id="d818d0d5310dc52a" memberName="regTypeCB"
            virtualName="" explicitFocusOrder="0" pos="640 280 128 16" editable="0"
            layout="33" items="DaS&#10;Soft Limiting&#10;Tikhonov" textWhenNonSelected="Soft Limiting"
            textWhenNoItems="(no choices)"/>
  <SLIDER name="new slider" id="9f4f4ac547d19161" memberName="regAmountSlider"
          virtualName="" explicitFocusOrder="0" pos="640 312 128 16" min="0.00000000000000000000"
          max="80.00000000000000000000" int="0.01000000000000000021" style="LinearHorizontal"
          textBoxPos="TextBoxRight" textBoxEditable="1" textBoxWidth="55"
          textBoxHeight="20" skewFactor="1.00000000000000000000" needsCallback="1"/>
  <COMBOBOX name="new combo box" id="44b90530e58253e" memberName="CHOrderingCB"
            virtualName="" explicitFocusOrder="0" pos="640 376 128 16" editable="0"
            layout="33" items="ACN" textWhenNonSelected="ACN" textWhenNoItems="(no choices)"/>
  <COMBOBOX name="new combo box" id="caeee0fc74db72a4" memberName="normalisationCB"
            virtualName="" explicitFocusOrder="0" pos="640 408 128 16" editable="0"
            layout="33" items="N3D&#10;SN3D" textWhenNonSelected="N3D" textWhenNoItems="(no choices)"/>
  <SLIDER name="new slider" id="1835b85fb537168c" memberName="maxFreqSlider"
          virtualName="" explicitFocusOrder="0" pos="360 408 120 16" min="5000.00000000000000000000"
          max="24000.00000000000000000000" int="1.00000000000000000000"
          style="LinearHorizontal" textBoxPos="TextBoxRight" textBoxEditable="1"
          textBoxWidth="55" textBoxHeight="20" skewFactor="1.00000000000000000000"
          needsCallback="1"/>
  <SLIDER name="new slider" id="ee4c42494881e7dc" memberName="gainSlider"
          virtualName="" explicitFocusOrder="0" pos="640 344 128 16" min="-60.00000000000000000000"
          max="60.00000000000000000000" int="0.01000000000000000021" style="LinearHorizontal"
          textBoxPos="TextBoxRight" textBoxEditable="1" textBoxWidth="55"
          textBoxHeight="20" skewFactor="1.00000000000000000000" needsCallback="1"/>
  <TOGGLEBUTTON name="new toggle button" id="3fdbf3711d00f4db" memberName="degRadTB"
                virtualName="" explicitFocusOrder="0" pos="186 198 23 24" buttonText=""
                connectedEdges="0" needsCallback="1" radioGroupId="0" state="0"/>
  <TEXTBUTTON name="new button" id="dde3f82641b3717c" memberName="textButton"
              virtualName="" explicitFocusOrder="0" pos="712 37 72 16" bgColOff="ff5c68a4"
              buttonText="Analyse" connectedEdges="0" needsCallback="1" radioGroupId="0"/>
  <COMBOBOX name="new combo box" id="fb3d8d6828195921" memberName="dispWindow"
            virtualName="" explicitFocusOrder="0" pos="642 37 63 16" editable="0"
            layout="33" items="EQ&#10;Coh&#10;L Diff" textWhenNonSelected="EQ"
            textWhenNoItems="(no choices)"/>
  <TEXTBUTTON name="new button" id="527e24c6748d02d4" memberName="tb_loadJSON"
              virtualName="" explicitFocusOrder="0" pos="148 39 34 14" bgColOff="ff14889e"
              buttonText="Import" connectedEdges="2" needsCallback="1" radioGroupId="0"/>
  <TEXTBUTTON name="new button" id="48c5d3526dcfe64f" memberName="tb_saveJSON"
              virtualName="" explicitFocusOrder="0" pos="182 39 34 14" bgColOff="ff224d97"
              bgColOn="ff181f9a" buttonText="Export" connectedEdges="1" needsCallback="1"
              radioGroupId="0"/>
  <COMBOBOX name="new combo box" id="a465903000494955" memberName="CBencodingOrder"
            virtualName="" explicitFocusOrder="0" pos="363 276 118 20" editable="0"
            layout="33" items="" textWhenNonSelected="Default" textWhenNoItems="(no choices)"/>
</JUCER_COMPONENT>

END_JUCER_METADATA
*/
#endif


//[EndFile] You can add extra defines here...
//[/EndFile]
