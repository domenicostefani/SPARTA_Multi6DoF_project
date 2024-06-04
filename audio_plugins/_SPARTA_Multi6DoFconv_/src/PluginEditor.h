/*
  ==============================================================================

  This is an automatically generated GUI class created by the Projucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Projucer version: 6.1.6

  ------------------------------------------------------------------------------

  The Projucer is part of the JUCE library.
  Copyright (c) 2020 - Raw Material Software Limited.

  ==============================================================================
*/

#pragma once

//[Headers]     -- You can add your own extra header files here --

#include "JuceHeader.h"
#include "PluginProcessor.h"
#include "../../resources/SPARTALookAndFeel.h"
#include "sceneView.h"

typedef enum _SPARTA_WARNINGS{
    k_warning_none,
    k_warning_sampleRate_missmatch,
    k_warning_irs_resampled,
    k_warning_nInputs_more_than_64,
    k_warning_nOutputs_more_than_64

}SPARTA_WARNINGS;

#define BEYOND_SAFE_CROSSFADE_FACTOR 5.0f
#define SHOW_DEBUG_DATETIME

//[/Headers]



//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Introjucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class PluginEditor  : public AudioProcessorEditor,
                      public Timer,
                      private FilenameComponentListener,
                      public juce::Slider::Listener,
                      public juce::ComboBox::Listener,
                      public juce::Button::Listener
{
public:
    //==============================================================================
    PluginEditor (PluginProcessor* ownerFilter);
    ~PluginEditor() override;

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.

    /* Refresh coordinate limits based on loaded sofa files*/
    void refreshCoords();

    bool getRefreshSceneViewWindow();

    void setRefreshSceneViewWindow(bool val);

    /* update first partition size combo box and max partition size combo box */
    void setPartComboboxes();

    //[/UserMethods]

    void paint (juce::Graphics& g) override;
    void resized() override;
    void sliderValueChanged (juce::Slider* sliderThatWasMoved) override;
    void comboBoxChanged (juce::ComboBox* comboBoxThatHasChanged) override;
    void buttonClicked (juce::Button* buttonThatWasClicked) override;



private:
    //[UserVariables]   -- You can add your own custom variables in this section.
    PluginProcessor* hVst;
    void* hTVC;
    void* hRot;
    void timerCallback() override;
    bool partitionComboboxesSet = false;

    float maximumSafeCrossfadeMS = 0;
    void updateCrossfadeRange();

    /* Look and Feel */
    SPARTALookAndFeel LAF;

    #if (MCFX_CONVOLVER_MODE == CROSSFADED_CONVOLVERS_MODE) && defined(MCFX_CONVOLVER_MODE) && defined(CROSSFADED_CONVOLVERS_MODE)
        juce::Colour crossfadeSldBackground = juce::Colours::black.withAlpha(0.5f);
        juce::Rectangle<int> crossfadeWarningArea;
    #endif

    /** Custom Look and feel for small text comboboxes */
    class CustomLookAndFeel : public LookAndFeel_V4 {
    private:
        float comboBoxTextWidthPercentage = 0.7f,  // Percentage of combobox width to use for text
            comboBoxMinArrowWidth = 18;            // Minimum width of the arrow in pixels
    public:
        CustomLookAndFeel() = default;
        void drawComboBox(Graphics& g, int width, int height, bool, int, int, int, int, ComboBox& box) {
            auto cornerSize = box.findParentComponentOfClass<ChoicePropertyComponent>() != nullptr ? 0.0f : 3.0f;
            juce::Rectangle<int> boxBounds(0, 0, width, height);

            g.setColour(box.findColour(ComboBox::backgroundColourId));
            g.fillRoundedRectangle(boxBounds.toFloat(), cornerSize);

            g.setColour(box.findColour(ComboBox::outlineColourId));
            g.drawRoundedRectangle(boxBounds.toFloat().reduced(0.5f, 0.5f), cornerSize, 1.0f);

            int arrowWidth = jmax((int)width * (1.0f - comboBoxTextWidthPercentage), comboBoxMinArrowWidth);
            juce::Rectangle<int> arrowZone(width - arrowWidth, 0, arrowWidth, height);
            juce::Path path;
            path.startNewSubPath((float)arrowZone.getX() + 3.0f, (float)arrowZone.getCentreY() - 2.0f);
            path.lineTo((float)arrowZone.getCentreX(), (float)arrowZone.getCentreY() + 3.0f);
            path.lineTo((float)arrowZone.getRight() - 3.0f, (float)arrowZone.getCentreY() - 2.0f);

            g.setColour(box.findColour(ComboBox::arrowColourId).withAlpha((box.isEnabled() ? 0.9f : 0.2f)));
            g.strokePath(path, PathStrokeType(1.0f));
        }

        Font getComboBoxFont(ComboBox& box) override {
            return {jmin(12.0f, (float)box.getHeight() * 0.85f)};
        }

        void positionComboBoxText(ComboBox& box, Label& label) override {
            int width = box.getWidth();
            int arrowWidth = jmax((int)width * (1.0f - comboBoxTextWidthPercentage), comboBoxMinArrowWidth);
            label.setBounds(1, 1,
                            width - arrowWidth,
                            box.getHeight() - 2);

            label.setFont(getComboBoxFont(box));
        }
    } smallComboBoxLookAndFeel;

    /* sofa loading */
    std::unique_ptr<juce::FilenameComponent> fileComp;
    SAF_TVCONV_ERROR_CODES tvConvError;

    /* sofa file loading */
     void filenameComponentChanged (FilenameComponent*) override  {
         partitionComboboxesSet = false;
         String directory = fileComp->getCurrentFile().getFullPathName();
         const char* new_cstring = (const char*)directory.toUTF8();
         tvconv_setSofaFilePath(hTVC, new_cstring);
         refreshCoords();

     }

    /* scene view window */
    std::unique_ptr<sceneView> sceneWindow;
    int refreshInterval             = 40; /*ms (40ms = 25 frames per second) if refreshDecimationFactor = 1 */
    bool refreshSceneViewWindow;
    int refreshDecimationCounter    = 1;
    int targetDecimatedRefreshRate  = 1;

    /* warnings */
    SPARTA_WARNINGS currentWarning;
    SharedResourcePointer<TooltipWindow> tipWindow;
    std::unique_ptr<juce::ComboBox> pluginDescription; /* Dummy combo box to provide plugin description tooltip */

    //[/UserVariables]

    //==============================================================================
    std::unique_ptr<juce::Label> label_hostBlockSize;
    std::unique_ptr<juce::Label> label_filterLength;
    std::unique_ptr<juce::Label> label_hostfs;
    std::unique_ptr<juce::Label> label_filterfs;
    std::unique_ptr<juce::Label> label_nIRpositions;
    std::unique_ptr<juce::Label> label_NInputs;
    std::unique_ptr<juce::Slider> SL_source_y;
    std::unique_ptr<juce::Slider> SL_source_z;
    std::unique_ptr<juce::Slider> SL_source_x;
    std::unique_ptr<juce::Slider> SL_receiver_x;
    std::unique_ptr<juce::Slider> SL_receiver_y;
    std::unique_ptr<juce::Slider> SL_receiver_z;
    std::unique_ptr<juce::TextEditor> te_oscport;
    std::unique_ptr<juce::ComboBox> CBviewMode;
    std::unique_ptr<juce::Slider> s_yaw;
    std::unique_ptr<juce::Slider> s_pitch;
    std::unique_ptr<juce::Slider> s_roll;
    std::unique_ptr<juce::ToggleButton> t_flipYaw;
    std::unique_ptr<juce::ToggleButton> t_flipPitch;
    std::unique_ptr<juce::ToggleButton> t_flipRoll;
    std::unique_ptr<juce::ToggleButton> TBenableRotation;
    std::unique_ptr<juce::Label> label_NOutputs;
    std::unique_ptr<juce::Label> label_NIRs;
    std::unique_ptr<juce::ComboBox> box_first_part;
    std::unique_ptr<juce::ComboBox> box_maxpart;
    std::unique_ptr<juce::Label> label_receiverIdx;
    std::unique_ptr<juce::Slider> SL_crossfadeTimeMs;
    std::unique_ptr<juce::TextButton> btn_doubleCrossfade;
    std::unique_ptr<juce::TextButton> btn_halveCrossfade;


    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PluginEditor)
};

//[EndFile] You can add extra defines here...
//[/EndFile]

