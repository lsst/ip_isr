<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>interpolateOverMaskedPixels.h</name>
    <path>/astro/users/nms/code/isrtrunk/include/lsst/ip/isr/</path>
    <filename>interpolate_over_masked_pixels_8h</filename>
    <namespace>lsst</namespace>
    <namespace>lsst::ip</namespace>
    <namespace>lsst::ip::isr</namespace>
    <member kind="typedef">
      <type>boost::uint16_t</type>
      <name>maskPixelType</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>b013f9a2cf7ac68a3e1ded9f5fdb4a25</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::MaskedImage&lt; ImageT, MaskT &gt;</type>
      <name>interpolateOverMaskedPixels</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>6baca3e106c27028ba74a1cae054c44b</anchor>
      <arglist>(lsst::afw::image::MaskedImage&lt; ImageT, MaskT &gt; &amp;chunkMaskedImage, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, std::string const interpMethod)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>isr.h</name>
    <path>/astro/users/nms/code/isrtrunk/include/lsst/ip/isr/</path>
    <filename>isr_8h</filename>
    <namespace>lsst</namespace>
    <namespace>lsst::ip</namespace>
    <namespace>lsst::ip::isr</namespace>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>assembleChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>9b419f80938aae7420089c588c530b8f</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Image&lt; ImageT &gt; &amp;chunkImage, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::afw::image::Mask&lt; MaskT &gt; &amp;badPixelMask)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>saturationCorrectionForChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>f06d4a775dd3589b7336828489cd3330</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>overscanCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>ba77af464fdfc73e33eaa6f844c12ea7</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>trimChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>c66e8aac545e16d0aade66d7f74b303d</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>biasCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>f1e2ab254fc715c0ae5e8696bf459631</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>darkCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>a78a461df074c4275c067ed601e8bf05</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>maskBadPixelsInChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>4ef577270911258a9c790e81a13f2f2d</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Mask&lt; MaskT &gt; const &amp;badPixelMask, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>linearizeChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>b374d32a981f98e94932fef7928fd86a</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::afw::function::Function2&lt; FunctionT &gt; const &amp;function, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>flatCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>b990b8654c9f9f0736e5bec417eb89ec</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>scatteredLightCorrection</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>a9ef8930973190ac37dcc23bc57d9f6b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>geometricDistortionCorrection</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>c0175fb10f9feb7ad895e65ac645ca8f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>pupilImageCorrection</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>d679ae45c17dc937e4a238805bfe5963</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>maskOpticalGhostsAndDiffractionSpikes</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>31d2c32ff7e95c27d91e68a65f324a6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>defringeChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>e07e4f2ad1f7f7234dac257f7a123f1e</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>interpolateOverMaskedPixels</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>aafffa951b366412690dff58bde9bec4</anchor>
      <arglist>(lsst::afw:image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>startPipeline.py</name>
    <path>/astro/users/nms/code/isrtrunk/pipeline/examples/</path>
    <filename>start_pipeline_8py</filename>
    <namespace>startPipeline</namespace>
    <member kind="function">
      <type>def</type>
      <name>startPipeline</name>
      <anchorfile>namespacestart_pipeline.html</anchorfile>
      <anchor>16e83f552f0604ace6c7907e7982252c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>parseNodeList</name>
      <anchorfile>namespacestart_pipeline.html</anchorfile>
      <anchor>a216a7a359beee36f3591dd9e90cf9af</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>__init__.py</name>
    <path>/astro/users/nms/code/isrtrunk/python/lsst/ip/</path>
    <filename>____init_____8py</filename>
    <namespace>ip</namespace>
  </compound>
  <compound kind="file">
    <name>MetadataManipStage.py</name>
    <path>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/pipeline/</path>
    <filename>_metadata_manip_stage_8py</filename>
    <namespace>MetadataManipStage</namespace>
    <class kind="class">MetadataManipStage::MetadataManipStage</class>
  </compound>
  <compound kind="file">
    <name>assembleChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>assemble_chunk_exposure_8cc</filename>
    <includes id="interpolate_over_masked_pixels_8h" name="interpolateOverMaskedPixels.h" local="no" imported="no">lsst/ip/isr/interpolateOverMaskedPixels.h</includes>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>assemble_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>s timing lsst::fw::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>assembleChunkExposure</name>
      <anchorfile>assemble_chunk_exposure_8cc.html</anchorfile>
      <anchor>e7a95ce584d3fc83bd441d7effd40b1b</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Image&lt; ImageT &gt; &amp;chunkImage, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::afw::image::Mask&lt; MaskT &gt; &amp;badPixelMask, std::string const interpMethod)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>biasCorrectChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>bias_correct_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>bias_correct_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>crosstalkCorrectChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>crosstalk_correct_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>crosstalk_correct_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>darkCurrentCorrectChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>dark_current_correct_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>dark_current_correct_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>defringeChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>defringe_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>defringe_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>flatCorrectChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>flat_correct_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>flat_correct_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>geometricDistortionCorrection.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>geometric_distortion_correction_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>geometric_distortion_correction_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>illuminationCorrection.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>illumination_correction_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>illumination_correction_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>instrumentSignatureRemovalController.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>instrument_signature_removal_controller_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>instrument_signature_removal_controller_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>interpolateOverMaskedPixels.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>interpolate_over_masked_pixels_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>interpolate_over_masked_pixels_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>lsst::fw::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>interpolateOverMaskedpixels</name>
      <anchorfile>interpolate_over_masked_pixels_8cc.html</anchorfile>
      <anchor>63800a869cd26c62b1d1c758132ca787</anchor>
      <arglist>(lsst::afw::Exposure&lt; ImageT, MaskT &gt; &amp;chunkMaskedImage, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, std::string interpMethod)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>linearizeChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>linearize_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>linearize_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>maskAdditionalArtifacts.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>mask_additional_artifacts_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>mask_additional_artifacts_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>overscanCorrectChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>overscan_correct_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>overscan_correct_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>pupilImageCorrection.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>pupil_image_correction_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>pupil_image_correction_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>saturationCorrectionForChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>saturation_correction_for_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>saturation_correction_for_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>saturationCorrectionForChunkExposure</name>
      <anchorfile>saturation_correction_for_chunk_exposure_8cc.html</anchorfile>
      <anchor>72ca3a29c66fd4371eb6a05802c3a7c4</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>scatteredLightCorrection.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>scattered_light_correction_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>scattered_light_correction_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>trimChunkExposure.cc</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>trim_chunk_exposure_8cc</filename>
    <includes id="isr_8h" name="isr.h" local="yes" imported="no">lsst/ip/isr/isr.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>LSST_MAX_TRACE</name>
      <anchorfile>trim_chunk_exposure_8cc.html</anchorfile>
      <anchor>38976cfd5857871a4654ca0b100ee08e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>ip</name>
    <filename>namespaceip.html</filename>
  </compound>
  <compound kind="namespace">
    <name>lsst</name>
    <filename>namespacelsst.html</filename>
    <namespace>lsst::ip</namespace>
  </compound>
  <compound kind="namespace">
    <name>lsst::ip</name>
    <filename>namespacelsst_1_1ip.html</filename>
    <namespace>lsst::ip::isr</namespace>
  </compound>
  <compound kind="namespace">
    <name>lsst::ip::isr</name>
    <filename>namespacelsst_1_1ip_1_1isr.html</filename>
    <member kind="typedef">
      <type>boost::uint16_t</type>
      <name>maskPixelType</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>b013f9a2cf7ac68a3e1ded9f5fdb4a25</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::MaskedImage&lt; ImageT, MaskT &gt;</type>
      <name>interpolateOverMaskedPixels</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>6baca3e106c27028ba74a1cae054c44b</anchor>
      <arglist>(lsst::afw::image::MaskedImage&lt; ImageT, MaskT &gt; &amp;chunkMaskedImage, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, std::string const interpMethod)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>assembleChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>9b419f80938aae7420089c588c530b8f</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Image&lt; ImageT &gt; &amp;chunkImage, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::afw::image::Mask&lt; MaskT &gt; &amp;badPixelMask)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>saturationCorrectionForChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>f06d4a775dd3589b7336828489cd3330</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>overscanCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>ba77af464fdfc73e33eaa6f844c12ea7</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>trimChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>c66e8aac545e16d0aade66d7f74b303d</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>biasCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>f1e2ab254fc715c0ae5e8696bf459631</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>darkCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>a78a461df074c4275c067ed601e8bf05</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>maskBadPixelsInChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>4ef577270911258a9c790e81a13f2f2d</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Mask&lt; MaskT &gt; const &amp;badPixelMask, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>linearizeChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>b374d32a981f98e94932fef7928fd86a</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::afw::function::Function2&lt; FunctionT &gt; const &amp;function, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>flatCorrectChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>b990b8654c9f9f0736e5bec417eb89ec</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>scatteredLightCorrection</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>a9ef8930973190ac37dcc23bc57d9f6b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>geometricDistortionCorrection</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>c0175fb10f9feb7ad895e65ac645ca8f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>pupilImageCorrection</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>d679ae45c17dc937e4a238805bfe5963</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>maskOpticalGhostsAndDiffractionSpikes</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>31d2c32ff7e95c27d91e68a65f324a6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>defringeChunkExposure</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>e07e4f2ad1f7f7234dac257f7a123f1e</anchor>
      <arglist>(lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::afw::image::Exposure&lt; ImageT, MaskT &gt; &amp;masterExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::daf::data::DataProperty::PtrType &amp;masterMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
    <member kind="function">
      <type>lsst::afw::image::Exposure&lt; ImageT, MaskT &gt;</type>
      <name>interpolateOverMaskedPixels</name>
      <anchorfile>namespacelsst_1_1ip_1_1isr.html</anchorfile>
      <anchor>aafffa951b366412690dff58bde9bec4</anchor>
      <arglist>(lsst::afw:image::Exposure&lt; ImageT, MaskT &gt; &amp;chunkExposure, lsst::daf::data::DataProperty::PtrType &amp;chunkMetaData, lsst::pex::policy::Policy &amp;policy)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>MetadataManipStage</name>
    <filename>namespace_metadata_manip_stage.html</filename>
    <class kind="class">MetadataManipStage::MetadataManipStage</class>
  </compound>
  <compound kind="class">
    <name>MetadataManipStage::MetadataManipStage</name>
    <filename>class_metadata_manip_stage_1_1_metadata_manip_stage.html</filename>
    <member kind="function">
      <type>def</type>
      <name>preprocess</name>
      <anchorfile>class_metadata_manip_stage_1_1_metadata_manip_stage.html</anchorfile>
      <anchor>e066b6e53ec8f6d899d3aea7eae80cb8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>activeClipboard</name>
      <anchorfile>class_metadata_manip_stage_1_1_metadata_manip_stage.html</anchorfile>
      <anchor>8493617fce6a3fafb8b244f80a9ffc69</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>startPipeline</name>
    <filename>namespacestart_pipeline.html</filename>
    <member kind="function">
      <type>def</type>
      <name>startPipeline</name>
      <anchorfile>namespacestart_pipeline.html</anchorfile>
      <anchor>16e83f552f0604ace6c7907e7982252c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>parseNodeList</name>
      <anchorfile>namespacestart_pipeline.html</anchorfile>
      <anchor>a216a7a359beee36f3591dd9e90cf9af</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/pipeline/examples/</name>
    <path>/astro/users/nms/code/isrtrunk/pipeline/examples/</path>
    <filename>dir_f5ec0b654371a70ca88f9cb8435beabd.html</filename>
    <file>startPipeline.py</file>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/include/</name>
    <path>/astro/users/nms/code/isrtrunk/include/</path>
    <filename>dir_f96beea7452e15a44a1b60d8464fbc9c.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/include/lsst/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/python/lsst/ip/</name>
    <path>/astro/users/nms/code/isrtrunk/python/lsst/ip/</path>
    <filename>dir_87da45fc59e0cdc03d7936e80d0c4700.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/</dir>
    <file>__init__.py</file>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/include/lsst/ip/</name>
    <path>/astro/users/nms/code/isrtrunk/include/lsst/ip/</path>
    <filename>dir_94b62009528b501ba3450e19e8e6919a.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/include/lsst/ip/isr/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/</name>
    <path>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/</path>
    <filename>dir_d31f46f29035f6d18dd6833cc3266e10.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/pipeline/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/include/lsst/ip/isr/</name>
    <path>/astro/users/nms/code/isrtrunk/include/lsst/ip/isr/</path>
    <filename>dir_f806a12b2f032300422a54c6caa41895.html</filename>
    <file>interpolateOverMaskedPixels.h</file>
    <file>isr.h</file>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/python/lsst/</name>
    <path>/astro/users/nms/code/isrtrunk/python/lsst/</path>
    <filename>dir_afb10e17c1a2879c53da1d254147fbd9.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/python/lsst/ip/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/include/lsst/</name>
    <path>/astro/users/nms/code/isrtrunk/include/lsst/</path>
    <filename>dir_cba010ac5ffb510b3a6034ae41c916bb.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/include/lsst/ip/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/pipeline/</name>
    <path>/astro/users/nms/code/isrtrunk/python/lsst/ip/isr/pipeline/</path>
    <filename>dir_049e4bc73b014aca8c8761b152449fb4.html</filename>
    <file>MetadataManipStage.py</file>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/pipeline/</name>
    <path>/astro/users/nms/code/isrtrunk/pipeline/</path>
    <filename>dir_7fde17d3085bc4772f630d743fb50db5.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/pipeline/examples/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/python/</name>
    <path>/astro/users/nms/code/isrtrunk/python/</path>
    <filename>dir_f9c231dd2f00d9dec55626b72106afac.html</filename>
    <dir>/astro/users/nms/code/isrtrunk/python/lsst/</dir>
  </compound>
  <compound kind="dir">
    <name>/astro/users/nms/code/isrtrunk/src/</name>
    <path>/astro/users/nms/code/isrtrunk/src/</path>
    <filename>dir_983150d5428f77374e6bdea86632e9e6.html</filename>
    <file>assembleChunkExposure.cc</file>
    <file>biasCorrectChunkExposure.cc</file>
    <file>crosstalkCorrectChunkExposure.cc</file>
    <file>darkCurrentCorrectChunkExposure.cc</file>
    <file>defringeChunkExposure.cc</file>
    <file>flatCorrectChunkExposure.cc</file>
    <file>geometricDistortionCorrection.cc</file>
    <file>illuminationCorrection.cc</file>
    <file>instrumentSignatureRemovalController.cc</file>
    <file>interpolateOverMaskedPixels.cc</file>
    <file>linearizeChunkExposure.cc</file>
    <file>maskAdditionalArtifacts.cc</file>
    <file>overscanCorrectChunkExposure.cc</file>
    <file>pupilImageCorrection.cc</file>
    <file>saturationCorrectionForChunkExposure.cc</file>
    <file>scatteredLightCorrection.cc</file>
    <file>trimChunkExposure.cc</file>
  </compound>
</tagfile>
