class IsrTaskLSSTConnections(pipeBase.PipelineTaskConnections,
                         dimensions={"instrument", "exposure", "detector"},
                         defaultTemplates={}):
    ccdExposure = cT.Input(
        name="raw",
        doc="Input exposure to process.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera to construct complete exposures.",
        dimensions=["instrument"],
        isCalibration=True,
    )
    #Non-linearity correction DM-36636

class isrTaskLSST(pipeBase.PipelineTask):
    ConfigClass = IsrTaskLSSTConfig
    _DefaultName = "isr"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("assembleCcd")
        self.makeSubtask("crosstalk")
        self.makeSubtask("strayLight")
        self.makeSubtask("fringe")
        self.makeSubtask("masking")
        self.makeSubtask("overscan")
        self.makeSubtask("vignette")
        self.makeSubtask("ampOffset")
        self.makeSubtask("deferredChargeCorrection")
        self.makeSubtask("isrStats")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        super().runQuantum(self, butlerQC, inputRefs, outputRefs)

    def validateInput(self,**kwargs):
        """
        This is a check that all the inputs required by the config
        are available.
        """

        inputMap = {'bias': self.config.doBias}
        for calibrationFile,configValue in inputMap.items():
            #TODO do checks for fringes, illumination correction, etc
            if inputMap[calibrationFile] is None and configValue:
                raise RuntimeError("Must supply ",calibrationFile)


    def diffNonLinearCorrection(self,ccdExposure,dnlLUT,**kwargs):
        #TODO DM 36636
        #isrFunctions.diffNonLinearCorrection
        pass

    def overscanCorrection(self, ccd, ccdExposure):
        #TODO DM 36637 for per amp

        overscans = []
        for amp in ccd:

            # Overscan correction on amp-by-amp basis.
            if amp.getRawHorizontalOverscanBBox().isEmpty():
                self.log.info("ISR_OSCAN: No overscan region.  Not performing overscan correction.")
                overscans.append(None)
            else:

                # Perform overscan correction on subregions.
                overscanResults = self.overscan.run(ccdExposure, amp)

                self.log.debug("Corrected overscan for amplifier %s.", amp.getName())
                if len(overscans) == 0:
                    ccdExposure.getMetadata().set('OVERSCAN', "Overscan corrected")

                overscans.append(overscanResults if overscanResults is not None else None)

        return overscans

    def snapCombine(self,**kwargs):
        #TODO DM 36638
        pass

    def gainNormalize(self,**kwargs):
        #TODO DM 36639
        gains = []
        readNoise = []

        return gains, readNoise

    def variancePlane(self, ccdExposure, ccd, overscans, gains, readNoises, **kwargs):
        for amp, gain, readNoise in zip(ccd, gains,readNoises):
            if ccdExposure.getBBox().contains(amp.getBBox()):
                self.log.debug("Constructing variance map for amplifer %s.", amp.getName())
                ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())

                isrFunctions.updateVariance(
                        maskedImage=ampExposure.getMaskedImage(),
                        gain=gain,
                        readNoise=readNoise,
                        )

#                if self.config.maskNegativeVariance:
#                    self.maskNegativeVariance(ccdExposure)

    def run(self,**kwargs):

        ccd = ccdExposure.getDetector()
        filterLabel = ccdExposure.getFilter()
        physicalFilter = isrFunctions.getPhysicalFilter(filterLabel, self.log)

        self.validateInput(**kwargs)
        if self.config.doDiffNonLinearCorrection:
            self.diffNonLinearCorrection(ccdExposure,dnlLUT,**kwargs)
        if self.config.doOverscan:
            overscans = self.overscanCorrection(self, ccd, ccdExposure)

        if self.config.doAssembleCcd:
            self.log.info("Assembling CCD from amplifiers.")
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warning("No WCS found in input exposure.")
            self.debugView(ccdExposure, "doAssembleCcd")

        if self.config.doSnapCombine:
            self.snapCombine(**kwargs)

        if self.config.doBias:
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.getMaskedImage(), bias.getMaskedImage(),
                                        trimToFit=self.config.doTrimToMatchCalib)
            self.debugView(ccdExposure, "doBias")

        if self.config.doDeferredCharge:
            self.log.info("Applying deferred charge/CTI correction.")
            self.deferredChargeCorrection.run(ccdExposure, deferredChargeCalib)
            self.debugView(ccdExposure, "doDeferredCharge")


        if self.config.doLinearize:
            self.log.info("Applying linearizer.")
            linearizer.applyLinearity(image=ccdExposure.getMaskedImage().getImage(),
                                      detector=ccd, log=self.log)

        if self.config.doGainNormalize:
            gains, readNoise = self.gainNormalize(**kwargs)

        if self.config.doVariance:
            self.variancePlane(gains, readNoise, **kwargs)

        if self.config.doCrosstalk:
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalk=crosstalk,
                               crosstalkSources=crosstalkSources, isTrimmed=True)
            self.debugView(ccdExposure, "doCrosstalk")











