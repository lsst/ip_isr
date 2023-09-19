__all__ = ["IsrTaskLSST", "IsrTaskLSSTConfig"]

import numpy
import math

from . import isrFunctions
from . import isrQa
from . import linearize
from .defects import Defects

from contextlib import contextmanager
from lsst.afw.cameraGeom import NullLinearityType
import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
from lsst.daf.butler import DimensionGraph
from lsst.meas.algorithms.detection import SourceDetectionTask

from .overscan import OverscanCorrectionTask
from .assembleCcdTask import AssembleCcdTask
from .deferredCharge import DeferredChargeTask
from .crosstalk import CrosstalkTask
from .masking import MaskingTask
from .isrStatistics import IsrStatisticsTask
from .isr import maskNans


def crosstalkSourceLookup(datasetType, registry, quantumDataId, collections):
    """Lookup function to identify crosstalkSource entries.

    This should return an empty list under most circumstances.  Only
    when inter-chip crosstalk has been identified should this be
    populated.

    Parameters
    ----------
    datasetType : `str`
        Dataset to lookup.
    registry : `lsst.daf.butler.Registry`
        Butler registry to query.
    quantumDataId : `lsst.daf.butler.ExpandedDataCoordinate`
        Data id to transform to identify crosstalkSources.  The
        ``detector`` entry will be stripped.
    collections : `lsst.daf.butler.CollectionSearch`
        Collections to search through.

    Returns
    -------
    results : `list` [`lsst.daf.butler.DatasetRef`]
        List of datasets that match the query that will be used as
        crosstalkSources.
    """
    newDataId = quantumDataId.subset(DimensionGraph(registry.dimensions, names=["instrument", "exposure"]))
    results = set(registry.queryDatasets(datasetType, collections=collections, dataId=newDataId,
                                         findFirst=True))
    # In some contexts, calling `.expanded()` to expand all data IDs in the
    # query results can be a lot faster because it vectorizes lookups.  But in
    # this case, expandDataId shouldn't need to hit the database at all in the
    # steady state, because only the detector record is unknown and those are
    # cached in the registry.
    return [ref.expanded(registry.expandDataId(ref.dataId, records=newDataId.records)) for ref in results]


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
    dnlLUT = cT.PrerequisiteInput(
        name="dnlLUT",
        doc="Look-up table for differential non-linearity.",
        storageClass="IsrCalib",
        dimensions=["instrument", "exposure", "detector"],
        isCalibration=True,
        # TODO DM 36636
    )
    bias = cT.PrerequisiteInput(
        name="bias",
        doc="Input bias calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    deferredChargeCalib = cT.PrerequisiteInput(
        name="cpCtiCalib",
        doc="Deferred charge/CTI correction dataset.",
        storageClass="IsrCalib",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    linearizer = cT.PrerequisiteInput(
        name='linearizer',
        storageClass="Linearizer",
        doc="Linearity correction calibration.",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can fall back to cameraGeom
    )
    ptc = cT.PrerequisiteInput(
        name="ptc",
        doc="Input Photon Transfer Curve dataset",
        storageClass="PhotonTransferCurveDataset",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    crosstalk = cT.PrerequisiteInput(
        name="crosstalk",
        doc="Input crosstalk object",
        storageClass="CrosstalkCalib",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can fall back to cameraGeom
    )
    crosstalkSources = cT.PrerequisiteInput(
        name="isrOverscanCorrected",
        doc="Overscan corrected input images.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
        deferLoad=True,
        multiple=True,
        lookupFunction=crosstalkSourceLookup,
        minimum=0,  # not needed for all instruments, no config to control this
    )
    defects = cT.PrerequisiteInput(
        name='defects',
        doc="Input defect tables.",
        storageClass="Defects",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    bfKernel = cT.PrerequisiteInput(
        name='brighterFatterKernel',
        doc="Newer complete kernel + gain solutions.",
        storageClass="BrighterFatterKernel",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can use either bfKernel or newBFKernel
    )
    dark = cT.PrerequisiteInput(
        name='dark',
        doc="Input dark calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    outputExposure = cT.Output(
        name='postISRCCD',
        doc="Output ISR processed exposure.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    preInterpExposure = cT.Output(
        name='preInterpISRCCD',
        doc="Output ISR processed exposure, with pixels left uninterpolated.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputStatistics = cT.Output(
        name="isrStatistics",
        doc="Output of additional statistics table.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "detector"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.doDiffNonLinearCorrection is not True:
            self.prerequisiteInputs.remove("dnlLUT")
        if config.doBias is not True:
            self.prerequisiteInputs.remove("bias")
        if config.doDeferredCharge is not True:
            self.prerequisiteInputs.remove("deferredChargeCalib")
        if config.doLinearize is not True:
            self.prerequisiteInputs.remove("linearizer")
        if config.usePtcGains is not True and config.usePtcReadNoise is not True:
            self.prerequisiteInputs.remove("ptc")
        if config.doCrosstalk is not True:
            self.prerequisiteInputs.remove("crosstalkSources")
            self.prerequisiteInputs.remove("crosstalk")
        if config.doDefect is not True:
            self.prerequisiteInputs.remove("defects")
        if config.doBrighterFatter is not True:
            self.prerequisiteInputs.remove("bfKernel")
        if config.doDark is not True:
            self.prerequisiteInputs.remove("dark")

        if config.doWrite is not True:
            self.outputs.remove("outputExposure")
            self.outputs.remove("preInterpExposure")

        if config.doSaveInterpPixels is not True:
            self.outputs.remove("preInterpExposure")

        if config.doCalculateStatistics is not True:
            self.outputs.remove("outputStatistics")


class IsrTaskLSSTConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=IsrTaskLSSTConnections):
    """Configuration parameters for IsrTaskLSST.

    Items are grouped in the order in which they are executed by the task.
    """
    datasetType = pexConfig.Field(
        dtype=str,
        doc="Dataset type for input data; users will typically leave this alone, "
        "but camera-specific ISR tasks will override it.",
        default="raw",
    )
    expectWcs = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Expect input science images to have a WCS (set False for e.g. spectrographs)."
    )
    qa = pexConfig.ConfigField(
        dtype=isrQa.IsrQaConfig,
        doc="QA related configuration options.",
    )
    doHeaderProvenance = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write calibration identifiers into output exposure header.",
    )

    # Differential non-linearity correction.
    doDiffNonLinearCorrection = pexConfig.Field(
        dtype=bool,
        doc="Do differential non-linearity correction?",
        default=True,
    )

    doOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do overscan subtraction?",
        default=True,
    )
    overscan = pexConfig.ConfigurableField(
        target=OverscanCorrectionTask,
        doc="Overscan subtraction task for image segments.",
    )

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











