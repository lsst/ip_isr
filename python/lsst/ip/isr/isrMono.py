# This file should never be included in an active stack.
# It is a one-off proof-of-concept example.
import lsst.pipe.base.connectionTypes as cT
import lsst.pipe.base as pipeBase

from . import isrTask


__all__ = ["IsrMonoTask", "IsrMonoTaskConfig",
           "BiasMonoTask", "BiasMonoTaskConfig",
           "DarkMonoTask", "DarkMonoTaskConfig",
           "FlatMonoTask", "FlatMonoTaskConfig"]


class BiasMonoTaskConnections(pipeBase.PipelineTaskConnections,
                             dimensions={"instrument", "visit", "detector"},
                             defaultTemplates={}):

    outputExposure = cT.Output(
        name="monoISRCCD",
        doc="Output ISR processed exposure; calibrations unrestricted",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )

    ccdExposure = cT.PrerequisiteInput(
        name="raw",
        doc="Input exposure to process.",
        storageClass="Exposure",
        dimensions=["instrument", "visit", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera to construct complete exposures.",
        dimensions=["instrument", "calibration_label"],
    )


class BiasMonoTaskConfig(isrTask.IsrTaskConfig,
                         pipelineConnections=BiasMonoTaskConnections):
    pass


class BiasMonoTask(isrTask.IsrTask):
    ConfigClass = BiasMonoTaskConfig
    _DefaultName = 'monoBias'

    pass


#########
class DarkMonoTaskConnections(BiasMonoTaskConnections):
    bias = cT.Input(
        name="bias",
        doc="Input bias calibration; unrestricted.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector"],
    )


class DarkMonoTaskConfig(isrTask.IsrTaskConfig,
                         pipelineConnections=DarkMonoTaskConnections):
    pass


class DarkMonoTask(isrTask.IsrTask):
    ConfigClass = DarkMonoTaskConfig
    _DefaultName = 'monoDark'

    pass


##########
class FlatMonoTaskConnections(DarkMonoTaskConnections):
    dark = cT.Input(
        name="dark",
        doc="Input dark calibration; unrestricted.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector"],
    )


class FlatMonoTaskConfig(isrTask.IsrTaskConfig,
                         pipelineConnections=FlatMonoTaskConnections):
    pass


class FlatMonoTask(isrTask.IsrTask):
    ConfigClass = FlatMonoTaskConfig
    _DefaultName = 'monoFlat'

    pass

###########
class IsrMonoTaskConnections(FlatMonoTaskConnections):
    flat = cT.Input(
        name="flat",
        doc="Input flat calibration; unrestricted.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector", "physical_filter"],
    )


class IsrMonoTaskConfig(isrTask.IsrTaskConfig,
                        pipelineConnections=IsrMonoTaskConnections):
    pass


class IsrMonoTask(isrTask.IsrTask):
    ConfigClass = IsrMonoTaskConfig
    _DefaultName = 'monoIsr'

    pass
