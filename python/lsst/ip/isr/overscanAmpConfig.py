import lsst.pex.config as pexConfig
import hashlib

from .overscan import SerialOverscanCorrectionTaskConfig, ParallelOverscanCorrectionTaskConfig


__all__ = [
    "OverscanAmpConfig",
    "OverscanDetectorConfig",
    "OverscanCameraConfig",
]


class OverscanAmpConfig(pexConfig.Config):
    """Overscan configurations applicable to a single amplifier."""
    doSerialOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do serial overscan subtraction?",
        default=True,
    )
    serialOverscanConfig = pexConfig.ConfigField(
        dtype=SerialOverscanCorrectionTaskConfig,
        doc="Serial overscan configuration.",
    )
    # TODO: Remove on DM-48394
    doParallelOverscanCrosstalk = pexConfig.Field(
        dtype=bool,
        doc="Apply crosstalk correction in parallel overscan region?",
        default=True,
        deprecated="This field is no longer used, and will be removed after v29.",
    )
    doParallelOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do parallel overscan subtraction?",
        default=True,
    )
    parallelOverscanConfig = pexConfig.ConfigField(
        dtype=ParallelOverscanCorrectionTaskConfig,
        doc="Parallel overscan configuration.",
    )
    saturation = pexConfig.Field(
        dtype=float,
        doc="The saturation level to use to override any detector/calibration product value "
            "(ignored if NaN). Units are ADU.",
        default=float("NaN"),
    )
    suspectLevel = pexConfig.Field(
        dtype=float,
        doc="The ``suspect`` level to use to override any detector/calibration product value "
            "(ignored if NaN). Units are ADU.",
        default=float("NaN"),
    )
    gain = pexConfig.Field(
        dtype=float,
        doc="The gain to use to override any calibration product value (ignored if NaN). "
            "Units are e-/ADU.",
        default=float("NaN"),
    )
    maskAmpAsBad = pexConfig.Field(
        dtype=bool,
        doc="Mask this amp as BAD (if doDefect=True)?",
        default=False,
    )

    def setDefaults(self):
        super().setDefaults()

        self.serialOverscanConfig.fitType = "MEDIAN_PER_ROW"
        self.serialOverscanConfig.leadingToSkip = 3
        self.serialOverscanConfig.trailingToSkip = 3
        self.serialOverscanConfig.overscanIsInt = False
        self.parallelOverscanConfig.fitType = "MEDIAN_PER_ROW"
        self.parallelOverscanConfig.leadingToSkip = 3
        self.parallelOverscanConfig.trailingToSkip = 3
        self.parallelOverscanConfig.overscanIsInt = False
        # We expect the parallel overscan to not deviate much
        # after serial overscan subtraction and crosstalk correction.
        self.parallelOverscanConfig.maxDeviation = 100.0

    @property
    def _stringForHash(self):
        """Turn this config into a simple string for hashing.

        Only essential data for tracking is returned.

        Returns
        -------
        stringForHash : `str`
        """
        stringForHash = (f"doSerial={self.doSerialOverscan} "
                         f"serialFitType={self.serialOverscanConfig.fitType} "
                         f"doParallel={self.doParallelOverscan} "
                         f"parallelFitType={self.parallelOverscanConfig.fitType}")
        return stringForHash


class OverscanDetectorConfig(pexConfig.Config):
    """Overscan configurations applicable to multiple amplifiers in
    a single detector.
    """
    ampRules = pexConfig.ConfigDictField(
        doc="Amplifier level rules for overscan, keyed by amp name.",
        keytype=str,
        itemtype=OverscanAmpConfig,
        default={},
    )
    defaultAmpConfig = pexConfig.ConfigField(
        dtype=OverscanAmpConfig,
        doc="Default configuration for amplifiers.",
    )
    integerDitherMode = pexConfig.ChoiceField(
        dtype=str,
        doc="Dithering mode to cancel integerization of counts.",
        default="SYMMETRIC",
        allowed={
            "POSITIVE": "Dithering is done with a uniform random in the range [0, 1).",
            "NEGATIVE": "Dithering is done with a uniform random in the range [-1, 0).",
            "SYMMETRIC": "Dithering is done with a uniform random in the range [-0.5, 0.5).",
            "NONE": "No dithering is performed.",
        },
    )
    itlDipMinHeight = pexConfig.Field(
        dtype=int,
        doc="Minimum height for a saturated footprint column to contribute to a dip.",
        default=50,
    )
    itlDipMinWidth = pexConfig.Field(
        dtype=int,
        doc="Minimum number of columns in a saturated footprint with idlDipMinHeight "
            "to contribute to a dip.",
        default=15,
    )
    itlDipMaxWidth = pexConfig.Field(
        dtype=int,
        doc="Maximum number of columns to use for a dip mask.",
        default=50,
    )
    itlDipWidthScale = pexConfig.Field(
        dtype=float,
        doc="Scaling factor to widen saturated core for dip masking.",
        default=1.5,
    )
    itlDipBackgroundFraction = pexConfig.Field(
        dtype=float,
        doc="Fraction of background (scaled by width) that is in the center of the dip. "
            "Only dips that are greater than itlDipMinSkyNoiseFraction will be masked. "
            "If equal to 0.0, dip masking will be skipped.",
        default=0.0,
    )
    itlDipMinBackgroundNoiseFraction = pexConfig.Field(
        dtype=float,
        doc="Only max model dip depth greater than this fraction of the approximate "
            "background noise will be masked.",
        default=0.5,
    )
    itlDipMaxColsPerImage = pexConfig.Field(
        dtype=int,
        doc="Maximum number of columns detected as ``dip`` columns before dip masking "
            "is disabled on the image.",
        default=500,
    )

    @property
    def doAnySerialOverscan(self):
        """Check if any of the amp configs have doSerialOverscan.

        Returns
        -------
        doAnySerialOverscan : `bool`
        """
        if self.defaultAmpConfig.doSerialOverscan:
            return True

        for _, ampRule in self.ampRules.items():
            if ampRule.doSerialOverscan:
                return True

        return False

    @property
    def doAnyParallelOverscan(self):
        """Check if any of the amp configs have doParallelOverscan.

        Returns
        -------
        doAnyParallelOverscan : `bool`
        """
        if self.defaultAmpConfig.doParallelOverscan:
            return True

        for _, ampRule in self.ampRules.items():
            if ampRule.doParallelOverscan:
                return True

        return False

    @property
    def badAmpsToMask(self):
        """Get all amp names that should be masked.

        Returns
        -------
        badAmpsToMask : `list` [`str`]
            List of bad amp names.
        """
        badAmpsToMask = []

        for ampKey, ampRule in self.ampRules.items():
            if ampRule.maskAmpAsBad:
                badAmpsToMask.append(ampKey)

        return badAmpsToMask

    def getOverscanAmpConfig(self, amplifier):
        """Get the OverscanAmpConfig for a specific amplifier.

        Parameters
        ----------
        amplifier : `lsst.afw.cameraGeom.Amplifier`

        Returns
        -------
        overscanAmpConfig : `lsst.ip.isr.overscanAmpConfig.OverscanAmpConfig`
        """
        ampKey = amplifier.getName()

        if ampKey in self.ampRules.keys():
            overscanAmpConfig = self.ampRules[ampKey]
        else:
            overscanAmpConfig = self.defaultAmpConfig

        return overscanAmpConfig

    @property
    def _stringForHash(self):
        """Turn this config into a simple string for hashing.

        Only the default and amps that are different than the
        default are used in the string representation.

        Returns
        -------
        stringForHash : `str`
        """
        defaultString = self.defaultAmpConfig._stringForHash

        stringForHash = f"default: {defaultString}"
        for ampName in self.ampRules:
            ampString = self.ampRules[ampName]._stringForHash
            if ampString != defaultString:
                stringForHash += f" {ampName}: {ampString}"

        return stringForHash

    @property
    def md5(self):
        """Compute the MD5 hash of this config (detector + amps).

        This can be used to ensure overscan configs are consistent.

        Returns
        -------
        md5Hash : `str`
        """
        return hashlib.md5(self._stringForHash.encode("UTF-8")).hexdigest()


class OverscanCameraConfig(pexConfig.Config):
    """Overscan configurations applicable to multiple detectors in
    a single camera.
    """
    detectorRules = pexConfig.ConfigDictField(
        doc="Detector level rules for overscan",
        keytype=str,
        itemtype=OverscanDetectorConfig,
        default={},
    )
    defaultDetectorConfig = pexConfig.ConfigField(
        dtype=OverscanDetectorConfig,
        doc="Default configuration for detectors.",
    )
    detectorRuleKeyType = pexConfig.ChoiceField(
        doc="Detector rule key type.",
        dtype=str,
        default="NAME",
        allowed={
            "NAME": "DetectorRules has a key that is the detector name.",
            "SERIAL": "DetectorRules has a key that is the detector serial number.",
            "ID": "DetectorRules has a key that is the detector id number.",
        },
    )

    @property
    def doAnySerialOverscan(self):
        """Check if any of the detector/amp configs have doSerialOverscan.

        Returns
        -------
        doAnySerialOverscan : `bool`
        """
        if self.defaultDetectorConfig.doAnySerialOverscan:
            return True

        for _, detectorRule in self.detectorRules.items():
            if detectorRule.doAnySerialOverscan:
                return True

        return False

    @property
    def doAnyParallelOverscan(self):
        """Check if any of the detector/amp configs have
        doParallelOverscan.

        Returns
        -------
        doAnyParallelOverscan : `bool`
        """

        if self.defaultDetectorConfig.doAnyParallelOverscan:
            return True

        for _, detectorRule in self.detectorRules.items():
            if detectorRule.doAnyParallelOverscan:
                return True

        return False

    def getOverscanDetectorConfig(self, detector):
        """Get the OverscanDetectorConfig for a specific detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`

        Returns
        -------
        overscanDetectorConfig : `OverscanDetectorConfig`
        """
        match self.detectorRuleKeyType:
            case "NAME":
                key = detector.getName()
            case "SERIAL":
                key = detector.getSerial()
            case "ID":
                key = str(detector.getId())

        if key in self.detectorRules.keys():
            overscanDetectorConfig = self.detectorRules[key]
        else:
            overscanDetectorConfig = self.defaultDetectorConfig

        return overscanDetectorConfig
