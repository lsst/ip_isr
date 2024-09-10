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
    doParallelOverscanCrosstalk = pexConfig.Field(
        dtype=bool,
        doc="Apply crosstalk correction in parallel overscan region?",
        default=True,
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
                         f"doParallelCrosstalk={self.doParallelOverscanCrosstalk} "
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
    def doAnyParallelOverscanCrosstalk(self):
        """Check if any of the amp configs have doParallelOverscanCrosstalk.

        Returns
        -------
        doAnyParallelOverscanCrosstalk : `bool`
        """
        if self.defaultAmpConfig.doParallelOverscanCrosstalk:
            return True

        for _, ampRule in self.ampRules.items():
            if ampRule.doParallelOverscanCrosstalk:
                return True

        return False

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

    @property
    def doAnyParallelOverscanCrosstalk(self):
        """Check if any of the detector/amp configs have
        doParallelOverscanCrosstalk.

        Returns
        -------
        doAnyParallelOverscanCrosstalk : `bool`
        """

        if self.defaultDetectorConfig.doAnyParallelOverscanCrosstalk:
            return True

        for _, detectorRule in self.detectorRules.items():
            if detectorRule.doAnyParallelOverscanCrosstalk:
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
