import lsst.pex.config as pexConfig

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

    def setDefaults(self):
        super().setDefaults()

        self.serialOverscanConfig.fitType = "MEDIAN_PER_ROW"
        self.serialOverscanConfig.leadingToSkip = 3
        self.serialOverscanConfig.trailingToSkip = 3
        self.parallelOverscanConfig.fitType = "MEDIAN_PER_ROW"
        self.parallelOverscanConfig.leadingToSkip = 3
        self.parallelOverscanConfig.trailingToSkip = 3


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
