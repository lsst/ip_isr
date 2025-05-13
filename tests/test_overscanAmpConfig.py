#
# LSST Data Management System
# Copyright 2023 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import copy
import math
import tempfile
import unittest

import lsst.utils.tests
import lsst.afw.cameraGeom as cameraGeom

from lsst.ip.isr.overscanAmpConfig import (
    OverscanAmpConfig,
    OverscanDetectorConfig,
    OverscanCameraConfig,
)


class OverscanAmpConfigTestCase(lsst.utils.tests.TestCase):
    def _makeCamera(self):
        self.detectors = {
            "detector0": (0, "0010"),
            "detector1": (1, "0011"),
            "detector2": (2, "0012"),
        }
        self.amps = ["amp0", "amp1", "amp2", "amp3"]

        cameraBuilder = cameraGeom.Camera.Builder("Fake Camera")
        for detector in self.detectors:
            detectorBuilder = cameraBuilder.add(detector, self.detectors[detector][0])
            detectorBuilder.setSerial(self.detectors[detector][1])

            for amp in self.amps:
                ampBuilder = cameraGeom.Amplifier.Builder()
                ampBuilder.setName(amp)
                detectorBuilder.append(ampBuilder)

        return cameraBuilder.finish()

    def _serializeAndReadConfig(self, configIn):
        # This bit of code serializes and reads the config
        # to ensure that everything here works as expected
        # with the nested dictionaries of configs.

        with tempfile.NamedTemporaryFile(suffix=".py") as f:
            configIn.save(f.name)
            configOut = OverscanCameraConfig()
            configOut.load(f.name)

        return configOut

    def _checkOverscanConfig(
            self,
            overscanAmpConfig,
            doSerialOverscan=True,
            doParallelOverscan=True,
            serialFitType="MEDIAN_PER_ROW",
            parallelFitType="MEDIAN_PER_ROW",
            saturation=float("NaN"),
            gain=float("NaN"),
            suspectLevel=float("NaN"),
    ):
        self.assertEqual(overscanAmpConfig.doSerialOverscan, doSerialOverscan)
        self.assertEqual(overscanAmpConfig.doParallelOverscan, doParallelOverscan)
        self.assertEqual(overscanAmpConfig.serialOverscanConfig.fitType, serialFitType)
        self.assertEqual(overscanAmpConfig.parallelOverscanConfig.fitType, parallelFitType)
        if math.isnan(saturation):
            self.assertTrue(math.isnan(overscanAmpConfig.saturation))
        else:
            self.assertEqual(overscanAmpConfig.saturation, saturation)
        if math.isnan(gain):
            self.assertTrue(math.isnan(overscanAmpConfig.gain))
        else:
            self.assertEqual(overscanAmpConfig.gain, gain)
        if math.isnan(suspectLevel):
            self.assertTrue(math.isnan(overscanAmpConfig.suspectLevel))
        else:
            self.assertEqual(overscanAmpConfig.suspectLevel, suspectLevel)

    def _checkAnyOverscanConfig(
            self,
            config,
            doSerialOverscan=True,
            doParallelOverscan=True,
    ):
        self.assertEqual(config.doAnySerialOverscan, doSerialOverscan)
        self.assertEqual(config.doAnyParallelOverscan, doParallelOverscan)

    def _checkDetectorOverscanConfig(
            self,
            overscanDetectorConfig,
            integerDitherMode="SYMMETRIC",
    ):
        self.assertEqual(overscanDetectorConfig.integerDitherMode, integerDitherMode)

    def testAmpConfigNoOverrides(self):
        camera = self._makeCamera()

        config = OverscanCameraConfig()

        config = self._serializeAndReadConfig(config)

        for detector in camera:
            for amp in detector:
                detectorConfig = config.getOverscanDetectorConfig(detector)
                ampConfig = detectorConfig.getOverscanAmpConfig(amp)
                self._checkOverscanConfig(ampConfig)

        self._checkAnyOverscanConfig(config)

    def testAmpConfigOverrideDetectorDefault(self):
        camera = self._makeCamera()

        overscanAmpConfig = OverscanAmpConfig(
            doSerialOverscan=False,
            doParallelOverscan=False,
            saturation=100_000.0,
            gain=1.7,
            suspectLevel=90_000.0,
        )

        overscanDetectorConfig = OverscanDetectorConfig(defaultAmpConfig=overscanAmpConfig)

        config = OverscanCameraConfig(defaultDetectorConfig=overscanDetectorConfig)

        config = self._serializeAndReadConfig(config)

        for detector in camera:
            for amp in detector:
                detectorConfig = config.getOverscanDetectorConfig(detector)
                overscanAmpConfig = detectorConfig.getOverscanAmpConfig(amp)
                self._checkOverscanConfig(
                    overscanAmpConfig,
                    doSerialOverscan=False,
                    doParallelOverscan=False,
                    saturation=100_000.0,
                    gain=1.7,
                    suspectLevel=90_000.0,
                )

        self._checkAnyOverscanConfig(
            config,
            doSerialOverscan=False,
            doParallelOverscan=False,
        )

    def testAmpConfigOverrideDetectorDefaultWithOneAmp(self):
        camera = self._makeCamera()

        overscanAmpConfigOverride = OverscanAmpConfig()
        overscanAmpConfigOverride.parallelOverscanConfig.fitType = "MEDIAN"

        overscanDetectorConfig = OverscanDetectorConfig()
        overscanDetectorConfig.ampRules["amp2"] = overscanAmpConfigOverride

        config = OverscanCameraConfig(defaultDetectorConfig=overscanDetectorConfig)

        config = self._serializeAndReadConfig(config)

        for detector in camera:
            for amp in detector:
                detectorConfig = config.getOverscanDetectorConfig(detector)
                ampConfig = detectorConfig.getOverscanAmpConfig(amp)
                if amp.getName() == "amp2":
                    self._checkOverscanConfig(ampConfig, parallelFitType="MEDIAN")
                else:
                    self._checkOverscanConfig(ampConfig)

        self._checkAnyOverscanConfig(config)

    def testAmpConfigOverrideOneDetector(self):
        camera = self._makeCamera()

        overscanAmpConfigOverride = OverscanAmpConfig(doParallelOverscan=False)
        overscanDetectorConfigOverride = OverscanDetectorConfig(
            defaultAmpConfig=overscanAmpConfigOverride,
            integerDitherMode="NONE",
        )

        for keyType in ["NAME", "SERIAL", "ID"]:
            config = OverscanCameraConfig()

            match keyType:
                case "NAME":
                    key = camera[1].getName()
                case "SERIAL":
                    key = camera[1].getSerial()
                case "ID":
                    key = str(camera[1].getId())

            config.detectorRuleKeyType = keyType
            config.detectorRules[key] = overscanDetectorConfigOverride

            config = self._serializeAndReadConfig(config)

            for detector in camera:
                detectorConfig = config.getOverscanDetectorConfig(detector)
                if detector.getName() == camera[1].getName():
                    self._checkDetectorOverscanConfig(detectorConfig, integerDitherMode="NONE")
                else:
                    self._checkDetectorOverscanConfig(detectorConfig)

                for amp in detector:
                    ampConfig = detectorConfig.getOverscanAmpConfig(amp)
                    if detector.getName() == camera[1].getName():
                        self._checkOverscanConfig(ampConfig, doParallelOverscan=False)
                    else:
                        self._checkOverscanConfig(ampConfig)

            self._checkAnyOverscanConfig(config)

    def testAmpConfigOverrideOneDetectorOneAmp(self):
        camera = self._makeCamera()

        overscanAmpConfigOverride = OverscanAmpConfig()
        overscanAmpConfigOverride.serialOverscanConfig.fitType = "MEDIAN"
        overscanDetectorConfigOverride = OverscanDetectorConfig()
        overscanDetectorConfigOverride.ampRules["amp3"] = overscanAmpConfigOverride

        config = OverscanCameraConfig()
        config.detectorRules["detector0"] = overscanDetectorConfigOverride

        config = self._serializeAndReadConfig(config)

        for detector in camera:
            for amp in detector:
                detectorConfig = config.getOverscanDetectorConfig(detector)
                ampConfig = detectorConfig.getOverscanAmpConfig(amp)
                if detector.getName() == "detector0" and amp.getName() == "amp3":
                    self._checkOverscanConfig(ampConfig, serialFitType="MEDIAN")
                else:
                    self._checkOverscanConfig(ampConfig)

        self._checkAnyOverscanConfig(config)

    def testAmpConfigMd5(self):
        # Check a default detectorConfig
        detectorConfig1 = OverscanDetectorConfig()
        configMd51 = detectorConfig1.md5

        # Make sure copying it has the same hash.
        detectorConfig2 = copy.copy(detectorConfig1)
        configMd52 = detectorConfig2.md5

        self.assertEqual(configMd51, configMd52)

        # Make a new one with an amp override.
        overscanAmpConfigOverride = OverscanAmpConfig()
        overscanAmpConfigOverride.parallelOverscanConfig.fitType = "MEDIAN"

        detectorConfig3 = OverscanDetectorConfig()
        detectorConfig3.ampRules["amp2"] = overscanAmpConfigOverride

        self.assertNotEqual(detectorConfig3.md5, detectorConfig1.md5)

        # Override another amp but with the default.  This should
        # give the same answer because amp overrides that match the default
        # are not hashed.
        detectorConfig4 = copy.copy(detectorConfig3)
        detectorConfig4.ampRules["amp3"] = OverscanAmpConfig()

        self.assertEqual(detectorConfig4.md5, detectorConfig3.md5)

    def testAmpConfigBadAmps(self):
        detectorConfig = OverscanDetectorConfig()
        overscanAmpConfigOverride = OverscanAmpConfig()
        overscanAmpConfigOverride.maskAmpAsBad = True

        badAmps = detectorConfig.badAmpsToMask
        self.assertEqual(len(badAmps), 0)

        detectorConfig.ampRules["amp1"] = overscanAmpConfigOverride
        detectorConfig.ampRules["amp2"] = overscanAmpConfigOverride

        badAmps = detectorConfig.badAmpsToMask
        self.assertEqual(badAmps, ["amp1", "amp2"])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
