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
            doParallelOverscanCrosstalk=True,
            doParallelOverscan=True,
            serialFitType="MEDIAN_PER_ROW",
            parallelFitType="MEDIAN_PER_ROW",
    ):
        self.assertEqual(overscanAmpConfig.doSerialOverscan, doSerialOverscan)
        self.assertEqual(overscanAmpConfig.doParallelOverscanCrosstalk, doParallelOverscanCrosstalk)
        self.assertEqual(overscanAmpConfig.doParallelOverscan, doParallelOverscan)
        self.assertEqual(overscanAmpConfig.serialOverscanConfig.fitType, serialFitType)
        self.assertEqual(overscanAmpConfig.parallelOverscanConfig.fitType, parallelFitType)

    def _checkAnyOverscanConfig(
            self,
            config,
            doSerialOverscan=True,
            doParallelOverscan=True,
            doParallelOverscanCrosstalk=True,
    ):
        self.assertEqual(config.doAnySerialOverscan, doSerialOverscan)
        self.assertEqual(config.doAnyParallelOverscan, doParallelOverscan)
        self.assertEqual(config.doAnyParallelOverscanCrosstalk, doParallelOverscanCrosstalk)

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
            doParallelOverscanCrosstalk=False,
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
                    doParallelOverscanCrosstalk=False,
                )

        self._checkAnyOverscanConfig(
            config,
            doSerialOverscan=False,
            doParallelOverscan=False,
            doParallelOverscanCrosstalk=False,
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

        overscanAmpConfigOverride = OverscanAmpConfig(doParallelOverscanCrosstalk=False)
        overscanDetectorConfigOverride = OverscanDetectorConfig(defaultAmpConfig=overscanAmpConfigOverride)

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
                for amp in detector:
                    detectorConfig = config.getOverscanDetectorConfig(detector)
                    ampConfig = detectorConfig.getOverscanAmpConfig(amp)
                    if detector.getName() == camera[1].getName():
                        self._checkOverscanConfig(ampConfig, doParallelOverscanCrosstalk=False)
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


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
