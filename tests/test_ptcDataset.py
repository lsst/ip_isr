# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import unittest
import tempfile
import copy
import logging

import numpy as np

import lsst.utils.tests

from lsst.ip.isr import PhotonTransferCurveDataset
import lsst.ip.isr.isrMock as isrMock


class PtcDatasetCases(lsst.utils.tests.TestCase):
    """Test that write/read methods of PhotonTransferCurveDataset work
    """
    def setUp(self):

        self.flatMean = 2000
        self.readNoiseAdu = 10
        mockImageConfig = isrMock.IsrMock.ConfigClass()

        # flatDrop is not really relevant as we replace the data
        # but good to note it in case we change how this image is made
        mockImageConfig.flatDrop = 0.99999
        mockImageConfig.isTrimmed = True

        self.flatExp1 = isrMock.FlatMock(config=mockImageConfig).run()
        self.flatExp2 = self.flatExp1.clone()
        (shapeY, shapeX) = self.flatExp1.getDimensions()

        self.flatWidth = np.sqrt(self.flatMean) + self.readNoiseAdu

        self.rng1 = np.random.RandomState(1984)
        flatData1 = self.rng1.normal(self.flatMean, self.flatWidth, (shapeX, shapeY))
        self.rng2 = np.random.RandomState(666)
        flatData2 = self.rng2.normal(self.flatMean, self.flatWidth, (shapeX, shapeY))

        self.flatExp1.image.array[:] = flatData1
        self.flatExp2.image.array[:] = flatData2

        self.flux = 1000.  # ADU/sec
        self.gain = 1.5  # e-/ADU
        self.noiseSq = 5*self.gain  # 7.5 (e-)^2
        self.c1 = 1./self.gain
        self.timeVec = np.arange(1., 101., 5)
        self.k2NonLinearity = -5e-6
        # quadratic signal-chain non-linearity
        muVec = self.flux*self.timeVec + self.k2NonLinearity*self.timeVec**2

        self.ampNames = [amp.getName() for amp in self.flatExp1.getDetector().getAmplifiers()]
        self.dataset = PhotonTransferCurveDataset(self.ampNames, " ")  # pack raw data for fitting
        self.covariancesSqrtWeights = {}
        for ampName in self.ampNames:  # just the expTimes and means here - vars vary per function
            self.dataset.rawExpTimes[ampName] = self.timeVec
            self.dataset.rawMeans[ampName] = muVec
            self.covariancesSqrtWeights[ampName] = []

    def _checkTypes(self, ptcDataset):
        """Check that all the types are correct for a ptc dataset."""
        for ampName in ptcDataset.ampNames:
            self.assertIsInstance(ptcDataset.expIdMask[ampName], np.ndarray)
            self.assertEqual(ptcDataset.expIdMask[ampName].dtype, bool)
            self.assertIsInstance(ptcDataset.rawExpTimes[ampName], np.ndarray)
            self.assertEqual(ptcDataset.rawExpTimes[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.rawMeans[ampName], np.ndarray)
            self.assertEqual(ptcDataset.rawMeans[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.rawVars[ampName], np.ndarray)
            self.assertEqual(ptcDataset.rawVars[ampName].dtype, np.float64)
            self.assertEqual(ptcDataset.rowMeanVariance[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.noiseList[ampName], np.ndarray)
            self.assertEqual(ptcDataset.noiseList[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.gain[ampName], float)
            self.assertIsInstance(ptcDataset.gainErr[ampName], float)
            self.assertIsInstance(ptcDataset.noise[ampName], float)
            self.assertIsInstance(ptcDataset.noiseErr[ampName], float)
            self.assertIsInstance(ptcDataset.histVars[ampName], np.ndarray)
            self.assertEqual(ptcDataset.histVars[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.histChi2Dofs[ampName], np.ndarray)
            self.assertEqual(ptcDataset.histChi2Dofs[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.kspValues[ampName], np.ndarray)
            self.assertEqual(ptcDataset.kspValues[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.ptcFitPars[ampName], np.ndarray)
            self.assertEqual(ptcDataset.ptcFitPars[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.ptcFitParsError[ampName], np.ndarray)
            self.assertEqual(ptcDataset.ptcFitParsError[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.ptcFitChiSq[ampName], float)
            self.assertIsInstance(ptcDataset.ptcTurnoff[ampName], float)
            self.assertIsInstance(ptcDataset.ptcTurnoffSamplingError[ampName], float)
            self.assertIsInstance(ptcDataset.covariances[ampName], np.ndarray)
            self.assertEqual(ptcDataset.covariances[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.covariancesModel[ampName], np.ndarray)
            self.assertEqual(ptcDataset.covariancesModel[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.covariancesSqrtWeights[ampName], np.ndarray)
            self.assertEqual(ptcDataset.covariancesSqrtWeights[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.aMatrix[ampName], np.ndarray)
            self.assertEqual(ptcDataset.aMatrix[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.bMatrix[ampName], np.ndarray)
            self.assertEqual(ptcDataset.bMatrix[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.noiseMatrix[ampName], np.ndarray)
            self.assertEqual(ptcDataset.noiseMatrix[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.covariancesModelNoB[ampName], np.ndarray)
            self.assertEqual(ptcDataset.covariancesModelNoB[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.aMatrixNoB[ampName], np.ndarray)
            self.assertEqual(ptcDataset.aMatrixNoB[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.noiseMatrixNoB[ampName], np.ndarray)
            self.assertEqual(ptcDataset.noiseMatrixNoB[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.finalVars[ampName], np.ndarray)
            self.assertEqual(ptcDataset.finalVars[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.finalModelVars[ampName], np.ndarray)
            self.assertEqual(ptcDataset.finalModelVars[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.finalMeans[ampName], np.ndarray)
            self.assertEqual(ptcDataset.finalMeans[ampName].dtype, np.float64)
            self.assertIsInstance(ptcDataset.photoCharges[ampName], np.ndarray)
            self.assertEqual(ptcDataset.photoCharges[ampName].dtype, np.float64)

        for key, value in ptcDataset.auxValues.items():
            self.assertIsInstance(value, np.ndarray)
            self.assertEqual(value.dtype, np.float64)

    def test_emptyPtcDataset(self):
        """Test an empty PTC dataset."""
        emptyDataset = PhotonTransferCurveDataset(
            self.ampNames,
            ptcFitType="PARTIAL",
        )
        self._checkTypes(emptyDataset)

        with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
            usedFilename = emptyDataset.writeText(f.name)
            fromText = PhotonTransferCurveDataset.readText(usedFilename)
        self.assertEqual(emptyDataset, fromText)
        self._checkTypes(emptyDataset)

        with tempfile.NamedTemporaryFile(suffix=".fits") as f:
            usedFilename = emptyDataset.writeFits(f.name)
            fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
        self.assertEqual(emptyDataset, fromFits)
        self._checkTypes(emptyDataset)

    def test_partialPtcDataset(self):
        """Test of a partial PTC dataset."""
        # Fill the dataset with made up data.
        nSideCovMatrix = 2
        nSideCovMatrixFullCovFit = 2

        partialDataset = PhotonTransferCurveDataset(
            self.ampNames,
            ptcFitType="PARTIAL",
            covMatrixSide=nSideCovMatrix,
            covMatrixSideFullCovFit=nSideCovMatrixFullCovFit
        )
        self._checkTypes(partialDataset)

        for ampName in partialDataset.ampNames:
            partialDataset.setAmpValuesPartialDataset(
                ampName,
                inputExpIdPair=(10, 11),
                rawExpTime=10.0,
                rawMean=10.0,
                rawVar=10.0,
            )

        for useAuxValues in [False, True]:
            if useAuxValues:
                partialDataset.setAuxValuesPartialDataset(
                    {
                        "CCOBCURR": 1.0,
                        "CCDTEMP": 0.0,
                    }
                )
            self._checkTypes(partialDataset)

            with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
                usedFilename = partialDataset.writeText(f.name)
                fromText = PhotonTransferCurveDataset.readText(usedFilename)
            self.assertEqual(fromText, partialDataset)
            self._checkTypes(fromText)

            with tempfile.NamedTemporaryFile(suffix=".fits") as f:
                usedFilename = partialDataset.writeFits(f.name)
                fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
            self.assertEqual(fromFits, partialDataset)
            self._checkTypes(fromFits)

    def test_ptcDatset(self):
        """Test of a full PTC dataset."""
        # Fill the dataset with made up data.
        nSignalPoints = 5
        nSideCovMatrixInput = 3  # Size of measured covariances

        for nSideCovMatrixFullCovFitInput in np.arange(1, nSideCovMatrixInput + 2):
            for fitType in ['POLYNOMIAL', 'EXPAPPROXIMATION', 'FULLCOVARIANCE']:
                localDataset = PhotonTransferCurveDataset(
                    self.ampNames,
                    ptcFitType=fitType,
                    covMatrixSide=nSideCovMatrixInput,
                    covMatrixSideFullCovFit=nSideCovMatrixFullCovFitInput,
                )
                nSideCovMatrix = localDataset.covMatrixSide
                nSideCovMatrixFullCovFit = localDataset.covMatrixSideFullCovFit
                localDataset.badAmps = [localDataset.ampNames[0], localDataset.ampNames[1]]
                for ampName in localDataset.ampNames:

                    localDataset.inputExpIdPairs[ampName] = [(1, 2)]*nSignalPoints
                    localDataset.expIdMask[ampName] = np.ones(nSignalPoints, dtype=bool)
                    localDataset.expIdMask[ampName][1] = False
                    localDataset.rawExpTimes[ampName] = np.arange(nSignalPoints, dtype=np.float64)
                    localDataset.rawMeans[ampName] = self.flux*np.arange(nSignalPoints)
                    localDataset.rawVars[ampName] = self.c1*self.flux*np.arange(nSignalPoints)
                    localDataset.photoCharges[ampName] = np.full(nSignalPoints, np.nan)
                    localDataset.gain[ampName] = self.gain
                    localDataset.gainErr[ampName] = 0.1
                    localDataset.noise[ampName] = self.noiseSq
                    localDataset.noiseErr[ampName] = 2.0
                    localDataset.histVars[ampName] = localDataset.rawVars[ampName]
                    localDataset.histChi2Dofs[ampName] = np.full(nSignalPoints, 1.0)
                    localDataset.kspValues[ampName] = np.full(nSignalPoints, 0.5)

                    localDataset.finalVars[ampName] = self.c1*self.flux*np.arange(nSignalPoints)
                    localDataset.finalModelVars[ampName] = np.full(nSignalPoints, 100.0)
                    localDataset.finalMeans[ampName] = self.flux*np.arange(nSignalPoints)

                    if fitType in ['POLYNOMIAL', 'EXPAPPROXIMATION', ]:
                        localDataset.ptcFitPars[ampName] = np.array([10.0, 1.5, 1e-6])
                        localDataset.ptcFitParsError[ampName] = np.array([1.0, 0.2, 1e-7])
                        localDataset.ptcFitChiSq[ampName] = 1.0
                        localDataset.ptcTurnoff[ampName] = localDataset.rawMeans[ampName][-1]
                        localDataset.ptcTurnoffSamplingError[ampName] = localDataset.ptcTurnoff[ampName]/100.

                        localDataset.covariances[ampName] = np.full(
                            (nSignalPoints, nSideCovMatrix, nSideCovMatrix), 105.0)
                        localDataset.covariancesModel[ampName] = np.full(
                            (nSignalPoints, nSideCovMatrixFullCovFit, nSideCovMatrixFullCovFit), np.nan)
                        localDataset.covariancesSqrtWeights[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                               nSideCovMatrix), 10.0)
                        localDataset.aMatrix[ampName] = np.full((nSideCovMatrixFullCovFit,
                                                                nSideCovMatrixFullCovFit), np.nan)
                        localDataset.bMatrix[ampName] = np.full((nSideCovMatrixFullCovFit,
                                                                nSideCovMatrixFullCovFit), np.nan)
                        localDataset.noiseMatrix[ampName] = np.full((nSideCovMatrixFullCovFit,
                                                                    nSideCovMatrixFullCovFit), np.nan)
                        localDataset.covariancesModelNoB[ampName] = np.full((nSignalPoints,
                                                                            nSideCovMatrixFullCovFit,
                                                                            nSideCovMatrixFullCovFit), np.nan)
                        localDataset.aMatrixNoB[ampName] = np.full(
                            (nSideCovMatrixFullCovFit, nSideCovMatrixFullCovFit), np.nan)
                        localDataset.noiseMatrixNoB[ampName] = np.full(
                            (nSideCovMatrixFullCovFit, nSideCovMatrixFullCovFit), np.nan)

                    if localDataset.ptcFitType in ['FULLCOVARIANCE', ]:
                        localDataset.ptcFitPars[ampName] = np.array([np.nan, np.nan])
                        localDataset.ptcFitParsError[ampName] = np.array([np.nan, np.nan])
                        localDataset.ptcFitChiSq[ampName] = np.nan
                        localDataset.ptcTurnoff[ampName] = np.nan
                        localDataset.ptcTurnoffSamplingError[ampName] = np.nan

                        localDataset.covariances[ampName] = np.full(
                            (nSignalPoints, nSideCovMatrix, nSideCovMatrix), 105.0)
                        localDataset.covariancesModel[ampName] = np.full(
                            (nSignalPoints, nSideCovMatrixFullCovFit, nSideCovMatrixFullCovFit), 100.0)
                        localDataset.covariancesSqrtWeights[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                               nSideCovMatrix), 10.0)
                        localDataset.aMatrix[ampName] = np.full((nSideCovMatrixFullCovFit,
                                                                nSideCovMatrixFullCovFit), 1e-6)
                        localDataset.bMatrix[ampName] = np.full((nSideCovMatrixFullCovFit,
                                                                nSideCovMatrixFullCovFit), 1e-7)
                        localDataset.noiseMatrix[ampName] = np.full((nSideCovMatrixFullCovFit,
                                                                    nSideCovMatrixFullCovFit), 3.0)
                        localDataset.covariancesModelNoB[ampName] = np.full((nSignalPoints,
                                                                            nSideCovMatrixFullCovFit,
                                                                            nSideCovMatrixFullCovFit), 15.0)
                        localDataset.aMatrixNoB[ampName] = np.full(
                            (nSideCovMatrixFullCovFit, nSideCovMatrixFullCovFit), 2e-6)
                        localDataset.noiseMatrixNoB[ampName] = np.full(
                            (nSideCovMatrixFullCovFit, nSideCovMatrixFullCovFit), 3.0)

                for useAuxValues in [False, True]:
                    if useAuxValues:
                        localDataset.auxValues = {
                            "CCOBCURR": np.ones(nSignalPoints),
                            "CCDTEMP": np.zeros(nSignalPoints),
                        }

                    self._checkTypes(localDataset)
                    with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
                        usedFilename = localDataset.writeText(f.name)
                        fromText = PhotonTransferCurveDataset.readText(usedFilename)
                    self.assertEqual(fromText, localDataset)
                    self._checkTypes(fromText)

                    with tempfile.NamedTemporaryFile(suffix=".fits") as f:
                        usedFilename = localDataset.writeFits(f.name)
                        fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
                    self.assertEqual(fromFits, localDataset)
                    self._checkTypes(fromFits)

    def test_getExpIdsUsed(self):
        localDataset = copy.copy(self.dataset)

        for pair in [(12, 34), (56, 78), (90, 10)]:
            localDataset.inputExpIdPairs["C:0,0"].append(pair)
        localDataset.expIdMask["C:0,0"] = np.array([True, False, True])
        self.assertTrue(np.all(localDataset.getExpIdsUsed("C:0,0") == [(12, 34), (90, 10)]))

        localDataset.expIdMask["C:0,0"] = np.array([True, False, True, True])  # wrong length now
        with self.assertRaises(AssertionError):
            localDataset.getExpIdsUsed("C:0,0")

    def test_getGoodAmps(self):
        dataset = self.dataset

        self.assertTrue(dataset.ampNames == self.ampNames)
        dataset.badAmps.append("C:0,1")
        self.assertTrue(dataset.getGoodAmps() == [amp for amp in self.ampNames if amp != "C:0,1"])

    def test_ptcDataset_pre_dm38309(self):
        """Test for PTC datasets created by cpSolvePtcTask prior to DM-38309.
        """
        localDataset = copy.copy(self.dataset)

        for pair in [[(12, 34)], [(56, 78)], [(90, 10)]]:
            localDataset.inputExpIdPairs["C:0,0"].append(pair)
        localDataset.expIdMask["C:0,0"] = np.array([True, False, True])

        with self.assertLogs("lsst.ip.isr.calibType", logging.WARNING) as cm:
            used = localDataset.getExpIdsUsed("C:0,0")
        self.assertIn("PTC file was written incorrectly", cm.output[0])

        self.assertTrue(np.all(used == [(12, 34), (90, 10)]))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
