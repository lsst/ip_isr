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

    def test_emptyPtcDataset(self):
        """Test an empty PTC dataset."""
        emptyDataset = PhotonTransferCurveDataset(
            self.ampNames,
            ptcFitType="PARTIAL",
        )

        with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
            usedFilename = emptyDataset.writeText(f.name)
            fromText = PhotonTransferCurveDataset.readText(usedFilename)
        self.assertEqual(emptyDataset, fromText)

        with tempfile.NamedTemporaryFile(suffix=".fits") as f:
            usedFilename = emptyDataset.writeFits(f.name)
            fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
        self.assertEqual(emptyDataset, fromFits)

    def test_partialPtcDataset(self):
        """Test of a partial PTC dataset."""
        # Fill the dataset with made up data.
        nSideCovMatrix = 2

        partialDataset = PhotonTransferCurveDataset(
            self.ampNames,
            ptcFitType="PARTIAL",
            covMatrixSide=nSideCovMatrix
        )

        for ampName in partialDataset.ampNames:
            partialDataset.setAmpValuesPartialDataset(
                ampName,
                inputExpIdPair=(10, 11),
                rawExpTime=10.0,
                rawMean=10.0,
                rawVar=10.0,
            )

        with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
            usedFilename = partialDataset.writeText(f.name)
            fromText = PhotonTransferCurveDataset.readText(usedFilename)
        self.assertEqual(partialDataset, fromText)

        with tempfile.NamedTemporaryFile(suffix=".fits") as f:
            usedFilename = partialDataset.writeFits(f.name)
            fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
        self.assertEqual(partialDataset, fromFits)

    def test_ptcDatset(self):
        """Test of a full PTC dataset."""
        # Fill the dataset with made up data.
        nSignalPoints = 5
        nSideCovMatrix = 2
        for fitType in ['POLYNOMIAL', 'EXPAPPROXIMATION', 'FULLCOVARIANCE']:
            localDataset = PhotonTransferCurveDataset(
                self.ampNames,
                ptcFitType=fitType,
                covMatrixSide=nSideCovMatrix,
            )
            localDataset.badAmps = [localDataset.ampNames[0], localDataset.ampNames[1]]
            for ampName in localDataset.ampNames:

                localDataset.inputExpIdPairs[ampName] = [(1, 2)]*nSignalPoints
                localDataset.expIdMask[ampName] = np.ones(nSignalPoints, dtype=bool)
                localDataset.expIdMask[ampName][1] = False
                localDataset.rawExpTimes[ampName] = np.arange(nSignalPoints)
                localDataset.rawMeans[ampName] = self.flux*np.arange(nSignalPoints)
                localDataset.rawVars[ampName] = self.c1*self.flux*np.arange(nSignalPoints)
                localDataset.photoCharges[ampName] = np.full(nSignalPoints, np.nan)
                localDataset.gain[ampName] = self.gain
                localDataset.gainErr[ampName] = 0.1
                localDataset.noise[ampName] = self.noiseSq
                localDataset.noiseErr[ampName] = 2.0

                localDataset.finalVars[ampName] = self.c1*self.flux*np.arange(nSignalPoints)
                localDataset.finalModelVars[ampName] = np.full(nSignalPoints, 100.0)
                localDataset.finalMeans[ampName] = self.flux*np.arange(nSignalPoints)

                if fitType in ['POLYNOMIAL', 'EXPAPPROXIMATION', ]:
                    localDataset.ptcFitPars[ampName] = np.array([10.0, 1.5, 1e-6])
                    localDataset.ptcFitParsError[ampName] = np.array([1.0, 0.2, 1e-7])
                    localDataset.ptcFitChiSq[ampName] = 1.0
                    localDataset.ptcTurnoff[ampName] = localDataset.rawMeans[ampName][-1]

                    localDataset.covariances[ampName] = np.full(
                        (nSignalPoints, nSideCovMatrix, nSideCovMatrix), 105.0)
                    localDataset.covariancesModel[ampName] = np.full(
                        (nSignalPoints, nSideCovMatrix, nSideCovMatrix), np.nan)
                    localDataset.covariancesSqrtWeights[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                           nSideCovMatrix), 10.0)
                    localDataset.aMatrix[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), np.nan)
                    localDataset.bMatrix[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), np.nan)
                    localDataset.covariancesModelNoB[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                        nSideCovMatrix), np.nan)
                    localDataset.aMatrixNoB[ampName] = np.full(
                        (nSideCovMatrix, nSideCovMatrix), np.nan)

                if localDataset.ptcFitType in ['FULLCOVARIANCE', ]:
                    localDataset.ptcFitPars[ampName] = np.array([np.nan, np.nan])
                    localDataset.ptcFitParsError[ampName] = np.array([np.nan, np.nan])
                    localDataset.ptcFitChiSq[ampName] = np.array([np.nan, np.nan])
                    localDataset.ptcTurnoff[ampName] = np.array([np.nan, np.nan])

                    localDataset.covariances[ampName] = np.full(
                        (nSignalPoints, nSideCovMatrix, nSideCovMatrix), 105.0)
                    localDataset.covariancesModel[ampName] = np.full(
                        (nSignalPoints, nSideCovMatrix, nSideCovMatrix), 100.0)
                    localDataset.covariancesSqrtWeights[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                           nSideCovMatrix), 10.0)
                    localDataset.aMatrix[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), 1e-6)
                    localDataset.bMatrix[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), 1e-7)
                    localDataset.covariancesModelNoB[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                        nSideCovMatrix), 15.0)
                    localDataset.aMatrixNoB[ampName] = np.full(
                        (nSideCovMatrix, nSideCovMatrix), 2e-6)

            with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
                usedFilename = localDataset.writeText(f.name)
                fromText = PhotonTransferCurveDataset.readText(usedFilename)
            self.assertEqual(localDataset, fromText)

            with tempfile.NamedTemporaryFile(suffix=".fits") as f:
                usedFilename = localDataset.writeFits(f.name)
                fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
            self.assertEqual(localDataset, fromFits)

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

        with self.assertWarnsRegex(RuntimeWarning, "PTC file was written incorrectly"):
            used = localDataset.getExpIdsUsed("C:0,0")

        self.assertTrue(np.all(used == [(12, 34), (90, 10)]))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
