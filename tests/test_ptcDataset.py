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

import numpy as np
import copy

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

        self.ampNames = [amp.getName() for amp in self.flatExp1.getDetector().getAmplifiers()]
        self.dataset = PhotonTransferCurveDataset(self.ampNames, " ")

    def test_ptcDatset(self):
        localDataset = copy.copy(self.dataset)

        # Fill the set up with made up data.
        nSignalPoints = 5
        nSideCovMatrix = 2
        localDataset.ptcFitType = 'POLYNOMIAL'
        localDataset.badAmps = [localDataset.ampNames[0], localDataset.ampNames[1]]

        for ampName in localDataset.ampNames:
            localDataset.inputExpIdPairs[ampName] = np.repeat(1, nSignalPoints)
            localDataset.expIdMask[ampName] = [True, False, True, True, False, True, False, True, True,
                                               True, True, False, True, False, True]
            localDataset.expIdMask[ampName] = np.arange(nSignalPoints)
            localDataset.rawExpTimes[ampName] = np.arange(nSignalPoints)
            localDataset.rawMeans[ampName] = self.flux*np.arange(nSignalPoints)
            localDataset.rawVars[ampName] = self.noiseSq + self.c1*self.flux*np.arange(nSignalPoints)
            localDataset.photoCharge[ampName] = np.repeat(np.nan, nSignalPoints)
            localDataset.gain[ampName] = self.gain
            localDataset.gainErr[ampName] = 0.1
            localDataset.noise[ampName] = self.noiseSq
            localDataset.noiseErr[ampName] = 2.0
            localDataset.ptcFitPars[ampName] = np.array([10.0, 1.5, 1e-6])
            localDataset.ptcFitParsError[ampName] = np.array([1.0, 0.2, 1e-7])
            localDataset.ptcFitChiSq[ampName] = 1.0
            localDataset.covariances[ampName] = np.full((nSignalPoints, nSideCovMatrix, nSideCovMatrix),
                                                        105.0)
            localDataset.covariancesModel[ampName] = np.full((nSignalPoints, nSideCovMatrix, nSideCovMatrix),
                                                             100.0)
            localDataset.covariancesSqrtWeights[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                   nSideCovMatrix), 10.0)
            localDataset.aMatrix[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), 1e-6)
            localDataset.bMatrix[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), np.nan)
            localDataset.covariancesNoB[ampName] = np.full((nSignalPoints, nSideCovMatrix, nSideCovMatrix),
                                                           np.nan)
            localDataset.covariancesModelNoB[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                nSideCovMatrix), np.nan)
            localDataset.covariancesSqrtWeightsNoB[ampName] = np.full((nSignalPoints, nSideCovMatrix,
                                                                      nSideCovMatrix), np.nan)
            localDataset.aMatrixNoB[ampName] = np.full((nSideCovMatrix, nSideCovMatrix), np.nan)
            localDataset.finalVars[ampName] = self.noiseSq + self.c1*self.flux*np.arange(nSignalPoints)
            localDataset.finalModelVars[ampName] = np.repeat(100.0, nSignalPoints)
            localDataset.finalMeans[ampName] = self.flux*np.arange(nSignalPoints)

        filename = tempfile.mktemp()
        usedFilename = localDataset.writeFits(filename + ".fits")
        fromFits = PhotonTransferCurveDataset.readFits(usedFilename)
        self.assertEqual(localDataset, fromFits)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
