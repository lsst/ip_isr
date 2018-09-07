#
# LSST Data Management System
# Copyright 2016 AURA/LSST.
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
import abc

import numpy as np

from lsst.pipe.base import Struct
from .applyLookupTable import applyLookupTable

__all__ = ["LinearizeBase", "LinearizeLookupTable", "LinearizeSquared"]


class LinearizeBase(metaclass=abc.ABCMeta):
    """Abstract base class functor for correcting non-linearity

    Subclasses must define __call__ and set class variable LinearityType to a string
    that will be used for linearity type in AmpInfoCatalog
    """
    LinearityType = None  # linearity type, a string used for AmpInfoCatalogs

    @abc.abstractmethod
    def __call__(self, image, detector, log=None):
        """Correct non-linearity

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            image to be corrected (an lsst.afw.image.Image)
        detector  detector information (an instance of lsst::afw::cameraGeom::Detector)
        log :
            logger (an lsst.log.Log), or None to disable logging;
            a warning is logged if amplifiers are skipped or other worrisome events occur

        Returns
        -------
        result : `Struct`
            return an lsst.pipe.base.Struct containing at least the following fields:

            - ``numAmps`` : number of amplifiers found
            - ``numLinearized`` :  number of amplifiers linearized

        Raises
        ------
        RuntimeError :
            if the linearity type is wrong
        a subclass of Exception
            if linearization fails for any other reason
        """
        pass

    def checkLinearityType(self, detector):
        """Verify that the linearity type is correct for this detector

        warning : only checks the first record of the amp info catalog

        Parameters
        ----------
        detector :  detector information (an instance of lsst::afw::cameraGeom::Detector)

        Raises
        ------
        RuntimeError
            if anything doesn't match
        """
        ampInfoType = detector.getAmpInfoCatalog()[0].getLinearityType()
        if self.LinearityType != ampInfoType:
            raise RuntimeError("Linearity types don't match: %s != %s" % (self.LinearityType, ampInfoType))


class LinearizeLookupTable(LinearizeBase):
    """Correct non-linearity with a persisted lookup table

    Notes
    -----
    for each i,j of image:

    .. code-block :: none

        rowInd = int(c0)
        colInd = int(c1 + uncorrImage[i,j])
        corrImage[i,j] = uncorrImage[i,j] + table[rowInd, colInd]

    where c0, c1 are collimation coefficients from the AmpInfoTable of the detector:
    - c0: row index; used to identify which row of the table to use (typically one per amplifier,
    though one can have multiple amplifiers use the same table)
    - c1: column index offset; added to the uncorrected image value before truncation;
    this supports tables that can handle negative image values; also, if the c1 ends with .5
    then the nearest index is used instead of truncating to the next smaller index

    In order to keep related data together, the coefficients are persisted along with the table.
    """
    LinearityType = "LookupTable"

    def __init__(self, table, detector):
        """Construct a LinearizeLookupTable

        Parameters
        ----------
        table :
            To avoid copying the table the last index should vary fastest (numpy default "C" order)
            lookup table; a 2-dimensional array of floats:

            - one row for each row index (value of coef[0] in the amp info catalog)
            - one column for each image value

        detector : `lsst::afw::cameraGeom::Detector`
            detector information (an instance of lsst::afw::cameraGeom::Detector);
            the name, serial, and amplifier linearization type and coefficients are saved

        Raises
        ------
        RuntimeError
            if table is not 2-dimensional,
            table has fewer columns than rows (indicating that the indices are swapped),
            or if any row index (linearity coefficient 0) is out of range
        """
        LinearizeBase.__init__(self)

        self._table = np.array(table, order="C")
        if len(table.shape) != 2:
            raise RuntimeError("table shape = %s; must have two dimensions" % (table.shape,))
        if table.shape[1] < table.shape[0]:
            raise RuntimeError("table shape = %s; indices are switched" % (table.shape,))

        self._detectorName = detector.getName()
        self._detectorSerial = detector.getSerial()
        self.checkLinearityType(detector)
        ampInfoCat = detector.getAmpInfoCatalog()
        rowIndList = []
        colIndOffsetList = []
        numTableRows = table.shape[0]
        for ampInfo in ampInfoCat:
            rowInd, colIndOffset = ampInfo.getLinearityCoeffs()[0:2]
            rowInd = int(rowInd)
            if rowInd < 0 or rowInd >= numTableRows:
                raise RuntimeError("Amplifier %s has rowInd=%s not in range[0, %s)" %
                                   (ampInfo.getName(), rowInd, numTableRows))
            rowIndList.append(int(rowInd))
            colIndOffsetList.append(colIndOffset)
        self._rowIndArr = np.array(rowIndList, dtype=int)
        self._colIndOffsetArr = np.array(colIndOffsetList)

    def __call__(self, image, detector, log=None):
        """Correct for non-linearity

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            image to be corrected (an lsst.afw.image.Image)
        detector : `lsst.afw.cameraGeom.Detector`
            detector info about image (an lsst.afw.cameraGeom.Detector);
            the name, serial and number of amplifiers must match persisted data;
            the bbox from each amplifier is read;
            the linearization coefficients are ignored in favor of the persisted values
        log : `lsst.log.Log`
            logger (an lsst.log.Log), or None to disable logging;
            a warning is logged if any pixels are out of range of their lookup table

        Returns
        -------
        Struct : `lsst.pipe.base.Struct`
            return an lsst.pipe.base.Struct containing:

            - ``numAmps`` : number of amplifiers found
            - ``numLinearized`` :  number of amplifiers linearized
            (always equal to numAmps for this linearizer)
            - ``numOutOfRange`` :  number of pixels out of range of their lookup table
            (summed across all amps)

        Raises
        ------
        RuntimeError
            if the linearity type is wrong or if the detector name, serial
            or number of amplifiers does not match the saved data
        """
        self.checkDetector(detector)
        ampInfoCat = detector.getAmpInfoCatalog()
        numOutOfRange = 0
        for ampInfo, rowInd, colIndOffset in zip(ampInfoCat, self._rowIndArr, self._colIndOffsetArr):
            bbox = ampInfo.getBBox()
            ampView = image.Factory(image, bbox)
            tableRow = self._table[rowInd, :]
            numOutOfRange += applyLookupTable(ampView, tableRow, colIndOffset)

        if numOutOfRange > 0 and log is not None:
            log.warn("%s pixels of detector \"%s\" were out of range of the linearization table",
                     numOutOfRange, detector.getName())
        numAmps = len(ampInfoCat)
        return Struct(
            numAmps=numAmps,
            numLinearized=numAmps,
            numOutOfRange=numOutOfRange,
        )

    def checkDetector(self, detector):
        """Check detector name and serial number, ampInfo table length and linearity type

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            detector info about image (an lsst.afw.cameraGeom.Detector);

        Raises
        ------
        RuntimeError
            if anything doesn't match
        """
        if self._detectorName != detector.getName():
            raise RuntimeError("Detector names don't match: %s != %s" %
                               (self._detectorName, detector.getName()))
        if self._detectorSerial != detector.getSerial():
            raise RuntimeError("Detector serial numbers don't match: %s != %s" %
                               (self._detectorSerial, detector.getSerial()))

        numAmps = len(detector.getAmpInfoCatalog())
        if numAmps != len(self._rowIndArr):
            raise RuntimeError("Detector number of amps = %s does not match saved value %s" %
                               (numAmps, len(self._rowIndArr)))
        self.checkLinearityType(detector)


class LinearizeSquared(LinearizeBase):
    """Correct non-linearity with a squared model

    corrImage = uncorrImage + c0*uncorrImage^2

    where c0 is linearity coefficient 0 in the AmpInfoCatalog of the detector
    """
    LinearityType = "Squared"

    def __call__(self, image, detector, log=None):
        """Correct for non-linearity

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            image to be corrected (an lsst.afw.image.Image)
        detector : `lsst.afw.cameraGeom.Detector`
            detector info about image (an lsst.afw.cameraGeom.Detector)
        log : `lsst.log.Log`
            logger (an lsst.log.Log), or None to disable logging;
            a warning is logged if any amplifiers are skipped because the square coefficient is 0

        Returns
        -------
        Struct : `lsst.pipe.base.Struct`
            return an lsst.pipe.base.Struct containing at least the following fields:

            - nAmps number of amplifiers found
            - nLinearized  number of amplifiers linearized

        Raises
        ------
        RuntimeError
            if the linearity type is wrong
        """
        self.checkLinearityType(detector)
        ampInfoCat = detector.getAmpInfoCatalog()
        numLinearized = 0
        for ampInfo in ampInfoCat:
            sqCoeff = ampInfo.getLinearityCoeffs()[0]
            if sqCoeff != 0:
                bbox = ampInfo.getBBox()
                ampArr = image.Factory(image, bbox).getArray()
                ampArr *= (1 + sqCoeff*ampArr)
                numLinearized += 1

        numAmps = len(ampInfoCat)
        if numAmps > numLinearized and log is not None:
            log.warn("%s of %s amps in detector \"%s\" were not linearized (coefficient = 0)",
                     numAmps - numLinearized, numAmps, detector.getName())
        return Struct(
            numAmps=numAmps,
            numLinearized=numLinearized,
        )
