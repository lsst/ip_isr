#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
"""
Define dataset class for MeasurePhotonTransferCurve task
"""

__all__ = ['PhotonTransferCurveDataset']

import numpy as np
import math
from astropy.table import Table

from lsst.ip.isr import IsrCalib


class PhotonTransferCurveDataset(IsrCalib):
    """A simple class to hold the output data from the PTC task.

    The dataset is made up of a dictionary for each item, keyed by the
    amplifiers' names, which much be supplied at construction time.
    New items cannot be added to the class to save accidentally saving to the
    wrong property, and the class can be frozen if desired.
    inputExpIdPairs records the exposures used to produce the data.
    When fitPtc() or fitCovariancesAstier() is run, a mask is built up, which
    is by definition always the same length as inputExpIdPairs, rawExpTimes,
    rawMeans and rawVars, and is a list of bools, which are incrementally set
    to False as points are discarded from the fits.
    PTC fit parameters for polynomials are stored in a list in ascending order
    of polynomial term, i.e. par[0]*x^0 + par[1]*x + par[2]*x^2 etc
    with the length of the list corresponding to the order of the polynomial
    plus one.

    Parameters
    ----------
    ampNames : `list`
        List with the names of the amplifiers of the detector at hand.
    ptcFitType : `str`, optional
        Type of model fitted to the PTC: "POLYNOMIAL", "EXPAPPROXIMATION",
        or "FULLCOVARIANCE".
    covMatrixSide : `int`, optional
        Maximum lag of measured covariances (size of square covariance
        matrices).
    covMatrixSideFullCovFit : `int, optional
        Maximum covariances lag for FULLCOVARIANCE fit. It should be less or
        equal than covMatrixSide.
    kwargs : `dict`, optional
        Other keyword arguments to pass to the parent init.

    Notes
    -----
    The stored attributes are:

    badAmps : `list` [`str`]
        List with bad amplifiers names.
    inputExpIdPairs : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the input exposures IDs.
    expIdMask : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the mask produced after
        outlier rejection. The mask produced by the "FULLCOVARIANCE"
        option may differ from the one produced in the other two PTC
        fit types.
    rawExpTimes : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the unmasked exposure times.
    rawMeans : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the unmasked average of the
        means of the exposures in each flat pair.
    rawVars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the variance of the
        difference image of the exposures in each flat pair.
    histVars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the variance of the
        difference image of the exposures in each flat pair estimated
        by fitting a Gaussian model.
    histChi2Dofs : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the chi-squared per degree
        of freedom fitting the difference image to a Gaussian model.
    kspValues : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the KS test p-value from
        fitting the difference image to a Gaussian model.
    gain : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the fitted gains.
    gainErr : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the errors on the
        fitted gains.
    noise : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the fitted noise.
    noiseErr : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the errors on the fitted
        noise.
    ptcFitPars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the fitted parameters of the
        PTC model for ptcFitTye in ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcFitParsError : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the errors on the fitted
        parameters of the PTC model for ptcFitTye in
        ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcFitChiSq : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the reduced chi squared
        of the fit for ptcFitTye in ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcTurnoff : `dict` [`str, `float`]
        Flux value (in ADU) where the variance of the PTC curve starts
        decreasing consistently.
    covariances : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing a list of measured
        covariances per mean flux.
    covariancesModel : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containinging covariances model
        (Eq. 20 of Astier+19) per mean flux.
    covariancesSqrtWeights : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containinging sqrt. of covariances
        weights.
    aMatrix : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "a" parameters from
        the model in Eq. 20 of Astier+19.
    bMatrix : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "b" parameters from
        the model in Eq. 20 of Astier+19.
    noiseMatrix : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "noise" parameters from
        the model in Eq. 20 of Astier+19.
    covariancesModelNoB : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing covariances model
        (with 'b'=0 in Eq. 20 of Astier+19)
        per mean flux.
    aMatrixNoB : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "a" parameters from the
        model in Eq. 20 of Astier+19
        (and 'b' = 0).
    noiseMatrixNoB : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "noise" parameters from
        the model in Eq. 20 of Astier+19, with 'b' = 0.
    finalVars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the masked variance of the
        difference image of each flat
        pair. If needed, each array will be right-padded with
        np.nan to match the length of rawExpTimes.
    finalModelVars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the masked modeled
        variance of the difference image of each flat pair. If needed, each
        array will be right-padded with np.nan to match the length of
        rawExpTimes.
    finalMeans : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the masked average of the
        means of the exposures in each flat pair. If needed, each array
        will be right-padded with np.nan to match the length of
        rawExpTimes.
    photoCharges : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the integrated photocharge
        for linearity calibration.
    auxValues : `dict`, [`str`, `np.ndarray`]
        Dictionary of per-detector auxiliary header values that can be used
        for PTC, linearity computation.

    Version 1.1 adds the `ptcTurnoff` attribute.
    Version 1.2 adds the `histVars`, `histChi2Dofs`, and `kspValues`
    attributes.
    Version 1.3 adds the `noiseMatrix` and `noiseMatrixNoB` attributes.
    Version 1.4 adds the `auxValues` attribute.
    Version 1.5 adds the `covMatrixSideFullCovFit` attribute.
    """

    _OBSTYPE = 'PTC'
    _SCHEMA = 'Gen3 Photon Transfer Curve'
    _VERSION = 1.5

    def __init__(self, ampNames=[], ptcFitType=None, covMatrixSide=1,
                 covMatrixSideFullCovFit=None, **kwargs):
        self.ptcFitType = ptcFitType
        self.ampNames = ampNames
        self.covMatrixSide = covMatrixSide
        if covMatrixSideFullCovFit is None:
            self.covMatrixSideFullCovFit = covMatrixSide
        else:
            self.covMatrixSideFullCovFit = covMatrixSideFullCovFit

        self.badAmps = []

        self.inputExpIdPairs = {ampName: [] for ampName in ampNames}
        self.expIdMask = {ampName: np.array([], dtype=bool) for ampName in ampNames}
        self.rawExpTimes = {ampName: np.array([]) for ampName in ampNames}
        self.rawMeans = {ampName: np.array([]) for ampName in ampNames}
        self.rawVars = {ampName: np.array([]) for ampName in ampNames}
        self.photoCharges = {ampName: np.array([]) for ampName in ampNames}

        self.gain = {ampName: np.nan for ampName in ampNames}
        self.gainErr = {ampName: np.nan for ampName in ampNames}
        self.noise = {ampName: np.nan for ampName in ampNames}
        self.noiseErr = {ampName: np.nan for ampName in ampNames}

        self.histVars = {ampName: np.array([]) for ampName in ampNames}
        self.histChi2Dofs = {ampName: np.array([]) for ampName in ampNames}
        self.kspValues = {ampName: np.array([]) for ampName in ampNames}

        self.ptcFitPars = {ampName: np.array([]) for ampName in ampNames}
        self.ptcFitParsError = {ampName: np.array([]) for ampName in ampNames}
        self.ptcFitChiSq = {ampName: np.nan for ampName in ampNames}
        self.ptcTurnoff = {ampName: np.nan for ampName in ampNames}

        self.covariances = {ampName: np.array([]) for ampName in ampNames}
        self.covariancesModel = {ampName: np.array([]) for ampName in ampNames}
        self.covariancesSqrtWeights = {ampName: np.array([]) for ampName in ampNames}
        self.aMatrix = {ampName: np.array([]) for ampName in ampNames}
        self.bMatrix = {ampName: np.array([]) for ampName in ampNames}
        self.noiseMatrix = {ampName: np.array([]) for ampName in ampNames}
        self.covariancesModelNoB = {ampName: np.array([]) for ampName in ampNames}
        self.aMatrixNoB = {ampName: np.array([]) for ampName in ampNames}
        self.noiseMatrixNoB = {ampName: np.array([]) for ampName in ampNames}

        self.finalVars = {ampName: np.array([]) for ampName in ampNames}
        self.finalModelVars = {ampName: np.array([]) for ampName in ampNames}
        self.finalMeans = {ampName: np.array([]) for ampName in ampNames}

        # Try this as a dict of arrays.
        self.auxValues = {}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['badAmps', 'inputExpIdPairs', 'expIdMask', 'rawExpTimes',
                                        'rawMeans', 'rawVars', 'gain', 'gainErr', 'noise', 'noiseErr',
                                        'ptcFitPars', 'ptcFitParsError', 'ptcFitChiSq', 'ptcTurnoff',
                                        'aMatrixNoB', 'covariances', 'covariancesModel',
                                        'covariancesSqrtWeights', 'covariancesModelNoB',
                                        'aMatrix', 'bMatrix', 'noiseMatrix', 'noiseMatrixNoB', 'finalVars',
                                        'finalModelVars', 'finalMeans', 'photoCharges', 'histVars',
                                        'histChi2Dofs', 'kspValues', 'auxValues'])

        self.updateMetadata(setCalibInfo=True, setCalibId=True, **kwargs)
        self._validateCovarianceMatrizSizes()

    def setAmpValuesPartialDataset(
            self,
            ampName,
            inputExpIdPair=(-1, -1),
            rawExpTime=np.nan,
            rawMean=np.nan,
            rawVar=np.nan,
            photoCharge=np.nan,
            expIdMask=False,
            covariance=None,
            covSqrtWeights=None,
            gain=np.nan,
            noise=np.nan,
            histVar=np.nan,
            histChi2Dof=np.nan,
            kspValue=0.0,
            auxValues=None,
    ):
        """
        Set the amp values for a partial PTC Dataset (from cpExtractPtcTask).

        Parameters
        ----------
        ampName : `str`
            Name of the amp to set the values.
        inputExpIdPair : `tuple` [`int`]
            Exposure IDs of input pair.
        rawExpTime : `float`, optional
            Exposure time for this exposure pair.
        rawMean : `float`, optional
            Average of the means of the exposures in this pair.
        rawVar : `float`, optional
            Variance of the difference of the exposures in this pair.
        photoCharge : `float`, optional
            Integrated photocharge for flat pair for linearity calibration.
        expIdMask : `bool`, optional
            Flag setting if this exposure pair should be used (True)
            or not used (False).
        covariance : `np.ndarray` or None, optional
            Measured covariance for this exposure pair.
        covSqrtWeights : `np.ndarray` or None, optional
            Measured sqrt of covariance weights in this exposure pair.
        gain : `float`, optional
            Estimated gain for this exposure pair.
        noise : `float`, optional
            Estimated read noise for this exposure pair.
        histVar : `float`, optional
            Variance estimated from fitting a histogram with a Gaussian model.
        histChi2Dof : `float`, optional
            Chi-squared per degree of freedom from Gaussian histogram fit.
        kspValue : `float`, optional
            KS test p-value from the Gaussian histogram fit.
        """
        nanMatrix = np.full((self.covMatrixSide, self.covMatrixSide), np.nan)
        nanMatrixFit = np.full((self.covMatrixSideFullCovFit,
                               self.covMatrixSideFullCovFit), np.nan)
        if covariance is None:
            covariance = nanMatrix
        if covSqrtWeights is None:
            covSqrtWeights = nanMatrix

        self.inputExpIdPairs[ampName] = [inputExpIdPair]
        self.rawExpTimes[ampName] = np.array([rawExpTime])
        self.rawMeans[ampName] = np.array([rawMean])
        self.rawVars[ampName] = np.array([rawVar])
        self.photoCharges[ampName] = np.array([photoCharge])
        self.expIdMask[ampName] = np.array([expIdMask])
        self.covariances[ampName] = np.array([covariance])
        self.covariancesSqrtWeights[ampName] = np.array([covSqrtWeights])
        self.gain[ampName] = gain
        self.noise[ampName] = noise
        self.histVars[ampName] = np.array([histVar])
        self.histChi2Dofs[ampName] = np.array([histChi2Dof])
        self.kspValues[ampName] = np.array([kspValue])

        # From FULLCOVARIANCE model
        self.covariancesModel[ampName] = np.array([nanMatrixFit])
        self.covariancesModelNoB[ampName] = np.array([nanMatrixFit])
        self.aMatrix[ampName] = nanMatrixFit
        self.bMatrix[ampName] = nanMatrixFit
        self.aMatrixNoB[ampName] = nanMatrixFit
        self.noiseMatrix[ampName] = nanMatrixFit
        self.noiseMatrixNoB[ampName] = nanMatrixFit

    def setAuxValuesPartialDataset(self, auxDict):
        """
        Set a dictionary of auxiliary values for a partial dataset.

        Parameters
        ----------
        auxDict : `dict` [`str`, `float`]
            Dictionary of float values.
        """
        for key, value in auxDict.items():
            self.auxValues[key] = np.atleast_1d(np.array(value, dtype=np.float64))

    def updateMetadata(self, **kwargs):
        """Update calibration metadata.
        This calls the base class's method after ensuring the required
        calibration keywords will be saved.

        Parameters
        ----------
        setDate : `bool`, optional
            Update the CALIBDATE fields in the metadata to the current
            time. Defaults to False.
        kwargs :
            Other keyword parameters to set in the metadata.
        """
        super().updateMetadata(PTC_FIT_TYPE=self.ptcFitType, **kwargs)

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a calibration from a dictionary of properties.
        Must be implemented by the specific calibration subclasses.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.PhotonTransferCurveDataset`
            Constructed calibration.

        Raises
        ------
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration.
        """
        calib = cls()
        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect Photon Transfer Curve dataset  supplied. "
                               f"Expected {calib._OBSTYPE}, found {dictionary['metadata']['OBSTYPE']}")
        calib.setMetadata(dictionary['metadata'])
        calib.ptcFitType = dictionary['ptcFitType']
        calib.covMatrixSide = dictionary['covMatrixSide']
        calib.covMatrixSideFullCovFit = dictionary['covMatrixSideFullCovFit']
        calib.badAmps = np.array(dictionary['badAmps'], 'str').tolist()
        calib.ampNames = []

        # The cov matrices are square
        covMatrixSide = calib.covMatrixSide
        covMatrixSideFullCovFit = calib.covMatrixSideFullCovFit
        # Number of final signal levels
        covDimensionsProduct = len(np.array(list(dictionary['covariances'].values())[0]).ravel())
        nSignalPoints = int(covDimensionsProduct/(covMatrixSide*covMatrixSide))

        for ampName in dictionary['ampNames']:
            calib.ampNames.append(ampName)
            calib.inputExpIdPairs[ampName] = dictionary['inputExpIdPairs'][ampName]
            calib.expIdMask[ampName] = np.array(dictionary['expIdMask'][ampName])
            calib.rawExpTimes[ampName] = np.array(dictionary['rawExpTimes'][ampName], dtype=np.float64)
            calib.rawMeans[ampName] = np.array(dictionary['rawMeans'][ampName], dtype=np.float64)
            calib.rawVars[ampName] = np.array(dictionary['rawVars'][ampName], dtype=np.float64)
            calib.gain[ampName] = float(dictionary['gain'][ampName])
            calib.gainErr[ampName] = float(dictionary['gainErr'][ampName])
            calib.noise[ampName] = float(dictionary['noise'][ampName])
            calib.noiseErr[ampName] = float(dictionary['noiseErr'][ampName])
            calib.histVars[ampName] = np.array(dictionary['histVars'][ampName], dtype=np.float64)
            calib.histChi2Dofs[ampName] = np.array(dictionary['histChi2Dofs'][ampName], dtype=np.float64)
            calib.kspValues[ampName] = np.array(dictionary['kspValues'][ampName], dtype=np.float64)
            calib.ptcFitPars[ampName] = np.array(dictionary['ptcFitPars'][ampName], dtype=np.float64)
            calib.ptcFitParsError[ampName] = np.array(dictionary['ptcFitParsError'][ampName],
                                                      dtype=np.float64)
            calib.ptcFitChiSq[ampName] = float(dictionary['ptcFitChiSq'][ampName])
            calib.ptcTurnoff[ampName] = float(dictionary['ptcTurnoff'][ampName])
            if nSignalPoints > 0:
                # Regular dataset
                calib.covariances[ampName] = np.array(dictionary['covariances'][ampName],
                                                      dtype=np.float64).reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide))
                calib.covariancesModel[ampName] = np.array(
                    dictionary['covariancesModel'][ampName],
                    dtype=np.float64).reshape(
                        (nSignalPoints, covMatrixSideFullCovFit, covMatrixSideFullCovFit))
                calib.covariancesSqrtWeights[ampName] = np.array(
                    dictionary['covariancesSqrtWeights'][ampName],
                    dtype=np.float64).reshape(
                        (nSignalPoints, covMatrixSide, covMatrixSide))
                calib.aMatrix[ampName] = np.array(dictionary['aMatrix'][ampName],
                                                  dtype=np.float64).reshape(
                    (covMatrixSideFullCovFit, covMatrixSideFullCovFit))
                calib.bMatrix[ampName] = np.array(dictionary['bMatrix'][ampName],
                                                  dtype=np.float64).reshape(
                    (covMatrixSideFullCovFit, covMatrixSideFullCovFit))
                calib.covariancesModelNoB[ampName] = np.array(
                    dictionary['covariancesModelNoB'][ampName], dtype=np.float64).reshape(
                        (nSignalPoints, covMatrixSideFullCovFit, covMatrixSideFullCovFit))
                calib.aMatrixNoB[ampName] = np.array(
                    dictionary['aMatrixNoB'][ampName],
                    dtype=np.float64).reshape((covMatrixSideFullCovFit, covMatrixSideFullCovFit))
                calib.noiseMatrix[ampName] = np.array(
                    dictionary['noiseMatrix'][ampName],
                    dtype=np.float64).reshape((covMatrixSideFullCovFit, covMatrixSideFullCovFit))
                calib.noiseMatrixNoB[ampName] = np.array(
                    dictionary['noiseMatrixNoB'][ampName],
                    dtype=np.float64).reshape((covMatrixSideFullCovFit, covMatrixSideFullCovFit))
            else:
                # Empty dataset
                calib.covariances[ampName] = np.array([], dtype=np.float64)
                calib.covariancesModel[ampName] = np.array([], dtype=np.float64)
                calib.covariancesSqrtWeights[ampName] = np.array([], dtype=np.float64)
                calib.aMatrix[ampName] = np.array([], dtype=np.float64)
                calib.bMatrix[ampName] = np.array([], dtype=np.float64)
                calib.covariancesModelNoB[ampName] = np.array([], dtype=np.float64)
                calib.aMatrixNoB[ampName] = np.array([], dtype=np.float64)
                calib.noiseMatrix[ampName] = np.array([], dtype=np.float64)
                calib.noiseMatrixNoB[ampName] = np.array([], dtype=np.float64)

            calib.finalVars[ampName] = np.array(dictionary['finalVars'][ampName], dtype=np.float64)
            calib.finalModelVars[ampName] = np.array(dictionary['finalModelVars'][ampName], dtype=np.float64)
            calib.finalMeans[ampName] = np.array(dictionary['finalMeans'][ampName], dtype=np.float64)
            calib.photoCharges[ampName] = np.array(dictionary['photoCharges'][ampName], dtype=np.float64)

        for key, value in dictionary['auxValues'].items():
            calib.auxValues[key] = np.atleast_1d(np.array(value, dtype=np.float64))

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.
        The dictionary should be able to be round-tripped through
        `fromDict`.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()

        outDict = dict()
        metadata = self.getMetadata()
        outDict['metadata'] = metadata

        def _dictOfArraysToDictOfLists(dictOfArrays):
            dictOfLists = {}
            for key, value in dictOfArrays.items():
                dictOfLists[key] = value.ravel().tolist()

            return dictOfLists

        outDict['ptcFitType'] = self.ptcFitType
        outDict['covMatrixSide'] = self.covMatrixSide
        outDict['covMatrixSideFullCovFit'] = self.covMatrixSideFullCovFit
        outDict['ampNames'] = self.ampNames
        outDict['badAmps'] = self.badAmps
        outDict['inputExpIdPairs'] = self.inputExpIdPairs
        outDict['expIdMask'] = _dictOfArraysToDictOfLists(self.expIdMask)
        outDict['rawExpTimes'] = _dictOfArraysToDictOfLists(self.rawExpTimes)
        outDict['rawMeans'] = _dictOfArraysToDictOfLists(self.rawMeans)
        outDict['rawVars'] = _dictOfArraysToDictOfLists(self.rawVars)
        outDict['gain'] = self.gain
        outDict['gainErr'] = self.gainErr
        outDict['noise'] = self.noise
        outDict['noiseErr'] = self.noiseErr
        outDict['histVars'] = self.histVars
        outDict['histChi2Dofs'] = self.histChi2Dofs
        outDict['kspValues'] = self.kspValues
        outDict['ptcFitPars'] = _dictOfArraysToDictOfLists(self.ptcFitPars)
        outDict['ptcFitParsError'] = _dictOfArraysToDictOfLists(self.ptcFitParsError)
        outDict['ptcFitChiSq'] = self.ptcFitChiSq
        outDict['ptcTurnoff'] = self.ptcTurnoff
        outDict['covariances'] = _dictOfArraysToDictOfLists(self.covariances)
        outDict['covariancesModel'] = _dictOfArraysToDictOfLists(self.covariancesModel)
        outDict['covariancesSqrtWeights'] = _dictOfArraysToDictOfLists(self.covariancesSqrtWeights)
        outDict['aMatrix'] = _dictOfArraysToDictOfLists(self.aMatrix)
        outDict['bMatrix'] = _dictOfArraysToDictOfLists(self.bMatrix)
        outDict['noiseMatrix'] = _dictOfArraysToDictOfLists(self.noiseMatrix)
        outDict['covariancesModelNoB'] = _dictOfArraysToDictOfLists(self.covariancesModelNoB)
        outDict['aMatrixNoB'] = _dictOfArraysToDictOfLists(self.aMatrixNoB)
        outDict['noiseMatrixNoB'] = _dictOfArraysToDictOfLists(self.noiseMatrixNoB)
        outDict['finalVars'] = _dictOfArraysToDictOfLists(self.finalVars)
        outDict['finalModelVars'] = _dictOfArraysToDictOfLists(self.finalModelVars)
        outDict['finalMeans'] = _dictOfArraysToDictOfLists(self.finalMeans)
        outDict['photoCharges'] = _dictOfArraysToDictOfLists(self.photoCharges)
        outDict['auxValues'] = _dictOfArraysToDictOfLists(self.auxValues)

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.
        This method uses the `fromDict` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables to use to construct the datasetPtc.

        Returns
        -------
        calib : `lsst.ip.isr.PhotonTransferCurveDataset`
            The calibration defined in the tables.
        """
        ptcTable = tableList[0]

        metadata = ptcTable.meta
        inDict = dict()
        inDict['metadata'] = metadata
        inDict['ampNames'] = []
        inDict['ptcFitType'] = []
        inDict['covMatrixSide'] = []
        inDict['covMatrixSideFullCovFit'] = []
        inDict['inputExpIdPairs'] = dict()
        inDict['expIdMask'] = dict()
        inDict['rawExpTimes'] = dict()
        inDict['rawMeans'] = dict()
        inDict['rawVars'] = dict()
        inDict['gain'] = dict()
        inDict['gainErr'] = dict()
        inDict['noise'] = dict()
        inDict['noiseErr'] = dict()
        inDict['histVars'] = dict()
        inDict['histChi2Dofs'] = dict()
        inDict['kspValues'] = dict()
        inDict['ptcFitPars'] = dict()
        inDict['ptcFitParsError'] = dict()
        inDict['ptcFitChiSq'] = dict()
        inDict['ptcTurnoff'] = dict()
        inDict['covariances'] = dict()
        inDict['covariancesModel'] = dict()
        inDict['covariancesSqrtWeights'] = dict()
        inDict['aMatrix'] = dict()
        inDict['bMatrix'] = dict()
        inDict['noiseMatrix'] = dict()
        inDict['covariancesModelNoB'] = dict()
        inDict['aMatrixNoB'] = dict()
        inDict['noiseMatrixNoB'] = dict()
        inDict['finalVars'] = dict()
        inDict['finalModelVars'] = dict()
        inDict['finalMeans'] = dict()
        inDict['badAmps'] = []
        inDict['photoCharges'] = dict()

        calibVersion = metadata['PTC_VERSION']
        if calibVersion == 1.0:
            cls().log.warning(f"Previous version found for PTC dataset: {calibVersion}. "
                              f"Setting 'ptcTurnoff' in all amps to last value in 'finalMeans'.")
        for record in ptcTable:
            ampName = record['AMPLIFIER_NAME']

            inDict['ptcFitType'] = record['PTC_FIT_TYPE']
            inDict['covMatrixSide'] = record['COV_MATRIX_SIDE']
            inDict['ampNames'].append(ampName)
            inDict['inputExpIdPairs'][ampName] = record['INPUT_EXP_ID_PAIRS'].tolist()
            inDict['expIdMask'][ampName] = record['EXP_ID_MASK']
            inDict['rawExpTimes'][ampName] = record['RAW_EXP_TIMES']
            inDict['rawMeans'][ampName] = record['RAW_MEANS']
            inDict['rawVars'][ampName] = record['RAW_VARS']
            inDict['gain'][ampName] = record['GAIN']
            inDict['gainErr'][ampName] = record['GAIN_ERR']
            inDict['noise'][ampName] = record['NOISE']
            inDict['noiseErr'][ampName] = record['NOISE_ERR']
            inDict['ptcFitPars'][ampName] = record['PTC_FIT_PARS']
            inDict['ptcFitParsError'][ampName] = record['PTC_FIT_PARS_ERROR']
            inDict['ptcFitChiSq'][ampName] = record['PTC_FIT_CHI_SQ']
            inDict['covariances'][ampName] = record['COVARIANCES']
            inDict['covariancesModel'][ampName] = record['COVARIANCES_MODEL']
            inDict['covariancesSqrtWeights'][ampName] = record['COVARIANCES_SQRT_WEIGHTS']
            inDict['aMatrix'][ampName] = record['A_MATRIX']
            inDict['bMatrix'][ampName] = record['B_MATRIX']
            inDict['covariancesModelNoB'][ampName] = record['COVARIANCES_MODEL_NO_B']
            inDict['aMatrixNoB'][ampName] = record['A_MATRIX_NO_B']
            inDict['finalVars'][ampName] = record['FINAL_VARS']
            inDict['finalModelVars'][ampName] = record['FINAL_MODEL_VARS']
            inDict['finalMeans'][ampName] = record['FINAL_MEANS']
            inDict['badAmps'] = record['BAD_AMPS'].tolist()
            inDict['photoCharges'][ampName] = record['PHOTO_CHARGE']
            if calibVersion == 1.0:
                mask = record['FINAL_MEANS'].mask
                array = record['FINAL_MEANS'][~mask]
                if len(array) > 0:
                    inDict['ptcTurnoff'][ampName] = record['FINAL_MEANS'][~mask][-1]
                else:
                    inDict['ptcTurnoff'][ampName] = np.nan
            else:
                inDict['ptcTurnoff'][ampName] = record['PTC_TURNOFF']
            if calibVersion < 1.2:
                inDict['histVars'][ampName] = np.array([np.nan])
                inDict['histChi2Dofs'][ampName] = np.array([np.nan])
                inDict['kspValues'][ampName] = np.array([0.0])
            else:
                inDict['histVars'][ampName] = record['HIST_VARS']
                inDict['histChi2Dofs'][ampName] = record['HIST_CHI2_DOFS']
                inDict['kspValues'][ampName] = record['KS_PVALUES']
            if calibVersion < 1.3:
                nanMatrix = np.full_like(inDict['aMatrix'][ampName], np.nan)
                inDict['noiseMatrix'][ampName] = nanMatrix
                inDict['noiseMatrixNoB'][ampName] = nanMatrix
            else:
                inDict['noiseMatrix'][ampName] = record['NOISE_MATRIX']
                inDict['noiseMatrixNoB'][ampName] = record['NOISE_MATRIX_NO_B']
            if calibVersion < 1.5:
                # Matched to `COV_MATRIX_SIDE`. Same for all amps.
                inDict['covMatrixSideFullCovFit'] = inDict['covMatrixSide']
            else:
                inDict['covMatrixSideFullCovFit'] = record['COV_MATRIX_SIDE_FULL_COV_FIT']

        inDict['auxValues'] = {}
        record = ptcTable[0]
        for col in record.columns.keys():
            if col.startswith('PTCAUX_'):
                parts = col.split('PTCAUX_')
                inDict['auxValues'][parts[1]] = record[col]

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the linearity calibration
            information.
        """
        tableList = []
        self.updateMetadata()

        badAmps = np.array(self.badAmps) if len(self.badAmps) else np.array([], dtype="U3")

        catalogList = []
        for ampName in self.ampNames:
            ampDict = {
                'AMPLIFIER_NAME': ampName,
                'PTC_FIT_TYPE': self.ptcFitType,
                'COV_MATRIX_SIDE': self.covMatrixSide,
                'COV_MATRIX_SIDE_FULL_COV_FIT': self.covMatrixSideFullCovFit,
                'INPUT_EXP_ID_PAIRS': self.inputExpIdPairs[ampName],
                'EXP_ID_MASK': self.expIdMask[ampName],
                'RAW_EXP_TIMES': self.rawExpTimes[ampName],
                'RAW_MEANS': self.rawMeans[ampName],
                'RAW_VARS': self.rawVars[ampName],
                'GAIN': self.gain[ampName],
                'GAIN_ERR': self.gainErr[ampName],
                'NOISE': self.noise[ampName],
                'NOISE_ERR': self.noiseErr[ampName],
                'HIST_VARS': self.histVars[ampName],
                'HIST_CHI2_DOFS': self.histChi2Dofs[ampName],
                'KS_PVALUES': self.kspValues[ampName],
                'PTC_FIT_PARS': np.array(self.ptcFitPars[ampName]),
                'PTC_FIT_PARS_ERROR': np.array(self.ptcFitParsError[ampName]),
                'PTC_FIT_CHI_SQ': self.ptcFitChiSq[ampName],
                'PTC_TURNOFF': self.ptcTurnoff[ampName],
                'A_MATRIX': self.aMatrix[ampName].ravel(),
                'B_MATRIX': self.bMatrix[ampName].ravel(),
                'A_MATRIX_NO_B': self.aMatrixNoB[ampName].ravel(),
                'NOISE_MATRIX': self.noiseMatrix[ampName].ravel(),
                'NOISE_MATRIX_NO_B': self.noiseMatrixNoB[ampName].ravel(),
                'BAD_AMPS': badAmps,
                'PHOTO_CHARGE': self.photoCharges[ampName],
                'COVARIANCES': self.covariances[ampName].ravel(),
                'COVARIANCES_MODEL': self.covariancesModel[ampName].ravel(),
                'COVARIANCES_SQRT_WEIGHTS': self.covariancesSqrtWeights[ampName].ravel(),
                'COVARIANCES_MODEL_NO_B': self.covariancesModelNoB[ampName].ravel(),
                'FINAL_VARS': self.finalVars[ampName],
                'FINAL_MODEL_VARS': self.finalModelVars[ampName],
                'FINAL_MEANS': self.finalMeans[ampName],
            }

            if self.auxValues:
                for key, value in self.auxValues.items():
                    ampDict[f"PTCAUX_{key}"] = value

            catalogList.append(ampDict)

        catalog = Table(catalogList)

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta
        tableList.append(catalog)

        return tableList

    def fromDetector(self, detector):
        """Read metadata parameters from a detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.detector`
            Input detector with parameters to use.

        Returns
        -------
        calib : `lsst.ip.isr.PhotonTransferCurveDataset`
            The calibration constructed from the detector.
        """

        pass

    def getExpIdsUsed(self, ampName):
        """Get the exposures used, i.e. not discarded, for a given amp.
        If no mask has been created yet, all exposures are returned.

        Parameters
        ----------
        ampName : `str`

        Returns
        -------
        expIdsUsed : `list` [`tuple`]
            List of pairs of exposure ids used in PTC.
        """
        if len(self.expIdMask[ampName]) == 0:
            return self.inputExpIdPairs[ampName]

        # if the mask exists it had better be the same length as the expIdPairs
        assert len(self.expIdMask[ampName]) == len(self.inputExpIdPairs[ampName])

        pairs = self.inputExpIdPairs[ampName]
        mask = self.expIdMask[ampName]
        # cast to bool required because numpy
        try:
            expIdsUsed = [(exp1, exp2) for ((exp1, exp2), m) in zip(pairs, mask) if m]
        except ValueError:
            self.log.warning("The PTC file was written incorrectly; you should rerun the "
                             "PTC solve task if possible.")
            expIdsUsed = []
            for pairList, m in zip(pairs, mask):
                if m:
                    expIdsUsed.append(pairList[0])

        return expIdsUsed

    def getGoodAmps(self):
        """Get the good amps from this PTC."""
        return [amp for amp in self.ampNames if amp not in self.badAmps]

    def getGoodPoints(self, ampName):
        """Get the good points used for a given amp in the PTC.

        Parameters
        ----------
        ampName : `str`
            Amplifier's name.

        Returns
        -------
        goodPoints : `np.ndarray`
            Boolean array of good points used in PTC.
        """
        return self.expIdMask[ampName]

    def validateGainNoiseTurnoffValues(self, ampName, doWarn=False):
        """Ensure the gain, read noise, and PTC turnoff have
        sensible values.

        Parameters
        ----------
        ampName : `str`
            Amplifier's name.
        """

        gain = self.gain[ampName]
        noise = self.noise[ampName]
        ptcTurnoff = self.ptcTurnoff[ampName]

        # Check if gain is not positive or is np.nan
        if not (isinstance(gain, (int, float)) and gain > 0) or math.isnan(gain):
            if doWarn:
                self.log.warning(f"Invalid gain value {gain}"
                                 " Setting to default: Gain=1")
            gain = 1

        # Check if noise is not positive or is np.nan
        if not (isinstance(noise, (int, float)) and noise > 0) or math.isnan(noise):
            if doWarn:
                self.log.warning(f"Invalid noise value: {noise}"
                                 " Setting to default: Noise=1")
            noise = 1

        # Check if ptcTurnoff is not positive or is np.nan
        if not (isinstance(ptcTurnoff, (int, float)) and ptcTurnoff > 0) or math.isnan(ptcTurnoff):
            if doWarn:
                self.log.warning(f"Invalid PTC turnoff value: {ptcTurnoff}"
                                 " Setting to default: PTC Turnoff=2e19")
            ptcTurnoff = 2e19

        self.gain[ampName] = gain
        self.noise[ampName] = noise
        self.ptcTurnoff[ampName] = ptcTurnoff

    def _validateCovarianceMatrizSizes(self):
        """Ensure  covMatrixSideFullCovFit <= covMatrixSide."""
        if self.covMatrixSideFullCovFit > self.covMatrixSide:
            self.log.warning("covMatrixSideFullCovFit > covMatrixSide "
                             f"({self.covMatrixSideFullCovFit} > {self.covMatrixSide})."
                             "Setting the former to the latter.")
            self.covMatrixSideFullCovFit = self.covMatrixSide
