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
import numpy as np
from astropy.table import Table

from lsst.ip.isr import IsrCalib

__all__ = ['PhotonTransferCurveDataset']


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

    ptcFitType : `str`
        Type of model fitted to the PTC: "POLYNOMIAL", "EXPAPPROXIMATION",
        or "FULLCOVARIANCE".

    covMatrixSide : `int`
        Maximum lag of covariances (size of square covariance matrices).

    kwargs : `dict`, optional
        Other keyword arguments to pass to the parent init.

    Notes
    -----
    The stored attributes are:
    badAmps : `list`
        List with bad amplifiers names.
    inputExpIdPairs : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the input exposures IDs.
    expIdMask : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the mask produced after
        outlier rejection. The mask produced by the "FULLCOVARIANCE"
        option may differ from the one produced in the other two PTC
        fit types.
    rawExpTimes : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the unmasked exposure times.
    rawMeans : `dict`, [`str`, `list`]
        Dictionary keyed by amp namescontaining the unmasked average of the
        means of the exposures in each flat pair.
    rawVars : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the variance of the
        difference image of the exposures in each flat pair.
    gain : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the fitted gains.
    gainErr : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the errors on the
        fitted gains.
    noise : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the fitted noise.
    noiseErr : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the errors on the fitted
        noise.
    ptcFitPars : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the fitted parameters of the
        PTC model for ptcFitTye in ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcFitParsError : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the errors on the fitted
        parameters of the PTC model for ptcFitTye in
        ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcFitChiSq : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the reduced chi squared
        of the fit for ptcFitTye in ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcTurnoff : `float`
        Flux value (in ADU) where the variance of the PTC curve starts
        decreasing consistently.
    covariances : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing a list of measured
        covariances per mean flux.
    covariancesModel : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containinging covariances model
        (Eq. 20 of Astier+19) per mean flux.
    covariancesSqrtWeights : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containinging sqrt. of covariances
        weights.
    aMatrix : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the "a" parameters from
        the model in Eq. 20 of Astier+19.
    bMatrix : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the "b" parameters from
        the model in Eq. 20 of Astier+19.
    covariancesModelNoB : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing covariances model
        (with 'b'=0 in Eq. 20 of Astier+19)
        per mean flux.
    aMatrixNoB : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the "a" parameters from the
        model in Eq. 20 of Astier+19
        (and 'b' = 0).
    finalVars : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the masked variance of the
        difference image of each flat
        pair. If needed, each array will be right-padded with
        np.nan to match the length of rawExpTimes.
    finalModelVars : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the masked modeled
        variance of the difference image of each flat pair. If needed, each
        array will be right-padded with np.nan to match the length of
        rawExpTimes.
    finalMeans : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing the masked average of the
        means of the exposures in each flat pair. If needed, each array
        will be right-padded with np.nan to match the length of
        rawExpTimes.
    photoCharge : `dict`, [`str`, `list`]
         Dictionary keyed by amp names containing the integrated photocharge
         for linearity calibration.

    Returns
    -------
    `lsst.cp.pipe.ptc.PhotonTransferCurveDataset`
        Output dataset from MeasurePhotonTransferCurveTask.

    Notes
    -----
    Version 1.1 adds the `ptcTurnoff` attribute.
    """

    _OBSTYPE = 'PTC'
    _SCHEMA = 'Gen3 Photon Transfer Curve'
    _VERSION = 1.1

    def __init__(self, ampNames=[], ptcFitType=None, covMatrixSide=1, **kwargs):

        self.ptcFitType = ptcFitType
        self.ampNames = ampNames
        self.covMatrixSide = covMatrixSide

        self.badAmps = [np.nan]

        self.inputExpIdPairs = {ampName: [] for ampName in ampNames}
        self.expIdMask = {ampName: [] for ampName in ampNames}
        self.rawExpTimes = {ampName: [] for ampName in ampNames}
        self.rawMeans = {ampName: [] for ampName in ampNames}
        self.rawVars = {ampName: [] for ampName in ampNames}
        self.photoCharge = {ampName: [] for ampName in ampNames}

        self.gain = {ampName: np.nan for ampName in ampNames}
        self.gainErr = {ampName: np.nan for ampName in ampNames}
        self.noise = {ampName: np.nan for ampName in ampNames}
        self.noiseErr = {ampName: np.nan for ampName in ampNames}

        self.ptcFitPars = {ampName: [] for ampName in ampNames}
        self.ptcFitParsError = {ampName: [] for ampName in ampNames}
        self.ptcFitChiSq = {ampName: np.nan for ampName in ampNames}
        self.ptcTurnoff = {ampName: np.nan for ampName in ampNames}

        self.covariances = {ampName: [] for ampName in ampNames}
        self.covariancesModel = {ampName: [] for ampName in ampNames}
        self.covariancesSqrtWeights = {ampName: [] for ampName in ampNames}
        self.aMatrix = {ampName: np.nan for ampName in ampNames}
        self.bMatrix = {ampName: np.nan for ampName in ampNames}
        self.covariancesModelNoB = {ampName: [] for ampName in ampNames}
        self.aMatrixNoB = {ampName: np.nan for ampName in ampNames}

        self.finalVars = {ampName: [] for ampName in ampNames}
        self.finalModelVars = {ampName: [] for ampName in ampNames}
        self.finalMeans = {ampName: [] for ampName in ampNames}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['badAmps', 'inputExpIdPairs', 'expIdMask', 'rawExpTimes',
                                        'rawMeans', 'rawVars', 'gain', 'gainErr', 'noise', 'noiseErr',
                                        'ptcFitPars', 'ptcFitParsError', 'ptcFitChiSq', 'ptcTurnoff',
                                        'aMatrixNoB', 'covariances', 'covariancesModel',
                                        'covariancesSqrtWeights', 'covariancesModelNoB',
                                        'aMatrix', 'bMatrix', 'finalVars', 'finalModelVars', 'finalMeans',
                                        'photoCharge'])

    def setAmpValues(self, ampName, inputExpIdPair=[(np.nan, np.nan)], expIdMask=[np.nan],
                     rawExpTime=[np.nan], rawMean=[np.nan], rawVar=[np.nan], photoCharge=[np.nan],
                     gain=np.nan, gainErr=np.nan, noise=np.nan, noiseErr=np.nan, ptcFitPars=[np.nan],
                     ptcFitParsError=[np.nan], ptcFitChiSq=np.nan, ptcTurnoff=np.nan, covArray=[],
                     covArrayModel=[], covSqrtWeights=[], aMatrix=[], bMatrix=[], covArrayModelNoB=[],
                     aMatrixNoB=[], finalVar=[np.nan], finalModelVar=[np.nan], finalMean=[np.nan]):
        """Function to initialize an amp of a PhotonTransferCurveDataset.

        Notes
        -----
        The parameters are all documented in `init`.
        """
        nanMatrix = np.full((self.covMatrixSide, self.covMatrixSide), np.nan)
        if len(covArray) == 0:
            covArray = [nanMatrix]
        if len(covArrayModel) == 0:
            covArrayModel = [nanMatrix]
        if len(covSqrtWeights) == 0:
            covSqrtWeights = [nanMatrix]
        if len(covArrayModelNoB) == 0:
            covArrayModelNoB = [nanMatrix]
        if len(aMatrix) == 0:
            aMatrix = nanMatrix
        if len(bMatrix) == 0:
            bMatrix = nanMatrix
        if len(aMatrixNoB) == 0:
            aMatrixNoB = nanMatrix

        self.inputExpIdPairs[ampName] = inputExpIdPair
        self.expIdMask[ampName] = expIdMask
        self.rawExpTimes[ampName] = rawExpTime
        self.rawMeans[ampName] = rawMean
        self.rawVars[ampName] = rawVar
        self.photoCharge[ampName] = photoCharge
        self.gain[ampName] = gain
        self.gainErr[ampName] = gainErr
        self.noise[ampName] = noise
        self.noiseErr[ampName] = noiseErr
        self.ptcFitPars[ampName] = ptcFitPars
        self.ptcFitParsError[ampName] = ptcFitParsError
        self.ptcFitChiSq[ampName] = ptcFitChiSq
        self.ptcTurnoff[ampName] = ptcTurnoff
        self.covariances[ampName] = covArray
        self.covariancesSqrtWeights[ampName] = covSqrtWeights
        self.covariancesModel[ampName] = covArrayModel
        self.covariancesModelNoB[ampName] = covArrayModelNoB
        self.aMatrix[ampName] = aMatrix
        self.bMatrix[ampName] = bMatrix
        self.aMatrixNoB[ampName] = aMatrixNoB
        self.ptcFitPars[ampName] = ptcFitPars
        self.ptcFitParsError[ampName] = ptcFitParsError
        self.ptcFitChiSq[ampName] = ptcFitChiSq
        self.finalVars[ampName] = finalVar
        self.finalModelVars[ampName] = finalModelVar
        self.finalMeans[ampName] = finalMean

    def updateMetadata(self, setDate=False, **kwargs):
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
        kwargs['PTC_FIT_TYPE'] = self.ptcFitType

        super().updateMetadata(setDate=setDate, **kwargs)

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
        calib : `lsst.ip.isr.CalibType`
            Constructed calibration.
        Raises
        ------
        RuntimeError :
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
        calib.badAmps = np.array(dictionary['badAmps'], 'str').tolist()
        # The cov matrices are square
        covMatrixSide = calib.covMatrixSide
        # Number of final signal levels
        covDimensionsProduct = len(np.array(list(dictionary['covariances'].values())[0]).ravel())
        nSignalPoints = int(covDimensionsProduct/(covMatrixSide*covMatrixSide))

        for ampName in dictionary['ampNames']:
            covsAmp = np.array(dictionary['covariances'][ampName]).reshape((nSignalPoints, covMatrixSide,
                                                                            covMatrixSide))

            # After cpPtcExtract runs in the PTC pipeline, the datasets
            # created ('PARTIAL' and 'DUMMY') have a single measurement.
            # Apply the maskign to the final ptcDataset, after running
            # cpPtcSolve.
            if len(covsAmp) > 1:
                # Masks for covariances padding in `toTable`
                maskCovsAmp = np.array([~np.isnan(entry).all() for entry in covsAmp])
                maskAmp = ~np.isnan(np.array(dictionary['finalMeans'][ampName]))
            else:
                maskCovsAmp = np.array([True])
                maskAmp = np.array([True])

            calib.ampNames.append(ampName)
            calib.inputExpIdPairs[ampName] = np.array(dictionary['inputExpIdPairs'][ampName]).tolist()
            calib.expIdMask[ampName] = np.array(dictionary['expIdMask'][ampName]).tolist()
            calib.rawExpTimes[ampName] = np.array(dictionary['rawExpTimes'][ampName]).tolist()
            calib.rawMeans[ampName] = np.array(dictionary['rawMeans'][ampName]).tolist()
            calib.rawVars[ampName] = np.array(dictionary['rawVars'][ampName]).tolist()
            calib.gain[ampName] = np.array(dictionary['gain'][ampName]).tolist()
            calib.gainErr[ampName] = np.array(dictionary['gainErr'][ampName]).tolist()
            calib.noise[ampName] = np.array(dictionary['noise'][ampName]).tolist()
            calib.noiseErr[ampName] = np.array(dictionary['noiseErr'][ampName]).tolist()
            calib.ptcFitPars[ampName] = np.array(dictionary['ptcFitPars'][ampName]).tolist()
            calib.ptcFitParsError[ampName] = np.array(dictionary['ptcFitParsError'][ampName]).tolist()
            calib.ptcFitChiSq[ampName] = np.array(dictionary['ptcFitChiSq'][ampName]).tolist()
            calib.ptcTurnoff[ampName] = np.array(dictionary['ptcTurnoff'][ampName]).tolist()
            calib.covariances[ampName] = covsAmp[maskCovsAmp].tolist()
            calib.covariancesModel[ampName] = np.array(
                dictionary['covariancesModel'][ampName]).reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide))[maskCovsAmp].tolist()
            calib.covariancesSqrtWeights[ampName] = np.array(
                dictionary['covariancesSqrtWeights'][ampName]).reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide))[maskCovsAmp].tolist()
            calib.aMatrix[ampName] = np.array(dictionary['aMatrix'][ampName]).reshape(
                (covMatrixSide, covMatrixSide)).tolist()
            calib.bMatrix[ampName] = np.array(dictionary['bMatrix'][ampName]).reshape(
                (covMatrixSide, covMatrixSide)).tolist()
            calib.covariancesModelNoB[ampName] = np.array(
                dictionary['covariancesModelNoB'][ampName]).reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide))[maskCovsAmp].tolist()
            calib.aMatrixNoB[ampName] = np.array(
                dictionary['aMatrixNoB'][ampName]).reshape((covMatrixSide, covMatrixSide)).tolist()
            calib.finalVars[ampName] = np.array(dictionary['finalVars'][ampName])[maskAmp].tolist()
            calib.finalModelVars[ampName] = np.array(dictionary['finalModelVars'][ampName])[maskAmp].tolist()
            calib.finalMeans[ampName] = np.array(dictionary['finalMeans'][ampName])[maskAmp].tolist()
            calib.photoCharge[ampName] = np.array(dictionary['photoCharge'][ampName]).tolist()
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

        outDict['ptcFitType'] = self.ptcFitType
        outDict['covMatrixSide'] = self.covMatrixSide
        outDict['ampNames'] = self.ampNames
        outDict['badAmps'] = self.badAmps
        outDict['inputExpIdPairs'] = self.inputExpIdPairs
        outDict['expIdMask'] = self.expIdMask
        outDict['rawExpTimes'] = self.rawExpTimes
        outDict['rawMeans'] = self.rawMeans
        outDict['rawVars'] = self.rawVars
        outDict['gain'] = self.gain
        outDict['gainErr'] = self.gainErr
        outDict['noise'] = self.noise
        outDict['noiseErr'] = self.noiseErr
        outDict['ptcFitPars'] = self.ptcFitPars
        outDict['ptcFitParsError'] = self.ptcFitParsError
        outDict['ptcFitChiSq'] = self.ptcFitChiSq
        outDict['ptcTurnoff'] = self.ptcTurnoff
        outDict['covariances'] = self.covariances
        outDict['covariancesModel'] = self.covariancesModel
        outDict['covariancesSqrtWeights'] = self.covariancesSqrtWeights
        outDict['aMatrix'] = self.aMatrix
        outDict['bMatrix'] = self.bMatrix
        outDict['covariancesModelNoB'] = self.covariancesModelNoB
        outDict['aMatrixNoB'] = self.aMatrixNoB
        outDict['finalVars'] = self.finalVars
        outDict['finalModelVars'] = self.finalModelVars
        outDict['finalMeans'] = self.finalMeans
        outDict['photoCharge'] = self.photoCharge

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
        calib : `lsst.cp.pipe.`
            The calibration defined in the tables.
        """
        ptcTable = tableList[0]

        metadata = ptcTable.meta
        inDict = dict()
        inDict['metadata'] = metadata
        inDict['ampNames'] = []
        inDict['ptcFitType'] = []
        inDict['covMatrixSide'] = []
        inDict['inputExpIdPairs'] = dict()
        inDict['expIdMask'] = dict()
        inDict['rawExpTimes'] = dict()
        inDict['rawMeans'] = dict()
        inDict['rawVars'] = dict()
        inDict['gain'] = dict()
        inDict['gainErr'] = dict()
        inDict['noise'] = dict()
        inDict['noiseErr'] = dict()
        inDict['ptcFitPars'] = dict()
        inDict['ptcFitParsError'] = dict()
        inDict['ptcFitChiSq'] = dict()
        inDict['ptcTurnoff'] = dict()
        inDict['covariances'] = dict()
        inDict['covariancesModel'] = dict()
        inDict['covariancesSqrtWeights'] = dict()
        inDict['aMatrix'] = dict()
        inDict['bMatrix'] = dict()
        inDict['covariancesModelNoB'] = dict()
        inDict['aMatrixNoB'] = dict()
        inDict['finalVars'] = dict()
        inDict['finalModelVars'] = dict()
        inDict['finalMeans'] = dict()
        inDict['badAmps'] = []
        inDict['photoCharge'] = dict()

        calibVersion = metadata['PTC_VERSION']
        if calibVersion == 1.0:
            cls().log.warning(f"Previous version found for PTC dataset: {calibVersion}. "
                              f"Setting 'ptcTurnoff' in all amps to last value in 'finalMeans'.")
        for record in ptcTable:
            ampName = record['AMPLIFIER_NAME']

            inDict['ptcFitType'] = record['PTC_FIT_TYPE']
            inDict['covMatrixSide'] = record['COV_MATRIX_SIDE']
            inDict['ampNames'].append(ampName)
            inDict['inputExpIdPairs'][ampName] = record['INPUT_EXP_ID_PAIRS']
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
            inDict['badAmps'] = record['BAD_AMPS']
            inDict['photoCharge'][ampName] = record['PHOTO_CHARGE']
            if calibVersion == 1.0:
                mask = record['FINAL_MEANS'].mask
                array = record['FINAL_MEANS'][~mask]
                if len(array) > 0:
                    inDict['ptcTurnoff'][ampName] = record['FINAL_MEANS'][~mask][-1]
                else:
                    inDict['ptcTurnoff'][ampName] = np.nan
            else:
                inDict['ptcTurnoff'][ampName] = record['PTC_TURNOFF']
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
        nPoints = []
        for i, ampName in enumerate(self.ampNames):
            nPoints.append(len(list(self.covariances.values())[i]))
        nSignalPoints = max(nPoints)
        nPadPoints = {}
        for i, ampName in enumerate(self.ampNames):
            nPadPoints[ampName] = nSignalPoints - len(list(self.covariances.values())[i])
        covMatrixSide = self.covMatrixSide

        catalog = Table([{'AMPLIFIER_NAME': ampName,
                          'PTC_FIT_TYPE': self.ptcFitType,
                          'COV_MATRIX_SIDE': self.covMatrixSide,
                          'INPUT_EXP_ID_PAIRS': self.inputExpIdPairs[ampName]
                         if len(self.expIdMask[ampName]) else np.nan,
                          'EXP_ID_MASK': self.expIdMask[ampName]
                         if len(self.expIdMask[ampName]) else np.nan,
                          'RAW_EXP_TIMES': np.array(self.rawExpTimes[ampName]).tolist()
                         if len(self.rawExpTimes[ampName]) else np.nan,
                          'RAW_MEANS': np.array(self.rawMeans[ampName]).tolist()
                         if len(self.rawMeans[ampName]) else np.nan,
                          'RAW_VARS': np.array(self.rawVars[ampName]).tolist()
                         if len(self.rawVars[ampName]) else np.nan,
                          'GAIN': self.gain[ampName],
                          'GAIN_ERR': self.gainErr[ampName],
                          'NOISE': self.noise[ampName],
                          'NOISE_ERR': self.noiseErr[ampName],
                          'PTC_FIT_PARS': np.array(self.ptcFitPars[ampName]).tolist(),
                          'PTC_FIT_PARS_ERROR': np.array(self.ptcFitParsError[ampName]).tolist(),
                          'PTC_FIT_CHI_SQ': self.ptcFitChiSq[ampName],
                          'PTC_TURNOFF': self.ptcTurnoff[ampName],
                          'COVARIANCES': np.pad(np.array(self.covariances[ampName]),
                                                ((0, nPadPoints[ampName]), (0, 0), (0, 0)),
                                                'constant', constant_values=np.nan).reshape(
                              nSignalPoints*covMatrixSide**2).tolist(),
                          'COVARIANCES_MODEL': np.pad(np.array(self.covariancesModel[ampName]),
                                                      ((0, nPadPoints[ampName]), (0, 0), (0, 0)),
                                                      'constant', constant_values=np.nan).reshape(
                              nSignalPoints*covMatrixSide**2).tolist(),
                          'COVARIANCES_SQRT_WEIGHTS': np.pad(np.array(self.covariancesSqrtWeights[ampName]),
                                                             ((0, nPadPoints[ampName]), (0, 0), (0, 0)),
                                                             'constant', constant_values=0.0).reshape(
                              nSignalPoints*covMatrixSide**2).tolist(),
                          'A_MATRIX': np.array(self.aMatrix[ampName]).reshape(covMatrixSide**2).tolist(),
                          'B_MATRIX': np.array(self.bMatrix[ampName]).reshape(covMatrixSide**2).tolist(),
                          'COVARIANCES_MODEL_NO_B':
                              np.pad(np.array(self.covariancesModelNoB[ampName]),
                                     ((0, nPadPoints[ampName]), (0, 0), (0, 0)),
                                     'constant', constant_values=np.nan).reshape(
                              nSignalPoints*covMatrixSide**2).tolist(),
                          'A_MATRIX_NO_B': np.array(self.aMatrixNoB[ampName]).reshape(
                              covMatrixSide**2).tolist(),
                          'FINAL_VARS': np.pad(np.array(self.finalVars[ampName]), (0, nPadPoints[ampName]),
                                               'constant', constant_values=np.nan).tolist(),
                          'FINAL_MODEL_VARS': np.pad(np.array(self.finalModelVars[ampName]),
                                                     (0, nPadPoints[ampName]),
                                                     'constant', constant_values=np.nan).tolist(),
                          'FINAL_MEANS': np.pad(np.array(self.finalMeans[ampName]),
                                                (0, nPadPoints[ampName]),
                                                'constant', constant_values=np.nan).tolist(),
                          'BAD_AMPS': np.array(self.badAmps).tolist() if len(self.badAmps) else np.nan,
                          'PHOTO_CHARGE': np.array(self.photoCharge[ampName]).tolist(),
                          } for ampName in self.ampNames])
        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta
        tableList.append(catalog)

        return(tableList)

    def getExpIdsUsed(self, ampName):
        """Get the exposures used, i.e. not discarded, for a given amp.
        If no mask has been created yet, all exposures are returned.
        """
        if len(self.expIdMask[ampName]) == 0:
            return self.inputExpIdPairs[ampName]

        # if the mask exists it had better be the same length as the expIdPairs
        assert len(self.expIdMask[ampName]) == len(self.inputExpIdPairs[ampName])

        pairs = self.inputExpIdPairs[ampName]
        mask = self.expIdMask[ampName]
        # cast to bool required because numpy
        return [(exp1, exp2) for ((exp1, exp2), m) in zip(pairs, mask) if bool(m) is True]

    def getGoodAmps(self):
        return [amp for amp in self.ampNames if amp not in self.badAmps]
