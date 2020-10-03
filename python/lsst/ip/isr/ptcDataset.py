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
        Dictionary keyed by amp names containing the errors on the fitted noise.
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
    covariancesNoB : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing a list of measured
        covariances per mean flux ('b'=0 in
        Astier+19).
    covariancesModelNoB : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing covariances model
        (with 'b'=0 in Eq. 20 of Astier+19)
        per mean flux.
    covariancesSqrtWeightsNoB : `dict`, [`str`, `list`]
        Dictionary keyed by amp names containing sqrt. of covariances weights
        ('b' = 0 in Eq. 20 of
        Astier+19).
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
    """

    _OBSTYPE = 'PTC'
    _SCHEMA = 'Gen3 Photon Transfer Curve'
    _VERSION = 1.0

    def __init__(self, ampNames=[], ptcFitType=None, **kwargs):

        self.ptcFitType = ptcFitType
        self.ampNames = ampNames
        self.badAmps = []

        self.inputExpIdPairs = {ampName: [] for ampName in ampNames}
        self.expIdMask = {ampName: [] for ampName in ampNames}
        self.rawExpTimes = {ampName: [] for ampName in ampNames}
        self.rawMeans = {ampName: [] for ampName in ampNames}
        self.rawVars = {ampName: [] for ampName in ampNames}
        self.photoCharge = {ampName: [] for ampName in ampNames}

        self.gain = {ampName: -1. for ampName in ampNames}
        self.gainErr = {ampName: -1. for ampName in ampNames}
        self.noise = {ampName: -1. for ampName in ampNames}
        self.noiseErr = {ampName: -1. for ampName in ampNames}

        self.ptcFitPars = {ampName: [] for ampName in ampNames}
        self.ptcFitParsError = {ampName: [] for ampName in ampNames}
        self.ptcFitChiSq = {ampName: [] for ampName in ampNames}

        self.covariances = {ampName: [] for ampName in ampNames}
        self.covariancesModel = {ampName: [] for ampName in ampNames}
        self.covariancesSqrtWeights = {ampName: [] for ampName in ampNames}
        self.aMatrix = {ampName: [] for ampName in ampNames}
        self.bMatrix = {ampName: [] for ampName in ampNames}
        self.covariancesNoB = {ampName: [] for ampName in ampNames}
        self.covariancesModelNoB = {ampName: [] for ampName in ampNames}
        self.covariancesSqrtWeightsNoB = {ampName: [] for ampName in ampNames}
        self.aMatrixNoB = {ampName: [] for ampName in ampNames}

        self.finalVars = {ampName: [] for ampName in ampNames}
        self.finalModelVars = {ampName: [] for ampName in ampNames}
        self.finalMeans = {ampName: [] for ampName in ampNames}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['badAmps', 'inputExpIdPairs', 'expIdMask', 'rawExpTimes',
                                        'rawMeans', 'rawVars', 'gain', 'gainErr', 'noise', 'noiseErr',
                                        'ptcFitPars', 'ptcFitParsError', 'ptcFitChiSq', 'aMatrixNoB',
                                        'covariances', 'covariancesModel', 'covariancesSqrtWeights',
                                        'covariancesNoB', 'covariancesModelNoB', 'covariancesSqrtWeightsNoB',
                                        'aMatrix', 'bMatrix', 'finalVars', 'finalModelVars', 'finalMeans',
                                        'photoCharge'])

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
        calib.badAmps = dictionary['badAmps']
        if calib.ptcFitType in ['FULLCOVARIANCE', ]:
            # Number of final signal levels
            nSignalPoints = len(list(dictionary['ampNames'].values())[0]['rawMeans'])
            # The cov matrices are square
            covMatrixSide = int(np.sqrt(len(
                                list(dictionary['ampNames'].values())[0]['covariances'])/nSignalPoints))
        for ampName in dictionary['ampNames']:
            calib.ampNames.append(ampName)
            calib.inputExpIdPairs.setdefault(ampName, dictionary['inputExpIdPairs'][ampName])
            calib.expIdMask.setdefault(ampName, dictionary['expIdMask'][ampName])
            calib.rawExpTimes.setdefault(ampName, dictionary['rawExpTimes'][ampName])
            calib.rawMeans.setdefault(ampName, dictionary['rawMeans'][ampName])
            calib.rawVars.setdefault(ampName, dictionary['rawVars'][ampName])
            calib.gain.setdefault(ampName, dictionary['gain'][ampName])
            calib.gainErr.setdefault(ampName, dictionary['gainErr'][ampName])
            calib.noise.setdefault(ampName, dictionary['noise'][ampName])
            calib.noiseErr.setdefault(ampName, dictionary['noiseErr'][ampName])
            calib.ptcFitPars.setdefault(ampName, dictionary['ptcFitPars'][ampName])
            calib.ptcFitParsError.setdefault(ampName, dictionary['ptcFitParsError'][ampName])
            calib.ptcFitChiSq.setdefault(ampName, dictionary['ptcFitChiSq'][ampName])
            if calib.ptcFitType in ['FULLCOVARIANCE', ]:
                calib.covariances.setdefault(ampName, dictionary['covariances'][ampName].reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide)))
                calib.covariancesModel.setdefault(ampName, dictionary['covariancesModel'][ampName].reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide)))
                calib.covariancesSqrtWeights.setdefault(
                    ampName, dictionary['covariancesSqrtWeights'][ampName].reshape(
                        (nSignalPoints, covMatrixSide, covMatrixSide)))
                calib.aMatrix.setdefault(ampName, dictionary['aMatrix'][ampName].reshape(
                    (covMatrixSide, covMatrixSide)))
                calib.bMatrix.setdefault(ampName, dictionary['bMatrix'][ampName].reshape(
                    (covMatrixSide, covMatrixSide)))
                calib.covariancesNoB.setdefault(ampName, dictionary['covariancesNoB'][ampName].reshape(
                    (nSignalPoints, covMatrixSide, covMatrixSide)))
                calib.covariancesModelNoB.setdefault(ampName,
                                                     dictionary['covariancesModelNoB'][ampName].reshape(
                                                         (nSignalPoints, covMatrixSide, covMatrixSide)))
                calib.covariancesSqrtWeightsNoB.setdefault(
                    ampName, dictionary['covariancesSqrtWeightsNoB'][ampName].reshape(
                        (nSignalPoints, covMatrixSide, covMatrixSide)))
                calib.aMatrixNoB.setdefault(ampName,
                                            dictionary['aMatrixNoB'][ampName].reshape(
                                                (covMatrixSide, covMatrixSide)))
            else:
                calib.covariances.setdefault(ampName, [])
                calib.covariancesModel.setdefault(ampName, [])
                calib.covariancesSqrtWeights.setdefault(ampName, [])
                calib.aMatrix.setdefault(ampName, [])
                calib.bMatrix.setdefault(ampName, [])
                calib.covariancesNoB.setdefault(ampName, [])
                calib.covariancesModelNoB.setdefault(ampName, [])
                calib.covariancesSqrtWeightsNoB.setdefault(ampName, [])
                calib.aMatrixNoB.setdefault(ampName, [])
            calib.finalVars.setdefault(ampName, dictionary['finalVars'][ampName])
            calib.finalModelVars.setdefault(ampName, dictionary['finalModelVars'][ampName])
            calib.finalMeans.setdefault(ampName, dictionary['finalMeans'][ampName])
            calib.photoCharge.setdefault(ampName, dictionary['photoCharge'][ampName])

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
        outDict['covariances'] = self.covariances
        outDict['covariancesModel'] = self.covariancesModel
        outDict['covariancesSqrtWeights'] = self.covariancesSqrtWeights
        outDict['aMatrix'] = self.aMatrix
        outDict['bMatrix'] = self.bMatrix
        outDict['covariancesNoB'] = self.covariancesNoB
        outDict['covariancesModelNoB'] = self.covariancesModelNoB
        outDict['covariancesSqrtWeightsNoB'] = self.covariancesSqrtWeightsNoB
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
        inDict['covariances'] = dict()
        inDict['covariancesModel'] = dict()
        inDict['covariancesSqrtWeights'] = dict()
        inDict['aMatrix'] = dict()
        inDict['bMatrix'] = dict()
        inDict['covariancesNoB'] = dict()
        inDict['covariancesModelNoB'] = dict()
        inDict['covariancesSqrtWeightsNoB'] = dict()
        inDict['aMatrixNoB'] = dict()
        inDict['finalVars'] = dict()
        inDict['finalModelVars'] = dict()
        inDict['finalMeans'] = dict()
        inDict['badAmps'] = dict()
        inDict['photoCharge'] = dict()

        for record in ptcTable:
            ampName = record['AMPLIFIER_NAME']

            inDict['ptcFitType'] = record['PTC_FIT_TYPE']
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
            inDict['covariancesNoB'][ampName] = record['COVARIANCES_NO_B']
            inDict['covariancesModelNoB'][ampName] = record['COVARIANCES_MODEL_NO_B']
            inDict['covariancesSqrtWeightsNoB'][ampName] = record['COVARIANCES_SQRT_WEIGHTS_NO_B']
            inDict['aMatrixNoB'][ampName] = record['A_MATRIX_NO_B']
            inDict['finalVars'][ampName] = record['FINAL_VARS']
            inDict['finalModelVars'][ampName] = record['FINAL_MODEL_VARS']
            inDict['finalMeans'][ampName] = record['FINAL_MEANS']
            inDict['badAmps'][ampName] = record['BAD_AMPS'] if 'BAD_AMPS' in record else []
            inDict['photoCharge'][ampName] = record['PHOTO_CHARGE']

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

        nSignalPoints = len(list(self.rawExpTimes.values())[0])
        covMatrixSide = list(self.aMatrix.values())[0].shape[0]

        catalog = Table([{'AMPLIFIER_NAME': ampName,
                          'PTC_FIT_TYPE': self.ptcFitType,
                          'INPUT_EXP_ID_PAIRS': self.inputExpIdPairs[ampName],
                          'EXP_ID_MASK': self.expIdMask[ampName],
                          'RAW_EXP_TIMES': self.rawExpTimes[ampName],
                          'RAW_MEANS': self.rawMeans[ampName],
                          'RAW_VARS': self.rawVars[ampName],
                          'GAIN': self.gain[ampName],
                          'GAIN_ERR': self.gainErr[ampName],
                          'NOISE': self.noise[ampName],
                          'NOISE_ERR': self.noiseErr[ampName],
                          'PTC_FIT_PARS': self.ptcFitPars[ampName],
                          'PTC_FIT_PARS_ERROR': self.ptcFitParsError[ampName],
                          'PTC_FIT_CHI_SQ': self.ptcFitChiSq[ampName],
                          'COVARIANCES': np.array(self.covariances[ampName]).reshape(
                              nSignalPoints*covMatrixSide**2),
                          'COVARIANCES_MODEL':
                              np.array(self.covariancesModel[ampName]).reshape(
                                  nSignalPoints*covMatrixSide**2),
                          'COVARIANCES_SQRT_WEIGHTS':
                              np.array(self.covariancesSqrtWeights[ampName]).reshape(
                                  nSignalPoints*covMatrixSide**2),
                          'A_MATRIX': np.array(self.aMatrix[ampName]).reshape(covMatrixSide**2),
                          'B_MATRIX': np.array(self.bMatrix[ampName]).reshape(covMatrixSide**2),
                          'COVARIANCES_NO_B':
                              np.array(self.covariancesNoB[ampName]).reshape(
                                  nSignalPoints*covMatrixSide**2),
                          'COVARIANCES_MODEL_NO_B':
                              np.array(self.covariancesModelNoB[ampName]).reshape(
                                  nSignalPoints*covMatrixSide**2),
                          'COVARIANCES_SQRT_WEIGHTS_NO_B':
                              np.array(self.covariancesSqrtWeightsNoB[ampName]).reshape(
                                  nSignalPoints*covMatrixSide**2),
                          'A_MATRIX_NO_B': np.array(self.aMatrixNoB[ampName]).reshape(covMatrixSide**2),
                          'FINAL_VARS': self.finalVars[ampName],
                          'FINAL_MODEL_VARS': self.finalModelVars[ampName],
                          'FINAL_MEANS': self.finalMeans[ampName],
                          'BAD_AMPS': self.badAmps if len(self.badAmps) else np.nan,
                          'PHOTO_CHARGE': self.photoCharge[ampName]
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
