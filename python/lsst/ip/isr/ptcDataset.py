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

import numbers
import numpy as np
import math
from astropy.table import Table
import numpy.polynomial.polynomial as poly
from scipy.signal import fftconvolve

from lsst.ip.isr import IsrCalib


def symmetrize(inputArray):
    """ Copy array over 4 quadrants prior to convolution.

    Parameters
    ----------
    inputarray : `numpy.array`
        Input array to symmetrize.

    Returns
    -------
    aSym : `numpy.array`
        Symmetrized array.
    """

    targetShape = list(inputArray.shape)
    r1, r2 = inputArray.shape[-1], inputArray.shape[-2]
    targetShape[-1] = 2*r1-1
    targetShape[-2] = 2*r2-1
    aSym = np.ndarray(tuple(targetShape))
    aSym[..., r2-1:, r1-1:] = inputArray
    aSym[..., r2-1:, r1-1::-1] = inputArray
    aSym[..., r2-1::-1, r1-1::-1] = inputArray
    aSym[..., r2-1::-1, r1-1:] = inputArray

    return aSym


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
    covMatrixSideFullCovFit : `int`, optional
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
    inputExpPairMjdStartList : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the start mjd from
        the first exposure in each flat pair.
    expIdMask : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the mask produced after
        outlier rejection. The mask produced by the "FULLCOVARIANCE"
        option may differ from the one produced in the other two PTC
        fit types.
    expIdRolloffMask : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the mask produced after
        outlier rejection and extended to fit the PTC rolloff.
    rawExpTimes : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the unmasked exposure times.
    rawMeans : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the unmasked average of the
        means of the exposures in each flat pair (units: adu).
    rawVars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the variance of the
        difference image of the exposures in each flat pair (units: adu^2).
    rawDeltas : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the scaled unmasked delta of
        the means of the exposures in each flat pair (mean2 - mean1)
        (units: adu).
    rowMeanVariance : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the variance of the
        means of the rows of the difference image of the exposures
        in each flat pair (units: adu^2).
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
        Dictionary keyed by amp names containing the fitted gains. May be
        adjusted by amp-offset gain ratios if configured in PTC solver.
    gainUnadjusted : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing unadjusted (raw) fit gain
        values. May be the same as gain values if amp-offset adjustment
        is not turned on.
    gainErr : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the errors on the
        fitted gains.
    gainList : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the gain estimated from
        each flat pair.
    overscanMedianLevelList : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the median overscan
        level from each input flat pair (units: adu).
    overscanMedian : `dict `, [`str`, `float`]
        Dictionary keyed by amp names containing the median of
        overscanMedianLevelList[expIdMask] (units: adu).
    overscanMedianSigma : `dict `, [`str`, `float`]
        Dictionary keyed by amp names containing the median absolute
        deviation of overscanMedianLevelList[expIdMask] (units: adu).
    noiseList : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the mean overscan
        standard deviation from each flat pair (units: adu).
    noise : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the fitted noise
        (units: electron).
    noiseErr : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the errors on the fitted
        noise (units: electron).
    ampOffsets : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing amp-to-amp offsets
        (units: adu).
    ptcFitPars : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the fitted parameters of the
        PTC model for ptcFitType in ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcFitParsError : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the errors on the fitted
        parameters of the PTC model for ptcFitType in
        ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcFitChiSq : `dict`, [`str`, `float`]
        Dictionary keyed by amp names containing the reduced chi squared
        of the fit for ptcFitType in ["POLYNOMIAL", "EXPAPPROXIMATION"].
    ptcTurnoff : `dict` [`str, `float`]
        Flux value (in adu) where the variance of the PTC curve starts
        decreasing consistently.
    ptcTurnoffSamplingError : `dict` [`str`, `float`]
        ``Sampling`` error on the ptcTurnoff, based on the flux sampling
        of the input PTC (units: adu).
    ptcRolloff : `dict` [`str, `float`]
        Flux value (in adu) where the variance of the PTC curve begins starts
        to move away from the PTC model by some threshold.
    ptcRolloffError : `dict` [`str`, `float`]
        Covariance error from the fit to the PTC rolloff (units: adu).
    ptcRolloff : `dict` [`str, `float`]
        Curvature parameter of the PTC rolloff.
    ptcRolloffError : `dict` [`str`, `float`]
        Fitting error of the curvature parameter of the PTC rolloff.
    nPixelCovariances : `dict`, [`str`, `int`]
        Dictionary keyed by amp names containing the number of pixels
        that were used to measure the covariances.
    covariances : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing a list of measured
        covariances per mean flux (units: adu^2).
    covariancesModel : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containinging covariances model
        (Eq. 20 of Astier+19) per mean flux (units: adu^2).
    covariancesSqrtWeights : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containinging sqrt. of covariances
        weights (units: 1/adu).
    aMatrix : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "a" parameters from
        the model in Eq. 20 of Astier+19 (units: 1/electron).
    bMatrix : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "b" parameters from
        the model in Eq. 20 of Astier+19 (units: 1/electron).
    noiseMatrix : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "noise" parameters from
        the model in Eq. 20 of Astier+19 (units: electron^2).
    covariancesModelNoB : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing covariances model
        (with 'b'=0 in Eq. 20 of Astier+19) per mean flux (units:
        adu^2). Will be removed after v29.
    aMatrixNoB : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "a" parameters from the
        model in Eq. 20 of Astier+19 (and 'b' = 0) (units: 1/electron).
        Will be removed after v29.
    noiseMatrixNoB : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the "noise" parameters from
        the model in Eq. 20 of Astier+19, with 'b' = 0 (units: electron^2).
        Will be removed after v29.
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
    photoChargeDeltas : `dict`, [`str`, `np.ndarray`]
        Dictionary keyed by amp names containing the delta for the integrated
        photocharge (photocharge2 - photocharge1).
    auxValues : `dict`, [`str`, `np.ndarray`]
        Dictionary of per-detector auxiliary header values that can be used
        for PTC, linearity computation.

    Version 1.1 adds the `ptcTurnoff` attribute.
    Version 1.2 adds the `histVars`, `histChi2Dofs`, and `kspValues`
        attributes.
    Version 1.3 adds the `noiseMatrix` and `noiseMatrixNoB` attributes.
    Version 1.4 adds the `auxValues` attribute.
    Version 1.5 adds the `covMatrixSideFullCovFit` attribute.
    Version 1.6 adds the `rowMeanVariance` attribute.
    Version 1.7 adds the `noiseList` attribute.
    Version 1.8 adds the `ptcTurnoffSamplingError` attribute.
    Version 1.9 standardizes PTC noise units to electron.
    Version 2.0 adds the `ampOffsets`, `gainUnadjusted`, and
        `gainList` attributes.
    Version 2.1 deprecates the `covariancesModelNoB`, `aMatrixNoB`, and
        `noiseMatrixNoB` attributes.
    Version 2.2 adds the `overscanMedianLevelList` and
        `inputExpPairMjdStartList` attributes.
    Version 2.3 adds the `overscanMedian` and
        `overscanMedianSigma` attrbutes.
    Version 2.4 adds the `nPixelCovariances` attribute.
    Version 2.5 adds the `rawDeltas` and `photoChargeDeltas` attributes.
    Version 2.6 adds the `expIdRolloffMask`, `ptcRolloff`, `ptcRolloffError`,
        `ptcRolloffTau`, and `ptcRolloffTauError` attributes.
    """

    _OBSTYPE = 'PTC'
    _SCHEMA = 'Gen3 Photon Transfer Curve'
    # When adding a new field to update the version, be sure to update the
    # following methods:
    #  * __init__()
    #  * fromDict()
    #  * toDict()
    #  * fromTable()
    #  * toTable()
    #  * setAmpValuesPartialDataset()
    #  * appendPartialPtc()
    #  * sort()
    #  * _checkTypes() in test_ptcDataset.py
    #  * test_ptcDataset() in test_ptcDataset.py
    #  * test_ptcDatasetSort in test_ptcDataset.py
    #  * test_ptcDatasetAppend in test_ptcDataset.py
    _VERSION = 2.6

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
        self.inputExpPairMjdStartList = {ampName: np.array([]) for ampName in ampNames}
        self.expIdMask = {ampName: np.array([], dtype=bool) for ampName in ampNames}
        self.expIdRolloffMask = {ampName: np.array([], dtype=bool) for ampName in ampNames}
        self.rawExpTimes = {ampName: np.array([]) for ampName in ampNames}
        self.rawMeans = {ampName: np.array([]) for ampName in ampNames}
        self.rawVars = {ampName: np.array([]) for ampName in ampNames}
        self.rawDeltas = {ampName: np.array([]) for ampName in ampNames}
        self.rowMeanVariance = {ampName: np.array([]) for ampName in ampNames}
        self.photoCharges = {ampName: np.array([]) for ampName in ampNames}
        self.photoChargeDeltas = {ampName: np.array([]) for ampName in ampNames}
        self.ampOffsets = {ampName: np.array([]) for ampName in ampNames}

        self.gain = {ampName: np.nan for ampName in ampNames}
        self.gainUnadjusted = {ampName: np.nan for ampName in ampNames}
        self.gainErr = {ampName: np.nan for ampName in ampNames}
        self.gainList = {ampName: np.array([]) for ampName in ampNames}
        self.overscanMedianLevelList = {ampName: np.array([]) for ampName in ampNames}
        self.overscanMedian = {ampName: np.nan for ampName in ampNames}
        self.overscanMedianSigma = {ampName: np.nan for ampName in ampNames}
        self.noiseList = {ampName: np.array([]) for ampName in ampNames}
        self.noise = {ampName: np.nan for ampName in ampNames}
        self.noiseErr = {ampName: np.nan for ampName in ampNames}

        self.histVars = {ampName: np.array([]) for ampName in ampNames}
        self.histChi2Dofs = {ampName: np.array([]) for ampName in ampNames}
        self.kspValues = {ampName: np.array([]) for ampName in ampNames}

        self.ptcFitPars = {ampName: np.array([]) for ampName in ampNames}
        self.ptcFitParsError = {ampName: np.array([]) for ampName in ampNames}
        self.ptcFitChiSq = {ampName: np.nan for ampName in ampNames}
        self.ptcTurnoff = {ampName: np.nan for ampName in ampNames}
        self.ptcTurnoffSamplingError = {ampName: np.nan for ampName in ampNames}
        self.ptcRolloff = {ampName: np.nan for ampName in ampNames}
        self.ptcRolloffError = {ampName: np.nan for ampName in ampNames}
        self.ptcRolloffTau = {ampName: np.nan for ampName in ampNames}
        self.ptcRolloffTauError = {ampName: np.nan for ampName in ampNames}

        self.nPixelCovariances = {ampName: -1 for ampName in ampNames}
        self.covariances = {ampName: np.array([]) for ampName in ampNames}
        self.covariancesModel = {ampName: np.array([]) for ampName in ampNames}
        self.covariancesSqrtWeights = {ampName: np.array([]) for ampName in ampNames}
        self.aMatrix = {ampName: np.array([]) for ampName in ampNames}
        self.bMatrix = {ampName: np.array([]) for ampName in ampNames}
        self.noiseMatrix = {ampName: np.array([]) for ampName in ampNames}

        self.finalVars = {ampName: np.array([]) for ampName in ampNames}
        self.finalModelVars = {ampName: np.array([]) for ampName in ampNames}
        self.finalMeans = {ampName: np.array([]) for ampName in ampNames}

        # Try this as a dict of arrays.
        self.auxValues = {}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['badAmps', 'inputExpIdPairs', 'inputExpPairMjdStartList',
                                        'expIdMask', 'expIdRolloffMask', 'rawExpTimes', 'rawMeans',
                                        'rawVars', 'rowMeanVariance', 'gain', 'gainErr', 'gainList',
                                        'noise', 'noiseErr', 'noiseList', 'overscanMedianLevelList',
                                        'overscanMedian', 'overscanMedianSigma', 'ptcFitPars',
                                        'ptcFitParsError', 'ptcFitChiSq', 'ptcTurnoff', 'covariances',
                                        'covariancesModel', 'covariancesSqrtWeights', 'aMatrix',
                                        'bMatrix', 'noiseMatrix', 'finalVars', 'finalModelVars',
                                        'finalMeans', 'photoCharges', 'histVars', 'histChi2Dofs',
                                        'kspValues', 'auxValues', 'ptcTurnoffSamplingError',
                                        'ptcRolloff', 'ptcRolloffError', 'ptcRolloffTau',
                                        'ptcRolloffTauError', 'ampOffsets', 'gainUnadjusted',
                                        'nPixelCovariances', 'rawDeltas', 'photoChargeDeltas'])

        self.updateMetadata(setCalibInfo=True, setCalibId=True, **kwargs)
        self._validateCovarianceMatrizSizes()

    def setAmpValuesPartialDataset(
            self,
            ampName,
            inputExpIdPair=(-1, -1),
            inputExpPairMjdStart=np.nan,
            rawExpTime=np.nan,
            rawMean=np.nan,
            rawVar=np.nan,
            rawDelta=np.nan,
            rowMeanVariance=np.nan,
            photoCharge=np.nan,
            photoChargeDelta=np.nan,
            ampOffset=np.nan,
            expIdMask=False,
            nPixelCovariance=-1,
            covariance=None,
            covSqrtWeights=None,
            gain=np.nan,
            noise=np.nan,
            overscanMedianLevel=np.nan,
            histVar=np.nan,
            histChi2Dof=np.nan,
            kspValue=0.0,
    ):
        """
        Set the amp values for a partial PTC Dataset (from cpExtractPtcTask).

        Parameters
        ----------
        ampName : `str`
            Name of the amp to set the values.
        inputExpIdPair : `tuple` [`int`]
            Exposure IDs of input pair.
        inputExpPairMjdStart : `float`, optional
            The start MJD of first exposure in the flat pair.
        rawExpTime : `float`, optional
            Exposure time for this exposure pair (units: sec).
        rawMean : `float`, optional
            Average of the means of the exposures in this pair
            (units: adu).
        rawVar : `float`, optional
            Variance of the difference of the exposures in this pair
            (units: adu^2).
        rawDelta : `float`, optional
            Delta of the means of the exposure in this pair
            (units: adu).
        rowMeanVariance : `float`, optional
            Variance of the means of the rows in the difference image
            of the exposures in this pair (units: adu^2).
        photoCharge : `float`, optional
            Integrated photocharge for flat pair for linearity calibration
            (arbitrary units).
        photoChargeDelta : `float`, optional
            Delta between integrated photocharge for the flat pair
            (arbitrary units).
        ampOffset : `float`, optional
            Amp offset for this amplifier.
        expIdMask : `bool`, optional
            Flag setting if this exposure pair should be used (True)
            or not used (False).
        nPixelCovariance : `int`, optional
            Number of pixels that went into the covariance measurement.
        covariance : `np.ndarray` or None, optional
            Measured covariance for this exposure pair (units: adu^2).
        covSqrtWeights : `np.ndarray` or None, optional
            Measured sqrt of covariance weights in this exposure pair
            (units: 1/adu).
        gain : `float`, optional
            Estimated gain for this exposure pair (units: electron/adu).
        noise : `float`, optional
            Estimated read noise for this exposure pair (units: electron).
        overscanMedianLevel : `float`, optional
            Average of the median overscan levels for this exposure pair.
            (units: adu)
        histVar : `float`, optional
            Variance estimated from fitting a histogram with a Gaussian model
            (units: adu).
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
        self.inputExpPairMjdStartList[ampName] = np.array([inputExpPairMjdStart])
        self.rawExpTimes[ampName] = np.array([rawExpTime])
        self.rawMeans[ampName] = np.array([rawMean])
        self.rawVars[ampName] = np.array([rawVar])
        self.rawDeltas[ampName] = np.array([rawDelta])
        self.rowMeanVariance[ampName] = np.array([rowMeanVariance])
        self.photoCharges[ampName] = np.array([photoCharge])
        self.photoChargeDeltas[ampName] = np.array([photoChargeDelta])
        self.ampOffsets[ampName] = np.array([ampOffset])
        self.expIdMask[ampName] = np.array([expIdMask])
        self.nPixelCovariances[ampName] = nPixelCovariance
        self.covariances[ampName] = np.array([covariance])
        self.covariancesSqrtWeights[ampName] = np.array([covSqrtWeights])
        self.gain[ampName] = gain
        self.gainUnadjusted[ampName] = gain
        self.gainList[ampName] = np.array([gain])
        self.noise[ampName] = noise
        self.noiseList[ampName] = np.array([noise])
        self.overscanMedianLevelList[ampName] = np.array([overscanMedianLevel])
        self.overscanMedian[ampName] = float(overscanMedianLevel)
        self.overscanMedianSigma[ampName] = float(0.0)
        self.histVars[ampName] = np.array([histVar])
        self.histChi2Dofs[ampName] = np.array([histChi2Dof])
        self.kspValues[ampName] = np.array([kspValue])

        # From FULLCOVARIANCE model
        self.covariancesModel[ampName] = np.array([nanMatrixFit])
        self.aMatrix[ampName] = nanMatrixFit
        self.bMatrix[ampName] = nanMatrixFit
        self.noiseMatrix[ampName] = nanMatrixFit

        # Filler values.
        self.finalVars[ampName] = np.array([np.nan])
        self.finalModelVars[ampName] = np.array([np.nan])
        self.finalMeans[ampName] = np.array([np.nan])

    def setAuxValuesPartialDataset(self, auxDict):
        """
        Set a dictionary of auxiliary values for a partial dataset.

        Parameters
        ----------
        auxDict : `dict` [`str`, `float`]
            Dictionary of float values.
        """
        for key, value in auxDict.items():
            if isinstance(value, numbers.Integral):
                self.auxValues[key] = np.atleast_1d(np.asarray(value).astype(np.int64))
            elif isinstance(value, (str, np.str_, np.bytes_)):
                self.auxValues[key] = np.atleast_1d(np.asarray(value))
            else:
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
            calib.inputExpPairMjdStartList[ampName] = np.array(
                dictionary['inputExpPairMjdStartList'][ampName],
            )
            calib.expIdMask[ampName] = np.array(dictionary['expIdMask'][ampName])
            calib.expIdRolloffMask[ampName] = np.array(dictionary['expIdRolloffMask'][ampName])
            calib.rawExpTimes[ampName] = np.array(dictionary['rawExpTimes'][ampName], dtype=np.float64)
            calib.rawMeans[ampName] = np.array(dictionary['rawMeans'][ampName], dtype=np.float64)
            calib.rawVars[ampName] = np.array(dictionary['rawVars'][ampName], dtype=np.float64)
            calib.rawDeltas[ampName] = np.array(dictionary['rawDeltas'][ampName], dtype=np.float64)
            calib.rowMeanVariance[ampName] = np.array(dictionary['rowMeanVariance'][ampName],
                                                      dtype=np.float64)
            calib.gain[ampName] = float(dictionary['gain'][ampName])
            calib.gainErr[ampName] = float(dictionary['gainErr'][ampName])
            calib.gainUnadjusted[ampName] = float(dictionary['gainUnadjusted'][ampName])
            calib.gainList[ampName] = np.array(dictionary['gainList'][ampName], dtype=np.float64)
            calib.noiseList[ampName] = np.array(dictionary['noiseList'][ampName], dtype=np.float64)
            calib.nPixelCovariances[ampName] = int(dictionary['nPixelCovariances'][ampName])
            calib.overscanMedianLevelList[ampName] = np.array(
                dictionary['overscanMedianLevelList'][ampName],
                dtype=np.float64,
            )
            calib.overscanMedian[ampName] = float(dictionary['overscanMedian'][ampName])
            calib.overscanMedianSigma[ampName] = float(dictionary['overscanMedianSigma'][ampName])
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
            calib.ptcTurnoffSamplingError[ampName] = float(dictionary['ptcTurnoffSamplingError'][ampName])
            calib.ptcRolloff[ampName] = float(dictionary['ptcRolloff'][ampName])
            calib.ptcRolloffError[ampName] = float(dictionary['ptcRolloffError'][ampName])
            calib.ptcRolloffTau[ampName] = float(dictionary['ptcRolloff'][ampName])
            calib.ptcRolloffTauError[ampName] = float(dictionary['ptcRolloffError'][ampName])
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
                calib.noiseMatrix[ampName] = np.array(
                    dictionary['noiseMatrix'][ampName],
                    dtype=np.float64).reshape((covMatrixSideFullCovFit, covMatrixSideFullCovFit))
            else:
                # Empty dataset
                calib.covariances[ampName] = np.array([], dtype=np.float64)
                calib.covariancesModel[ampName] = np.array([], dtype=np.float64)
                calib.covariancesSqrtWeights[ampName] = np.array([], dtype=np.float64)
                calib.aMatrix[ampName] = np.array([], dtype=np.float64)
                calib.bMatrix[ampName] = np.array([], dtype=np.float64)
                calib.noiseMatrix[ampName] = np.array([], dtype=np.float64)

            calib.finalVars[ampName] = np.array(dictionary['finalVars'][ampName], dtype=np.float64)
            calib.finalModelVars[ampName] = np.array(dictionary['finalModelVars'][ampName], dtype=np.float64)
            calib.finalMeans[ampName] = np.array(dictionary['finalMeans'][ampName], dtype=np.float64)
            calib.photoCharges[ampName] = np.array(dictionary['photoCharges'][ampName], dtype=np.float64)
            calib.photoChargeDeltas[ampName] = np.array(
                dictionary['photoChargeDeltas'][ampName],
                dtype=np.float64,
            )
            calib.ampOffsets[ampName] = np.array(dictionary['ampOffsets'][ampName], dtype=np.float64)

        for key, value in dictionary['auxValues'].items():
            if isinstance(value[0], numbers.Integral):
                calib.auxValues[key] = np.atleast_1d(np.asarray(value).astype(np.int64))
            elif isinstance(value[0], (str, np.str_, np.bytes_)):
                calib.auxValues[key] = np.atleast_1d(np.asarray(value))
            else:
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
        outDict['inputExpPairMjdStartList'] = _dictOfArraysToDictOfLists(self.inputExpPairMjdStartList)
        outDict['expIdMask'] = _dictOfArraysToDictOfLists(self.expIdMask)
        outDict['expIdRolloffMask'] = _dictOfArraysToDictOfLists(self.expIdRolloffMask)
        outDict['rawExpTimes'] = _dictOfArraysToDictOfLists(self.rawExpTimes)
        outDict['rawMeans'] = _dictOfArraysToDictOfLists(self.rawMeans)
        outDict['rawVars'] = _dictOfArraysToDictOfLists(self.rawVars)
        outDict['rawDeltas'] = _dictOfArraysToDictOfLists(self.rawDeltas)
        outDict['rowMeanVariance'] = _dictOfArraysToDictOfLists(self.rowMeanVariance)
        outDict['gain'] = self.gain
        outDict['gainErr'] = self.gainErr
        outDict['gainUnadjusted'] = self.gainUnadjusted
        outDict['gainList'] = _dictOfArraysToDictOfLists(self.gainList)
        outDict['noiseList'] = _dictOfArraysToDictOfLists(self.noiseList)
        outDict['overscanMedianLevelList'] = _dictOfArraysToDictOfLists(self.overscanMedianLevelList)
        outDict['overscanMedian'] = self.overscanMedian
        outDict['overscanMedianSigma'] = self.overscanMedianSigma
        outDict['noise'] = self.noise
        outDict['noiseErr'] = self.noiseErr
        outDict['histVars'] = self.histVars
        outDict['histChi2Dofs'] = self.histChi2Dofs
        outDict['kspValues'] = self.kspValues
        outDict['ptcFitPars'] = _dictOfArraysToDictOfLists(self.ptcFitPars)
        outDict['ptcFitParsError'] = _dictOfArraysToDictOfLists(self.ptcFitParsError)
        outDict['ptcFitChiSq'] = self.ptcFitChiSq
        outDict['ptcTurnoff'] = self.ptcTurnoff
        outDict['ptcTurnoffSamplingError'] = self.ptcTurnoffSamplingError
        outDict['ptcRolloff'] = self.ptcRolloff
        outDict['ptcRolloffError'] = self.ptcRolloffError
        outDict['ptcRolloffTau'] = self.ptcRolloff
        outDict['ptcRolloffTauError'] = self.ptcRolloffTauError
        outDict['nPixelCovariances'] = self.nPixelCovariances
        outDict['covariances'] = _dictOfArraysToDictOfLists(self.covariances)
        outDict['covariancesModel'] = _dictOfArraysToDictOfLists(self.covariancesModel)
        outDict['covariancesSqrtWeights'] = _dictOfArraysToDictOfLists(self.covariancesSqrtWeights)
        outDict['aMatrix'] = _dictOfArraysToDictOfLists(self.aMatrix)
        outDict['bMatrix'] = _dictOfArraysToDictOfLists(self.bMatrix)
        outDict['noiseMatrix'] = _dictOfArraysToDictOfLists(self.noiseMatrix)
        outDict['finalVars'] = _dictOfArraysToDictOfLists(self.finalVars)
        outDict['finalModelVars'] = _dictOfArraysToDictOfLists(self.finalModelVars)
        outDict['finalMeans'] = _dictOfArraysToDictOfLists(self.finalMeans)
        outDict['photoCharges'] = _dictOfArraysToDictOfLists(self.photoCharges)
        outDict['photoChargeDeltas'] = _dictOfArraysToDictOfLists(self.photoChargeDeltas)
        outDict['ampOffsets'] = _dictOfArraysToDictOfLists(self.ampOffsets)
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
        inDict['inputExpPairMjdStartList'] = dict()
        inDict['expIdMask'] = dict()
        inDict['expIdRolloffMask'] = dict()
        inDict['rawExpTimes'] = dict()
        inDict['rawMeans'] = dict()
        inDict['rawVars'] = dict()
        inDict['rawDeltas'] = dict()
        inDict['rowMeanVariance'] = dict()
        inDict['gain'] = dict()
        inDict['gainErr'] = dict()
        inDict['gainUnadjusted'] = dict()
        inDict['gainList'] = dict()
        inDict['noiseList'] = dict()
        inDict['overscanMedianLevelList'] = dict()
        inDict['overscanMedian'] = dict()
        inDict['overscanMedianSigma'] = dict()
        inDict['noise'] = dict()
        inDict['noiseErr'] = dict()
        inDict['histVars'] = dict()
        inDict['histChi2Dofs'] = dict()
        inDict['kspValues'] = dict()
        inDict['ptcFitPars'] = dict()
        inDict['ptcFitParsError'] = dict()
        inDict['ptcFitChiSq'] = dict()
        inDict['ptcTurnoff'] = dict()
        inDict['ptcTurnoffSamplingError'] = dict()
        inDict['ptcRolloff'] = dict()
        inDict['ptcRolloffError'] = dict()
        inDict['ptcRolloffTau'] = dict()
        inDict['ptcRolloffTauError'] = dict()
        inDict['nPixelCovariances'] = dict()
        inDict['covariances'] = dict()
        inDict['covariancesModel'] = dict()
        inDict['covariancesSqrtWeights'] = dict()
        inDict['aMatrix'] = dict()
        inDict['bMatrix'] = dict()
        inDict['noiseMatrix'] = dict()
        inDict['finalVars'] = dict()
        inDict['finalModelVars'] = dict()
        inDict['finalMeans'] = dict()
        inDict['badAmps'] = []
        inDict['photoCharges'] = dict()
        inDict['photoChargeDeltas'] = dict()
        inDict['ampOffsets'] = dict()

        # TODO: DM-47610, remove after v29
        inDict['noiseMatrixNoB'] = dict()
        inDict['covariancesModelNoB'] = dict()
        inDict['aMatrixNoB'] = dict()

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
            elif calibVersion >= 1.3 and calibVersion < 2.1:
                inDict['noiseMatrix'][ampName] = record['NOISE_MATRIX']
                inDict['noiseMatrixNoB'][ampName] = record['NOISE_MATRIX_NO_B']
            else:
                inDict['noiseMatrix'][ampName] = record['NOISE_MATRIX']
                nanMatrix = np.full_like(inDict['aMatrix'][ampName], np.nan)
                inDict['noiseMatrixNoB'][ampName] = nanMatrix
            if calibVersion < 1.5:
                # Matched to `COV_MATRIX_SIDE`. Same for all amps.
                inDict['covMatrixSideFullCovFit'] = inDict['covMatrixSide']
            else:
                inDict['covMatrixSideFullCovFit'] = record['COV_MATRIX_SIDE_FULL_COV_FIT']
            if calibVersion < 1.6:
                inDict['rowMeanVariance'][ampName] = np.full((len(inDict['expIdMask'][ampName]),), np.nan)
            else:
                inDict['rowMeanVariance'][ampName] = record['ROW_MEAN_VARIANCE']
            if calibVersion < 1.7:
                inDict['noiseList'][ampName] = np.full_like(inDict['rawMeans'][ampName], np.nan)
            else:
                inDict['noiseList'][ampName] = record['NOISE_LIST']
            if calibVersion < 1.8:
                inDict['ptcTurnoffSamplingError'][ampName] = np.nan
            else:
                inDict['ptcTurnoffSamplingError'][ampName] = record['PTC_TURNOFF_SAMPLING_ERROR']
            if calibVersion < 1.9 and inDict['ptcFitType'] == "FULLCOVARIANCE":
                # Before version 1.9, the noise stored in the PTC was in
                # units of electron^2 only if ptcFitType == FULLCOVARIANCE.
                # After version 1.9, we standardized the
                # PhotonTransferCurveDataset noise units to electron to fix
                # this bug. If a user tries to use an earlier version of
                # PTC with this fit type, we must be sure to do the
                # calculations properly. More information about this noise
                # issue can be found in DM-45976.
                if ampName == inDict['ampNames'][0]:
                    cls().log.info(f"Input PTC VERSION ({calibVersion}) < 1.9 and"
                                   " ptcFitType == FULLCOVARIANCE. Applying fix for"
                                   f" the DM-45976 noise issue.")
                # The noiseErr calculation was accidentally correct in the
                # previous version, so we only need to upday the noise
                # attribute.
                inDict['noise'][ampName] = np.sqrt(record['noise'][ampName])
            if calibVersion < 2.0:
                inDict['ampOffsets'][ampName] = np.full_like(inDict['rawMeans'][ampName], np.nan)
                inDict['gainUnadjusted'][ampName] = record['GAIN']
                inDict['gainList'][ampName] = np.full_like(inDict['rawMeans'][ampName], np.nan)
            else:
                inDict['ampOffsets'][ampName] = record['AMP_OFFSETS']
                inDict['gainUnadjusted'][ampName] = record['GAIN_UNADJUSTED']
                inDict['gainList'][ampName] = record['GAIN_LIST']
            if calibVersion < 2.1:
                inDict['covariancesModelNoB'][ampName] = record['COVARIANCES_MODEL_NO_B']
                inDict['aMatrixNoB'][ampName] = record['A_MATRIX_NO_B']
            else:
                nanMatrixList = np.full_like(inDict['covariances'][ampName], np.nan)
                inDict['covariancesModelNoB'][ampName] = nanMatrixList
                nanMatrix = np.full_like(inDict['aMatrix'][ampName], np.nan)
                inDict['aMatrixNoB'][ampName] = nanMatrix
            if calibVersion < 2.2:
                inDict['inputExpPairMjdStartList'][ampName] = np.full_like(
                    inDict['rawMeans'][ampName],
                    np.nan,
                )
                inDict['overscanMedianLevelList'][ampName] = np.full_like(
                    inDict['rawMeans'][ampName],
                    np.nan,
                )
            else:
                inDict['inputExpPairMjdStartList'][ampName] = record['INPUT_EXP_PAIR_MJD_START']
                inDict['overscanMedianLevelList'][ampName] = record['OVERSCAN_MEDIAN_LIST']
            if calibVersion < 2.3:
                inDict['overscanMedian'][ampName] = np.nan
                inDict['overscanMedianSigma'][ampName] = np.nan
            else:
                inDict['overscanMedian'][ampName] = record['OVERSCAN_MEDIAN']
                inDict['overscanMedianSigma'][ampName] = record['OVERSCAN_MEDIAN_SIGMA']
            if calibVersion < 2.4:
                inDict['nPixelCovariances'][ampName] = -1
            else:
                inDict['nPixelCovariances'][ampName] = record['NPIXEL_COVARIANCES']
            if calibVersion < 2.5:
                inDict['rawDeltas'][ampName] = np.full_like(
                    inDict['rawMeans'][ampName],
                    np.nan,
                )
                inDict['photoChargeDeltas'][ampName] = np.full_like(
                    inDict['rawMeans'][ampName],
                    np.nan,
                )
            else:
                inDict['rawDeltas'][ampName] = record['RAW_DELTAS']
                inDict['photoChargeDeltas'][ampName] = record['PHOTO_CHARGE_DELTAS']
            if calibVersion < 2.6:
                inDict['ptcRolloff'][ampName] = np.nan
                inDict['ptcRolloffError'][ampName] = np.nan
                inDict['ptcRolloffTau'][ampName] = np.nan
                inDict['ptcRolloffTauError'][ampName] = np.nan
                inDict['expIdRolloffMask'][ampName] = np.full_like(
                    inDict['expIdMask'][ampName],
                    False,
                )
            else:
                inDict['ptcRolloff'][ampName] = record['PTC_ROLLOFF']
                inDict['ptcRolloffError'][ampName] = record['PTC_ROLLOFF_ERROR']
                inDict['ptcRolloffTau'][ampName] = record['PTC_ROLLOFF_TAU']
                inDict['ptcRolloffTauError'][ampName] = record['PTC_ROLLOFF_TAU_ERROR']
                inDict['expIdRolloffMask'][ampName] = record['EXP_ID_ROLLOFF_MASK']

        inDict['auxValues'] = {}
        record = ptcTable[0]
        for col in record.columns.keys():
            if col.startswith('PTCAUX_'):
                parts = col.split('PTCAUX_')
                if isinstance(record[col][0], np.bytes_):
                    # Convert to a unicode string because astropy fits doesn't
                    # round-trip properly
                    inDict['auxValues'][parts[1]] = record[col].astype(np.str_)
                else:
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
                'INPUT_EXP_PAIR_MJD_START': self.inputExpPairMjdStartList[ampName],
                'EXP_ID_MASK': self.expIdMask[ampName],
                'EXP_ID_ROLLOFF_MASK': self.expIdRolloffMask[ampName],
                'RAW_EXP_TIMES': self.rawExpTimes[ampName],
                'RAW_MEANS': self.rawMeans[ampName],
                'RAW_VARS': self.rawVars[ampName],
                'RAW_DELTAS': self.rawDeltas[ampName],
                'ROW_MEAN_VARIANCE': self.rowMeanVariance[ampName],
                'GAIN': self.gain[ampName],
                'GAIN_ERR': self.gainErr[ampName],
                'GAIN_UNADJUSTED': self.gainUnadjusted[ampName],
                'GAIN_LIST': self.gainList[ampName],
                'OVERSCAN_MEDIAN_LIST': self.overscanMedianLevelList[ampName],
                'OVERSCAN_MEDIAN': self.overscanMedian[ampName],
                'OVERSCAN_MEDIAN_SIGMA': self.overscanMedianSigma[ampName],
                'NOISE_LIST': self.noiseList[ampName],
                'NOISE': self.noise[ampName],
                'NOISE_ERR': self.noiseErr[ampName],
                'HIST_VARS': self.histVars[ampName],
                'HIST_CHI2_DOFS': self.histChi2Dofs[ampName],
                'KS_PVALUES': self.kspValues[ampName],
                'PTC_FIT_PARS': np.array(self.ptcFitPars[ampName]),
                'PTC_FIT_PARS_ERROR': np.array(self.ptcFitParsError[ampName]),
                'PTC_FIT_CHI_SQ': self.ptcFitChiSq[ampName],
                'PTC_TURNOFF': self.ptcTurnoff[ampName],
                'PTC_TURNOFF_SAMPLING_ERROR': self.ptcTurnoffSamplingError[ampName],
                'PTC_ROLLOFF': self.ptcRolloff[ampName],
                'PTC_ROLLOFF_ERROR': self.ptcRolloffError[ampName],
                'PTC_ROLLOFF_TAU': self.ptcRolloffTau[ampName],
                'PTC_ROLLOFF_TAU_ERROR': self.ptcRolloffTauError[ampName],
                'A_MATRIX': self.aMatrix[ampName].ravel(),
                'B_MATRIX': self.bMatrix[ampName].ravel(),
                'NOISE_MATRIX': self.noiseMatrix[ampName].ravel(),
                'BAD_AMPS': badAmps,
                'PHOTO_CHARGE': self.photoCharges[ampName],
                'PHOTO_CHARGE_DELTAS': self.photoChargeDeltas[ampName],
                'AMP_OFFSETS': self.ampOffsets[ampName],
                'NPIXEL_COVARIANCES': self.nPixelCovariances[ampName],
                'COVARIANCES': self.covariances[ampName].ravel(),
                'COVARIANCES_MODEL': self.covariancesModel[ampName].ravel(),
                'COVARIANCES_SQRT_WEIGHTS': self.covariancesSqrtWeights[ampName].ravel(),
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

    def appendPartialPtc(self, partialPtc):
        """Append a partial PTC dataset to this dataset.

        Parameters
        ----------
        partialPtc : `lsst.ip.isr.PhotonTransferCurveDataset`
            Partial PTC to append. Should only have one element.
        """
        if self.ampNames != partialPtc.ampNames:
            raise ValueError("partialPtc has mis-matched amps.")
        if len(partialPtc.rawMeans[self.ampNames[0]]) != 1 or partialPtc.ptcFitType != "PARTIAL":
            raise ValueError("partialPtc does not appear to be the correct format.")

        # Record the initial length of the PTC, for checking auxValues.
        initialLength = len(self.expIdMask[self.ampNames[0]])

        for key, value in partialPtc.auxValues.items():
            if key in self.auxValues:
                self.auxValues[key] = np.append(self.auxValues[key], value)
            elif initialLength == 0:
                # This is the first partial, so we can set the dict key.
                self.auxValues[key] = value
            else:
                raise ValueError(f"partialPtc has mismatched auxValue key {key}.")

        for ampName in self.ampNames:
            if initialLength == 0:
                # This is the first partial, so we can set the dict key.
                self.nPixelCovariances[ampName] = partialPtc.nPixelCovariances[ampName]
            elif partialPtc.nPixelCovariances[ampName] != self.nPixelCovariances[ampName]:
                raise ValueError(f"partialPtc has mismatched nPixelCovariances for amp {ampName}.")

        for ampName in self.ampNames:
            # The partial dataset consists of lists of values for each
            # quantity. In the case of the input exposure pairs and the
            # input exposure MJDs, this is a list of tuples. In all cases
            # we only want the first (and only) element of the list.
            self.inputExpIdPairs[ampName].append(partialPtc.inputExpIdPairs[ampName][0])
            self.inputExpPairMjdStartList[ampName] = np.append(
                self.inputExpPairMjdStartList[ampName],
                partialPtc.inputExpPairMjdStartList[ampName][0],
            )
            self.expIdMask[ampName] = np.append(self.expIdMask[ampName],
                                                partialPtc.expIdMask[ampName][0])
            self.expIdRolloffMask[ampName] = np.append(self.expIdRolloffMask[ampName],
                                                       partialPtc.expIdRolloffMask[ampName][0])
            self.rawExpTimes[ampName] = np.append(self.rawExpTimes[ampName],
                                                  partialPtc.rawExpTimes[ampName][0])
            self.rawMeans[ampName] = np.append(self.rawMeans[ampName],
                                               partialPtc.rawMeans[ampName][0])
            self.rawVars[ampName] = np.append(self.rawVars[ampName],
                                              partialPtc.rawVars[ampName][0])
            self.rawDeltas[ampName] = np.append(self.rawDeltas[ampName],
                                                partialPtc.rawDeltas[ampName][0])
            self.rowMeanVariance[ampName] = np.append(self.rowMeanVariance[ampName],
                                                      partialPtc.rowMeanVariance[ampName][0])
            self.photoCharges[ampName] = np.append(self.photoCharges[ampName],
                                                   partialPtc.photoCharges[ampName][0])
            self.photoChargeDeltas[ampName] = np.append(self.photoChargeDeltas[ampName],
                                                        partialPtc.photoChargeDeltas[ampName][0])
            self.ampOffsets[ampName] = np.append(self.ampOffsets[ampName],
                                                 partialPtc.ampOffsets[ampName][0])
            self.histVars[ampName] = np.append(self.histVars[ampName],
                                               partialPtc.histVars[ampName][0])
            self.histChi2Dofs[ampName] = np.append(self.histChi2Dofs[ampName],
                                                   partialPtc.histChi2Dofs[ampName][0])
            self.kspValues[ampName] = np.append(self.kspValues[ampName],
                                                partialPtc.kspValues[ampName][0])
            self.gainList[ampName] = np.append(self.gainList[ampName],
                                               partialPtc.gain[ampName])
            self.overscanMedianLevelList[ampName] = np.append(
                self.overscanMedianLevelList[ampName],
                partialPtc.overscanMedianLevelList[ampName][0],
            )
            self.noiseList[ampName] = np.append(self.noiseList[ampName],
                                                partialPtc.noise[ampName])
            self.finalVars[ampName] = np.append(self.finalVars[ampName],
                                                partialPtc.finalVars[ampName][0])
            self.finalModelVars[ampName] = np.append(self.finalModelVars[ampName],
                                                     partialPtc.finalModelVars[ampName][0])
            self.finalMeans[ampName] = np.append(self.finalMeans[ampName],
                                                 partialPtc.finalMeans[ampName][0])
            self.covariances[ampName] = np.append(
                self.covariances[ampName].ravel(),
                partialPtc.covariances[ampName].ravel()
            ).reshape(
                (
                    len(self.rawExpTimes[ampName]),
                    self.covMatrixSide,
                    self.covMatrixSide,
                )
            )
            self.covariancesSqrtWeights[ampName] = np.append(
                self.covariancesSqrtWeights[ampName].ravel(),
                partialPtc.covariancesSqrtWeights[ampName].ravel()
            ).reshape(
                (
                    len(self.rawExpTimes[ampName]),
                    self.covMatrixSide,
                    self.covMatrixSide,
                )
            )
            self.covariancesModel[ampName] = np.append(
                self.covariancesModel[ampName].ravel(),
                partialPtc.covariancesModel[ampName].ravel()
            ).reshape(
                (
                    len(self.rawExpTimes[ampName]),
                    self.covMatrixSide,
                    self.covMatrixSide,
                )
            )

    def sort(self, sortIndex):
        """Sort the components of the PTC by a given sort index.

        The PTC is sorted in-place.

        Parameters
        ----------
        sortIndex : `list` or `np.ndarray`
            The sorting index, which must be the same length as
            the number of elements of the PTC.
        """
        index = np.atleast_1d(sortIndex)

        # First confirm everything matches.
        for ampName in self.ampNames:
            if len(index) != len(self.rawExpTimes[ampName]):
                raise ValueError(
                    f"Length of sortIndex ({len(index)}) does not match number of PTC "
                    f"elements ({len(self.rawExpTimes[ampName])})",
                )

        # Note that gain, gainUnadjusted, gainErr, noise, noiseErr,
        # ptcTurnoff, ptcTurnoffSamplingError, and the full covariance fit
        # parameters are global and not sorted by input pair.

        for ampName in self.ampNames:
            self.inputExpIdPairs[ampName] = np.array(
                self.inputExpIdPairs[ampName]
            )[index].tolist()
            self.inputExpPairMjdStartList[ampName] = self.inputExpPairMjdStartList[ampName][index]

            self.expIdMask[ampName] = self.expIdMask[ampName][index]
            self.expIdRolloffMask[ampName] = self.expIdRolloffMask[ampName][index]
            self.rawExpTimes[ampName] = self.rawExpTimes[ampName][index]
            self.rawMeans[ampName] = self.rawMeans[ampName][index]
            self.rawVars[ampName] = self.rawVars[ampName][index]
            self.rawDeltas[ampName] = self.rawDeltas[ampName][index]
            self.rowMeanVariance[ampName] = self.rowMeanVariance[ampName][index]
            self.photoCharges[ampName] = self.photoCharges[ampName][index]
            self.photoChargeDeltas[ampName] = self.photoChargeDeltas[ampName][index]
            self.ampOffsets[ampName] = self.ampOffsets[ampName][index]

            self.gainList[ampName] = self.gainList[ampName][index]
            self.noiseList[ampName] = self.noiseList[ampName][index]

            self.overscanMedianLevelList[ampName] = self.overscanMedianLevelList[ampName][index]

            self.histVars[ampName] = self.histVars[ampName][index]
            self.histChi2Dofs[ampName] = self.histChi2Dofs[ampName][index]
            self.kspValues[ampName] = self.kspValues[ampName][index]

            self.covariances[ampName] = self.covariances[ampName][index]
            self.covariancesSqrtWeights[ampName] = self.covariancesSqrtWeights[ampName][index]
            self.covariancesModel[ampName] = self.covariancesModel[ampName][index]

            self.finalVars[ampName] = self.finalVars[ampName][index]
            self.finalModelVars[ampName] = self.finalModelVars[ampName][index]
            self.finalMeans[ampName] = self.finalMeans[ampName][index]

        # Sort the auxiliary values which are not stored per-amp.
        for key, value in self.auxValues.items():
            self.auxValues[key] = value[index]

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

    def evalPtcModel(self, mu):
        """Computes the covariance model at specific signal levels.

        Parameters
        ----------
        mu : `numpy.array`, (N,)
            List of mean signals in ADU.

        Raises
        ------
        RuntimeError
            Raised if ptcFitType is invalid.

        Returns
        -------
        covModel : `numpy.array`, (N, M, M)
            Covariances model at mu (in ADU^2).

        Notes
        -----
        Computes the covModel for all mu, and it returns
        cov[N, M, M], where the variance model is cov[:,0,0].
        Both mu and cov are in ADUs and ADUs squared. This
        routine evaulates the n-degree polynomial model (defined
        by polynomialFitDegree) if self.ptcFitType == POLYNOMIAL,
        the approximation in Eq. 16 of Astier+19 (1905.08677)
        if self.ptcFitType == EXPAPPROXIMATION, and Eq. 20 of
        Astier+19 if self.ptcFitType == FULLCOVARIANCE.

        The POLYNOMIAL model and the EXPAPPROXIMATION model
        (Eq. 16 of Astier+19) are only approximations for the
        variance (cov[0,0]), so the function returns covModel
        of shape (N,), representing an array of [C_{00}]
        if self.ptcFitType == EXPAPPROXIMATION or
        self.ptcFitType == POLYNOMAIL.
        """

        ampNames = self.ampNames
        covModel = {ampName: np.array([]) for ampName in ampNames}

        if self.ptcFitType == "POLYNOMIAL":
            pars = self.ptcFitPars

            for ampName in ampNames:
                c00 = poly.polyval(mu, [*pars[ampName]])
                covModel[ampName] = c00

        elif self.ptcFitType == "EXPAPPROXIMATION":
            pars = self.ptcFitPars

            for ampName in ampNames:
                a00, gain, noise = pars[ampName]
                f1 = 0.5/(a00*gain*gain)*(np.exp(2*a00*mu*gain)-1)
                f2 = noise/(gain*gain)
                c00 = f1 + f2
                covModel[ampName] = c00

        elif self.ptcFitType in ["FULLCOVARIANCE", "FULLCOVARIANCE_NO_B"]:
            for ampName in ampNames:
                noiseMatrix = self.noiseMatrix[ampName]
                gain = self.gain[ampName]
                aMatrix = self.aMatrix[ampName]
                bMatrix = self.bMatrix[ampName]
                cMatrix = aMatrix*bMatrix

                matrixSideFit = self.covMatrixSideFullCovFit
                sa = (matrixSideFit, matrixSideFit)

                # pad a with zeros and symmetrize
                aEnlarged = np.zeros((int(sa[0]*1.5)+1, int(sa[1]*1.5)+1))
                aEnlarged[0:sa[0], 0:sa[1]] = aMatrix
                aSym = symmetrize(aEnlarged)

                # pad c with zeros and symmetrize
                cEnlarged = np.zeros((int(sa[0]*1.5)+1, int(sa[1]*1.5)+1))
                cEnlarged[0:sa[0], 0:sa[1]] = cMatrix

                cSym = symmetrize(cEnlarged)
                a2 = fftconvolve(aSym, aSym, mode='same')
                a3 = fftconvolve(a2, aSym, mode='same')
                ac = fftconvolve(aSym, cSym, mode='same')
                (xc, yc) = np.unravel_index(np.abs(aSym).argmax(), a2.shape)

                a1 = aMatrix[np.newaxis, :, :]
                a2 = a2[np.newaxis, xc:xc + matrixSideFit, yc:yc + matrixSideFit]
                a3 = a3[np.newaxis, xc:xc + matrixSideFit, yc:yc + matrixSideFit]
                ac = ac[np.newaxis, xc:xc + matrixSideFit, yc:yc + matrixSideFit]
                c1 = cMatrix[np.newaxis, ::]

                # assumes that mu is 1d
                bigMu = mu[:, np.newaxis, np.newaxis]*gain
                # c(=a*b in Astier+19) also has a contribution to the last
                # term, that is absent for now.

                covModel[ampName] = (bigMu/(gain*gain)*(a1*bigMu+2./3.*(bigMu*bigMu)*(a2 + c1)
                                     + (1./3.*a3 + 5./6.*ac)*(bigMu*bigMu*bigMu))
                                     + noiseMatrix[np.newaxis, :, :]/gain**2)

                # add the Poisson term, and the read out noise (variance)
                covModel[ampName][:, 0, 0] += mu/gain
        else:
            raise RuntimeError("Cannot compute PTC  model for "
                               "ptcFitType %s." % self.ptcFitType)

        return covModel

    def _validateCovarianceMatrizSizes(self):
        """Ensure  covMatrixSideFullCovFit <= covMatrixSide."""
        if self.covMatrixSideFullCovFit > self.covMatrixSide:
            self.log.warning("covMatrixSideFullCovFit > covMatrixSide "
                             f"({self.covMatrixSideFullCovFit} > {self.covMatrixSide})."
                             "Setting the former to the latter.")
            self.covMatrixSideFullCovFit = self.covMatrixSide
