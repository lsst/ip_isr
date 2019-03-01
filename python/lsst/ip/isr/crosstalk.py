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
Apply intra-CCD crosstalk corrections
"""

import lsst.afw.math
import lsst.afw.table
import lsst.afw.detection
from lsst.pex.config import Config, Field, ChoiceField
from lsst.pipe.base import Task

__all__ = ["CrosstalkConfig", "CrosstalkTask", "subtractCrosstalk", "writeCrosstalkCoeffs",
           "NullCrosstalkTask"]


class CrosstalkConfig(Config):
    """Configuration for intra-CCD crosstalk removal"""
    minPixelToMask = Field(
        dtype=float,
        doc="Set crosstalk mask plane for pixels over this value.",
        default=45000
    )
    crosstalkMaskPlane = Field(
        dtype=str,
        doc="Name for crosstalk mask plane.",
        default="CROSSTALK"
    )
    crosstalkBackgroundMethod = ChoiceField(
        dtype=str,
        doc="Type of background subtraction to use when applying correction.",
        default="None",
        allowed={
            "None": "Do no background subtraction.",
            "AMP": "Subtract amplifier-by-amplifier background levels.",
            "DETECTOR": "Subtract detector level background."
        },
    )


class CrosstalkTask(Task):
    """Apply intra-CCD crosstalk correction"""
    ConfigClass = CrosstalkConfig
    _DefaultName = 'isrCrosstalk'

    def prepCrosstalk(self, dataRef):
        """Placeholder for crosstalk preparation method, e.g., for inter-CCD crosstalk.

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            Butler reference of the detector data to be processed.

        See also
        --------
        lsst.obs.decam.crosstalk.DecamCrosstalkTask.prepCrosstalk
        """
        return

    def run(self, exposure, crosstalkSources=None, isTrimmed=False):
        """Apply intra-CCD crosstalk correction

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure for which to remove crosstalk.
        crosstalkSources : `defaultdict`, optional
            Image data and crosstalk coefficients from other CCDs/amps that are
            sources of crosstalk in exposure.
            The default for intra-CCD crosstalk here is None.
        isTrimmed : `bool`
            The image is already trimmed.
            This should no longer be needed once DM-15409 is resolved.

        Raises
        ------
        RuntimeError
            Raised if called for a detector that does not have a
            crosstalk correction
        """
        detector = exposure.getDetector()
        if not detector.hasCrosstalk():
            raise RuntimeError("Attempted to correct crosstalk without crosstalk coefficients")
        self.log.info("Applying crosstalk correction")
        subtractCrosstalk(exposure, minPixelToMask=self.config.minPixelToMask,
                          crosstalkStr=self.config.crosstalkMaskPlane, isTrimmed=isTrimmed,
                          backgroundMethod=self.config.crosstalkBackgroundMethod)


# Flips required to get the corner to the lower-left
# (an arbitrary choice; flips are relative, so the choice of reference here is not important)
X_FLIP = {lsst.afw.table.LL: False, lsst.afw.table.LR: True,
          lsst.afw.table.UL: False, lsst.afw.table.UR: True}
Y_FLIP = {lsst.afw.table.LL: False, lsst.afw.table.LR: False,
          lsst.afw.table.UL: True, lsst.afw.table.UR: True}


class NullCrosstalkTask(CrosstalkTask):
    def run(self, exposure, crosstalkSources=None):
        self.log.info("Not performing any crosstalk correction")


def extractAmp(image, amp, corner, isTrimmed=False):
    """Return an image of the amp

    The returned image will have the amp's readout corner in the
    nominated `corner`.

    Parameters
    ----------
    image : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage`
        Image containing the amplifier of interest.
    amp : `lsst.afw.table.AmpInfoRecord`
        Amplifier information.
    corner : `lsst.afw.table.ReadoutCorner` or `None`
        Corner in which to put the amp's readout corner, or `None` for
        no flipping.
    isTrimmed : `bool`
        The image is already trimmed.
        This should no longer be needed once DM-15409 is resolved.

    Returns
    -------
    output : `lsst.afw.image.Image`
        Image of the amplifier in the standard configuration.
    """
    output = image[amp.getBBox() if isTrimmed else amp.getRawDataBBox()]
    ampCorner = amp.getReadoutCorner()
    # Flipping is necessary only if the desired configuration doesn't match what we currently have
    xFlip = X_FLIP[corner] ^ X_FLIP[ampCorner]
    yFlip = Y_FLIP[corner] ^ Y_FLIP[ampCorner]
    return lsst.afw.math.flipImage(output, xFlip, yFlip)


def calculateBackground(mi, badPixels=["BAD"]):
    """Calculate median background in image

    Getting a great background model isn't important for crosstalk correction,
    since the crosstalk is at a low level. The median should be sufficient.

    Parameters
    ----------
    mi : `lsst.afw.image.MaskedImage`
        MaskedImage for which to measure background.
    badPixels : `list` of `str`
        Mask planes to ignore.

    Returns
    -------
    bg : `float`
        Median background level.
    """
    mask = mi.getMask()
    stats = lsst.afw.math.StatisticsControl()
    stats.setAndMask(mask.getPlaneBitMask(badPixels))
    return lsst.afw.math.makeStatistics(mi, lsst.afw.math.MEDIAN, stats).getValue()


def subtractCrosstalk(exposure, badPixels=["BAD"], minPixelToMask=45000,
                      crosstalkStr="CROSSTALK", isTrimmed=False,
                      backgroundMethod="None"):
    """Subtract the intra-CCD crosstalk from an exposure

    We set the mask plane indicated by ``crosstalkStr`` in a target amplifier
    for pixels in a source amplifier that exceed `minPixelToMask`. Note that
    the correction is applied to all pixels in the amplifier, but only those
    that have a substantial crosstalk are masked with ``crosstalkStr``.

    The uncorrected image is used as a template for correction. This is good
    enough if the crosstalk is small (e.g., coefficients < ~ 1e-3), but if it's
    larger you may want to iterate.

    This method needs unittests (DM-18876), but such testing requires
    DM-18610 to allow the test detector to have the crosstalk
    parameters set.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure for which to subtract crosstalk.
    badPixels : `list` of `str`
        Mask planes to ignore.
    minPixelToMask : `float`
        Minimum pixel value (relative to the background level) in
        source amplifier for which to set ``crosstalkStr`` mask plane
        in target amplifier.
    crosstalkStr : `str`
        Mask plane name for pixels greatly modified by crosstalk.
    isTrimmed : `bool`
        The image is already trimmed.
        This should no longer be needed once DM-15409 is resolved.
    backgroundMethod : `str`
        Method used to subtract the background.  "AMP" uses
        amplifier-by-amplifier background levels, "DETECTOR" uses full
        exposure/maskedImage levels.  Any other value results in no
        background subtraction.
    """
    mi = exposure.getMaskedImage()
    mask = mi.getMask()

    ccd = exposure.getDetector()
    numAmps = len(ccd)
    coeffs = ccd.getCrosstalk()
    assert coeffs.shape == (numAmps, numAmps)

    # Set background level based on the requested method.  The
    # thresholdBackground holds the offset needed so that we only mask
    # pixels high relative to the background, not in an absolute
    # sense.
    thresholdBackground = calculateBackground(mi, badPixels)

    backgrounds = [0.0 for amp in ccd]
    if backgroundMethod is None:
        pass
    elif backgroundMethod == "AMP":
        backgrounds = [calculateBackground(mi[amp.getBBox()], badPixels) for amp in ccd]
    elif backgroundMethod == "DETECTOR":
        backgrounds = [calculateBackground(mi, badPixels) for amp in ccd]

    # Set the crosstalkStr bit for the bright pixels (those which will have significant crosstalk correction)
    crosstalkPlane = mask.addMaskPlane(crosstalkStr)
    footprints = lsst.afw.detection.FootprintSet(mi, lsst.afw.detection.Threshold(minPixelToMask +
                                                                                  thresholdBackground))
    footprints.setMask(mask, crosstalkStr)
    crosstalk = mask.getPlaneBitMask(crosstalkStr)

    # Do pixel level crosstalk correction.
    subtrahend = mi.Factory(mi.getBBox())
    subtrahend.set((0, 0, 0))
    for ii, iAmp in enumerate(ccd):
        iImage = subtrahend[iAmp.getBBox() if isTrimmed else iAmp.getRawDataBBox()]
        for jj, jAmp in enumerate(ccd):
            if ii == jj:
                assert coeffs[ii, jj] == 0.0
            if coeffs[ii, jj] == 0.0:
                continue

            jImage = extractAmp(mi, jAmp, iAmp.getReadoutCorner(), isTrimmed)
            jImage.getMask().getArray()[:] &= crosstalk  # Remove all other masks
            jImage -= backgrounds[jj]

            iImage.scaledPlus(coeffs[ii, jj], jImage)

    # Set crosstalkStr bit only for those pixels that have been significantly modified (i.e., those
    # masked as such in 'subtrahend'), not necessarily those that are bright originally.
    mask.clearMaskPlane(crosstalkPlane)
    mi -= subtrahend  # also sets crosstalkStr bit for bright pixels


def writeCrosstalkCoeffs(outputFileName, coeff, det=None, crosstalkName="Unknown", indent=2):
    """Write a yaml file containing the crosstalk coefficients

    The coeff array is indexed by [i, j] where i and j are amplifiers
    corresponding to the amplifiers in det

    Parameters
    ----------
    outputFileName : `str`
        Name of output yaml file
    coeff : `numpy.array(namp, namp)`
        numpy array of coefficients
    det : `lsst.afw.cameraGeom.Detector`
        Used to provide the list of amplifier names;
        if None use ['0', '1', ...]
    ccdType : `str`
        Name of CCD, used to index the yaml file
        If all CCDs are identical could be the type (e.g. ITL)
    indent : `int`
        Indent width to use when writing the yaml file
    """

    if det is None:
        ampNames = [str(i) for i in range(coeff.shape[0])]
    else:
        ampNames = [a.getName() for a in det]

    assert coeff.shape == (len(ampNames), len(ampNames))

    dIndent = indent
    indent = 0
    with open(outputFileName, "w") as fd:
        print(indent*" " + "crosstalk :", file=fd)
        indent += dIndent
        print(indent*" " + "%s :" % crosstalkName, file=fd)
        indent += dIndent

        for i, ampNameI in enumerate(ampNames):
            print(indent*" " + "%s : {" % ampNameI, file=fd)
            indent += dIndent
            print(indent*" ", file=fd, end='')

            for j, ampNameJ in enumerate(ampNames):
                print("%s : %11.4e, " % (ampNameJ, coeff[i, j]), file=fd,
                      end='\n' + indent*" " if j%4 == 3 else '')
            print("}", file=fd)

            indent -= dIndent
