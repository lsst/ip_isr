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
#

__all__ = ["IsrQaFlatnessConfig", "IsrQaConfig", "makeThumbnail"]

import lsst.afw.display.rgb as afwRGB
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig


class IsrQaFlatnessConfig(pexConfig.Config):
    meshX = pexConfig.Field(
        dtype=int,
        doc="Mesh size in X for flatness statistics",
        default=256,
    )
    meshY = pexConfig.Field(
        dtype=int,
        doc="Mesh size in Y for flatness statistics",
        default=256,
    )
    doClip = pexConfig.Field(
        dtype=bool,
        doc="Clip outliers for flatness statistics?",
        default=True,
    )
    clipSigma = pexConfig.Field(
        dtype=float,
        doc="Number of sigma deviant a pixel must be to be clipped from flatness statistics.",
        default=3.0,
    )
    nIter = pexConfig.Field(
        dtype=int,
        doc="Number of iterations used for outlier clipping in flatness statistics.",
        default=3,
    )


class IsrQaConfig(pexConfig.Config):
    saveStats = pexConfig.Field(
        dtype=bool,
        doc="Calculate ISR statistics while processing?",
        default=True,
    )

    flatness = pexConfig.ConfigField(
        dtype=IsrQaFlatnessConfig,
        doc="Flatness statistics configuration.",
    )

    doWriteOss = pexConfig.Field(
        dtype=bool,
        doc="Write overscan subtracted image?",
        default=False,
    )
    doThumbnailOss = pexConfig.Field(
        dtype=bool,
        doc="Write overscan subtracted thumbnail?",
        default=False,
    )

    doWriteFlattened = pexConfig.Field(
        dtype=bool,
        doc="Write image after flat-field correction?",
        default=False,
    )
    doThumbnailFlattened = pexConfig.Field(
        dtype=bool,
        doc="Write thumbnail after flat-field correction?",
        default=False,
    )

    thumbnailBinning = pexConfig.Field(
        dtype=int,
        doc="Thumbnail binning factor.",
        default=4,
    )
    thumbnailStdev = pexConfig.Field(
        dtype=float,
        doc="Number of sigma below the background to set the thumbnail minimum.",
        default=3.0,
    )
    thumbnailRange = pexConfig.Field(
        dtype=float,
        doc="Total range in sigma for thumbnail mapping.",
        default=5.0,
    )
    thumbnailQ = pexConfig.Field(
        dtype=float,
        doc="Softening parameter for thumbnail mapping.",
        default=20.0,
    )
    thumbnailSatBorder = pexConfig.Field(
        dtype=int,
        doc="Width of border around saturated pixels in thumbnail.",
        default=2,
    )


def makeThumbnail(exposure, isrQaConfig=None):
    """Create a snapshot thumbnail from input exposure.

    The output thumbnail image is constructed based on the parameters
    in the configuration file.  Currently, the asinh mapping is the
    only mapping method used.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The exposure to be converted into a thumbnail.
    isrQaConfig : `Config`, optional
        Configuration object containing all parameters to control the
        thumbnail generation.

    Returns
    -------
    rgbImage : `numpy.ndarray`
        Binned and scaled version of the exposure, converted to an
        integer array to allow it to be written as PNG.
    """
    if isrQaConfig is not None:
        binning = isrQaConfig.thumbnailBinning
        binnedImage = afwMath.binImage(exposure.getMaskedImage(), binning, binning, afwMath.MEAN)

        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setAndMask(binnedImage.getMask().getPlaneBitMask(["SAT", "BAD", "INTRP"]))
        stats = afwMath.makeStatistics(binnedImage,
                                       afwMath.MEDIAN | afwMath.STDEVCLIP | afwMath.MAX, statsCtrl)

        low = stats.getValue(afwMath.MEDIAN) - isrQaConfig.thumbnailStdev*stats.getValue(afwMath.STDEVCLIP)

        if isrQaConfig.thumbnailSatBorder:
            afwRGB.replaceSaturatedPixels(binnedImage, binnedImage, binnedImage,
                                          isrQaConfig.thumbnailSatBorder, stats.getValue(afwMath.MAX))

        asinhMap = afwRGB.AsinhMapping(low, isrQaConfig.thumbnailRange, Q=isrQaConfig.thumbnailQ)
        rgbImage = asinhMap.makeRgbImage(binnedImage)

        return rgbImage
