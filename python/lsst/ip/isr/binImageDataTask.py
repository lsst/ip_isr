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

__all__ = ["BinImageDataTask", "BinImageDataConfig", "binImageData"]

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
from lsst.utils.timer import timeMethod


class BinImageDataConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "exposure", "detector"),
    defaultTemplates={"inputName": "postISRCCD", "outputName": "postISRCCDBin"},
):

    inputData = cT.Input(
        name="{inputName}",
        doc="Input image data to bin.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputData = cT.Output(
        name="{outputName}",
        doc="Binned image data.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )

    def __init__(self, *, config=None):
        """Customize the connections and storageClass for a specific
        instance. This enables both to be dynamically set at runtime,
        allowing BinImageDataTask to work with different types of
        image and image-like data.

        Parameters
        ----------
        config : `BinExposureConfig`
            A config for `BinExposureTask`.
        """
        super().__init__(config=config)
        if config and config.inputDimensions != self.inputData.dimensions:
            self.dimensions.clear()
            self.dimensions.update(config.inputDimensions)
            self.inputData = cT.Input(
                name=self.inputData.name,
                doc=self.inputData.doc,
                storageClass=self.inputData.storageClass,
                dimensions=frozenset(config.inputDimensions),
            )
            self.outputData = cT.Output(
                name=self.outputData.name,
                doc=self.outputData.doc,
                storageClass=self.outputData.storageClass,
                dimensions=frozenset(config.inputDimensions),
            )
        if config and config.inputStorageClass != self.inputData.storageClass:
            self.inputData = cT.Input(
                name=self.inputData.name,
                doc=self.inputData.doc,
                storageClass=config.inputStorageClass,
                dimensions=self.inputData.dimensions,
            )
            self.outputData = cT.Output(
                name=self.outputData.name,
                doc=self.outputData.doc,
                storageClass=config.inputStorageClass,
                dimensions=self.outputData.dimensions,
            )


class BinImageDataConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=BinImageDataConnections
):
    """Config for BinImageDataTask"""

    inputDimensions = pexConfig.ListField(
        # Sort to ensure default order is consistent between runs
        default=sorted(BinImageDataConnections.dimensions),
        dtype=str,
        doc="Override for the dimensions of the input and output data.",
    )
    inputStorageClass = pexConfig.Field(
        default="ExposureF",
        dtype=str,
        doc=(
            "Override the storageClass of the input and output data. "
            "Must be of type `Image`, `MaskedImage`, or `Exposure`, "
            "or one of their subtypes."
        ),
    )
    binFactor = pexConfig.Field(
        dtype=int,
        doc="Binning factor applied to both spatial dimensions.",
        default=8,
        check=lambda x: x > 1,
    )


class BinImageDataTask(pipeBase.PipelineTask):
    """Perform an nxn binning of an image or image-like dataset.

    The binning factor is the same in both spatial dimensions (i.e.,
    an nxn binning is performed). In the case of MaskedImages and Exposures,
    each of the input image planes are binned by the same factor.
    """

    ConfigClass = BinImageDataConfig
    _DefaultName = "binImageData"

    @timeMethod
    def run(self, inputData, binFactor=None):
        """Perform an nxn binning of image and image-like data.

        Parameters:
        -----------
        inputData : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage` or
                    `lsst.afw.image.Exposure` or one of their sub-types.
            Data to spatially bin
        binFactor : `int`, optional.
            nxn binning factor. If not provided then self.config.binFactor
            is used.

        Returns:
        --------
        result : `lsst.pipe.base.Struct`
            Results as a struct with attributes:

            ``outputData``
               Binned data (`lsst.afw.image.Image` or
               `lsst.afw.image.MaskedImage` or `lsst.afw.image.Exposure`
               or one of their sub-types. The type matches that of the input.).
        """
        if not binFactor:
            binFactor = self.config.binFactor
        return pipeBase.Struct(outputData=binImageData(inputData, binFactor))


def binImageData(inputData, binFactor=8):
    """Bin image and image-like data to reduce its spatial dimensions.

    Performs an nxn binning of the input data, reducing both spatial
    dimensions of each of the input image data by the provided
    factor.

    Parameters:
    -----------
    inputData:  `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage` or
                `lsst.afw.image.Exposure` or one of their sub-types.
        Input data to bin.
    binFactor:  `int`
        Binning factor to apply to each input exposure's image data.
        Default 8.

    Returns:
    --------
    binnedImage or binnedExposure: `lsst.afw.image.Image` or
        `lsst.afw.image.MaskedImage` or `lsst.afw.image.Exposure` or one of
        their sub-types.
        Binned version of input image.

    Raises
    ------
    TypeError
        Raised if either the binning factor is not of type `int`, or if the
        input data to be binned is not of type `lsst.afw.image.Exposure`
        or one of its sub-types.
    """

    if not isinstance(binFactor, int):
        raise TypeError("binFactor must be of type int")

    if isinstance(inputData, afwImage.Exposure):
        inputImage = inputData.getMaskedImage()
        isExposure = True
    elif isinstance(inputData, (afwImage.Image, afwImage.MaskedImage)):
        inputImage = inputData
        isExposure = False
    else:
        message = (
            "inputData must be of type `lsst.afw.image.Image`, `lsst.afw.MaskedImage`, "
            "or `lsst.afw.image.Exposure`, or one of their sub-types."
        )
        raise TypeError(message)

    binnedImage = afwMath.binImage(inputImage, binFactor)

    if isExposure:
        binnedExposure = afwImage.makeExposure(binnedImage)
        binnedExposure.setInfo(inputData.getInfo())
        return binnedExposure
    else:
        return binnedImage
