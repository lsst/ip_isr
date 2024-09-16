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

__all__ = ["BinExposureTask", "BinExposureConfig", "binExposure"]

import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
from lsst.utils.timer import timeMethod
import lsst.afw.image as afwImage


class BinExposureConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "exposure", "detector"),
    defaultTemplates={"inputName": "postISRCCD", "outputName": "postISRCCDBin"}
):

    inputExposure = cT.Input(
        name="{inputName}",
        doc="Input exposure to bin.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )
    binnedExposure = cT.Output(
        name="{outputName}",
        doc="Binned exsposure.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )

    def __init__(self, *, config=None):
        """Customize the connections for a specific instance.
        Parameters
        ----------
        config : `BinExposureConfig`
            A config for `BinExposureTask` or one of its subclasses.
        """
        super().__init__(config=config)
        if config and config.exposureDimensions != self.inputExposure.dimensions:
            self.dimensions.clear()
            self.dimensions.update(config.exposureDimensions)
            self.inputExposure = cT.Input(
                name=self.inputExposure.name,
                doc=self.inputExposure.doc,
                storageClass=self.inputExposure.storageClass,
                dimensions=frozenset(config.exposureDimensions),
            )
            self.binnedExposure = cT.Output(
                name=self.binnedExposure.name,
                doc=self.binnedExposure.doc,
                storageClass=self.binnedExposure.storageClass,
                dimensions=frozenset(config.exposureDimensions),
            )
        if config and config.exposureStorageClass != self.inputExposure.storageClass:
            self.inputExposure = cT.Input(
                name=self.inputExposure.name,
                doc=self.inputExposure.doc,
                storageClass=config.exposureStorageClass,
                dimensions=self.inputExposure.dimensions,
            )
            self.binnedExposure = cT.Output(
                name=self.binnedExposure.name,
                doc=self.binnedExposure.doc,
                storageClass=config.exposureStorageClass,
                dimensions=self.binnedExposure.dimensions,
            )


class BinExposureConfig(
    pipeBase.PipelineTaskConfig,
    pipelineConnections=BinExposureConnections
):
    exposureDimensions = pexConfig.ListField(
        # Sort to ensure default order is consistent between runs
        default=sorted(BinExposureConnections.dimensions),
        dtype=str,
        doc="Override for the dimensions of the input and binned exposures",
    )
    exposureStorageClass = pexConfig.Field(
        default='ExposureF',
        dtype=str,
        doc="Override the storageClass of the input and binned exposures"
    )
    binFactor = pexConfig.Field(
        dtype=int,
        doc="Binning factor for binned exposure.",
        default=8,
        check=lambda x: x > 1,
    )


class BinExposureTask(pipeBase.PipelineTask):
    ConfigClass = BinExposureConfig
    _DefaultName = "binExposure"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self, inputExposure, binFactor=None):
        if not binFactor:
            binFactor = self.config.binFactor
        return pipeBase.Struct(
            binnedExposure=binExposure(inputExposure, binFactor)
        )


def binExposure(inputExp, binFactor=8):
    '''Bin an exposure to reduce its spatial dimensions.

    Takes an input exposure and bins its pixels, reducing the spatial
    dimensions of each of the input exposure's image data by the
    given factor.

    Parameters:
    -----------
    inputExposure: `lsst.afw.image.Exposure`
        Input exposure data to bin.
    binFactor: `int`
        Binning factor to apply to each input exposure's image data.
        Default 8.

    Returns:
    --------
    binnedExp: `lsst.afw.image.Exposure`
        Binned version of input image.
    '''

    if not isinstance(binFactor, int):
        raise TypeError('binFactor must be of type int')
    if not isinstance(inputExp, afwImage.Exposure):
        raise TypeError('inputExp must be of type lsst.afw.image.Exposure')

    binned = inputExp.getMaskedImage()
    binned = afwMath.binImage(binned, binFactor)
    binnedExp = afwImage.makeExposure(binned)

    binnedExp.setInfo(inputExp.getInfo())

    return binnedExp
