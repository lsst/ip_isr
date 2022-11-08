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

__all__ = ["FreqDomainMetricsConfig", "FreqDomainMetricsTask"]

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import Exposure
import numpy as np
from scipy.signal.windows import hamming, hann, gaussian
from typing import Callable, Union, Tuple
from itertools import permutations


class FreqDomainMetricsConfig(pexConfig.Config):
    """Options for the calculation of frequency domain metrics of post-ISR images"""

    windowRadial = pexConfig.Field(
        doc="construct radial type window rather than separable",
        dtype=bool,
        default=False,
    )
    windowType = pexConfig.ChoiceField(
        doc="type of window function to use",
        dtype=str,
        default="HAMMING",
        allowed={
            "HAMMING": "A Hamming type window",
            "HANN": "a Hann type window",
            "GAUSSIAN": "a Gaussian window",
        },
    )
    transformDimsType = pexConfig.ChoiceField(
        doc="the type of Fourier Transform data to calculate",
        dtype=str,
        default="1DSLICE",
        allowed={
            "1DSLICE": "calculate the full 2D transform but retain only 1D axis-wise slices through it",
#Note: for now disallowing 2D until we have a more sensible serialization
#            "2D": "calculate and retain the full 2D transform",
        },
    )

    doMeanSubtract = pexConfig.Field(
        doc="subtract the mean before transform (to get canonical DC component value",
        dtype=bool,
        default=True,
    )

    useProjSliceTheorem = pexConfig.Field(
        doc="use the projection slice theorem to calculate sliced transform",
        dtype=bool,
        default=False,
    )

    doCommonModePSD = pexConfig.Field(
        doc="calculate the common mode PSD of all amplifier channels in addition to individual amps",
        dtype=bool,
        default=True)

    saveRawFTs = pexConfig.Field(
        doc="save the raw FT calculations (i.e. as well as the common mode PSD or other metrics",
        dtype=bool,
        default=True)
    

class FreqDomainMetricsTask(pipeBase.Task):
    """Task which calculates frequency domain metrics in post ISR images.

    At present it simply calculates configurable Fourier Transforms of images, retaining both phase and amplitude information. This is useful for further post-processing, for example to extrace power spectral densities for noise analysis purposes, or to obtain correlative maps of noise patterns between different amplifier channels

    Parameters
    ----------



    """
    _DefaultName = "freqdomainmetrics"
    ConfigClass = FreqDomainMetricsConfig

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, exposure: Exposure):
        """

        Parameters
        ----------
        """

        if exposure.detector is None:
            raise RuntimeError(
                "require a detector definition to calculate per-amplifier frequency domain metrics!"
            )

        raw_ft_data = {}
        for amplifier in exposure.detector:
            box = amplifier.getBBox()
            name: str = amplifier.getName()
            input_arr = exposure.image.subset(box).array
            if self.config.doMeanSubtract:
                input_arr -= np.mean(input_arr)
            raw_ft_data[name] = self._do_fft_calc(input_arr)

        output_data = {}
        if self.config.doCommonModePSD:
            cmpsd = self.calcCommonModePSD(raw_ft_data, exposure.detector)
            output_data["cmPSD"] = cmpsd

        if self.config.saveRawFTs:
            output_data |= raw_ft_data
        return output_data

    def calcCommonModePSD(self, data: dict, detector) -> np.ndarray:
        cmodepsd = None
        for amplifier in detector:
            name: str = amplifier.getName()
            thisdat = data[name]
            if not cmodepsd:
                cmodepsd = thisdat
            else:
                match thisdat:
                    case (*1dslcs,):
                        cmodepsd = tuple(a + b for a,b in zip(cmodepsd, 1dslcs))
                    case np.ndarray():
                        cmodepsd += thisdat

        match cmodepsd:
            case(*1dslcs,):
                return tuple(np.abs(_) for _ in 1dslcs)
            case np.ndarray():
                return np.abs(cmodepsd)
            
                
    
    def _do_fft_calc(self, inp_arr: np.ndarray) -> Union[np.ndarray, Tuple[np.ndarray]]:
        # Note: for now, all possibilities use the 2D window. In future,
        # we may add a combination where a 1D window function would be appropriate
        wdw = self._calc_wdw(inp_arr.shape)
        d = inp_arr * wdw

        if self.config.transformDimsType == "1DSLICE":
            if self.config.useProjSliceTheorem:
                projs = (np.sum(d, axis=_) for _ in range(len(d.shape)))
                return tuple(np.fft.rfft(_) for _ in projs)
            else:
                trans = np.fft.rfft2(d)
                baseslc = [0] + [slice(None, None, None)] * (len(trans.shape) - 1)
                return tuple((trans[_] for _ in permutations(baseslc, 2)))
        else:
            # For now, the only other case is the 2D case.
            return np.fft.rfft2(d)

    def _calc_wdw(self, shp: tuple) -> np.ndarray:
        if self.config.windowRadial:
            raise NotImplementedError("radial windows are not supported yet!")

        # Note that the below should never fail because the value bounds should be
        # enforced by some kind of pexConfig machinery
        funclookup = {"HAMMING": hamming, "HANN": hann, "GAUSSIAN": gaussian}
        wdwfunc: Callable = funclookup[self.config.windowType]

        match shp:
            case (x,):
                return wdwfunc(x)
            case (*dims,):
                return np.outer(*(wdwfunc(_) for _ in dims))
        raise RuntimeError("invalid case!")
