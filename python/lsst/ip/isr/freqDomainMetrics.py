import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import Exposure
import numpy as np
from scipy.signal.windows import hamming, hann, gaussian

class FreqDomainMetricsConfig(pexConfig.Config):
    """Options for the calculation of frequency domain metrics of post-ISR images"""
     windowRadial = pexConfig.Field(
         doc="construct radial type window rather than separable",
         dtype=bool,
         default=False
         )
     windowType = pexConfig.Field(
         doc="type of window function to use",
         dtype=str,
         default="HAMMING",
         allowed={
             "HAMMING" : "A Hamming type window",
             "HANN" : "a Hann type window",
             "GAUSSIAN": "a Gaussian window"}
         )
     transformDimsType = pexConfig.field(
         doc="the type of Fourier Transform data to calculate",
         dtype=str,
         default="1DSLICE",
         allowed={
             "1DSLICE" : "calculate the full 2D transform but retain only 1D axis-wise slices through it",
             "2D" : "calculate and retain the full 2D transform",
             "1DPROJECT" : "calculate the full 2D transform but retain only 1D axis-parallel projections through it"}
         )
     useProjSliceTheorem = pexConfig.field(
         doc="use the projection slice theorem to compute 1D slices",
         dtype=bool,
         default=False)

class FreqDomainMetricsTask(pipeBase.Task):
    """Task which calculates frequency domain metrics in post ISR images.

    At present it simply calculates configurable Fourier Transforms of images, retaining both phase and amplitude information. This is useful for further post-processing, for example to extrace power spectral densities for noise analysis purposes, or to obtain correlative maps of noise patterns between different amplifier channels

    Parameters
    ----------

    
    
    """
    ConfigClass = FreqDomainMetricsConfig

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, exposure: Exposure):
        """

        Parameters
        ----------
        """ 

        det = exposure.getDetector()

        if not det:
            raise RuntimeError("require a detector definition to calculate per-amplifier frequency domain metrics!")


        for amplifier in det:
            box = amplifier.getBBox()
            
