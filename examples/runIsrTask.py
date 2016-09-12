from __future__ import print_function
from lsst.ip.isr import IsrTask
import lsst.afw.display.ds9 as ds9
import exampleUtils
import sys
import numpy


def runIsr():
    '''Run the task to do ISR on a ccd'''

    #Create the isr task with modified config
    isrConfig = IsrTask.ConfigClass()
    isrConfig.doBias = False #We didn't make a zero frame
    isrConfig.doDark = True
    isrConfig.doFlat = True
    isrConfig.doFringe = False #There is no fringe frame for this example

    isrConfig.assembleCcd.setGain = False
    isrTask = IsrTask(config=isrConfig)

    #Make raw, flat and dark exposures
    DARKVAL = 2. #e-/sec
    OSCAN = 1000. #DN
    GRADIENT = .10
    EXPTIME = 15 #seconds
    DARKEXPTIME = 40. #seconds

    darkExposure = exampleUtils.makeDark(DARKVAL, DARKEXPTIME)
    flatExposure = exampleUtils.makeFlat(GRADIENT)
    rawExposure = exampleUtils.makeRaw(DARKVAL, OSCAN, GRADIENT, EXPTIME)

    output = isrTask.run(rawExposure, dark=darkExposure, flat=flatExposure)
    return output.exposure

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of IsrTask")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display the result?", default=False)
    parser.add_argument('--write', '-w', action="store_true", help="Write the result?", default=False)
    args = parser.parse_args()

    if args.debug:
        try:
            import debug
        except ImportError as e:
            print(e, file=sys.stderr)

    exposure = runIsr()

    if args.ds9:
        im = exposure.getMaskedImage().getImage()
        im_median = numpy.median(im.getArray())
        ds9.mtv(im)
        ds9.scale(min=im_median*0.90, max=im_median*1.1, type='SQRT')

    if args.write:
        exposure.writeFits("postISRCCD.fits")
