from lsst.ip.isr import IsrTask
import lsst.afw.display.ds9 as ds9
import exampleUtils

def runIsr():
    '''Run the task to do ISR on a ccd'''

    #Create the isr task with modified config
    isrConfig = IsrTask.ConfigClass()
    isrConfig.doBias = False #We didn't make a zero frame
    isrConfig.doDark = True
    isrConfig.doFlat = True
    isrConfig.doFringe = False #There is no fringe frame for this example
    isrConfig.doWrite = True

    isrConfig.assembleCcd.doRenorm = False #We'll take care of gain in the flats
    isrConfig.assembleCcd.setGain = False 
    isrTask = IsrTask(config=isrConfig)

    sensorRef = exampleUtils.FakeDataRef()
    output = isrTask.run(sensorRef)
    return output.exposure

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of IsrTask")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--displayResult', '-R', action="store_true", help="Display the result?", default=False)
    args = parser.parse_args()

    if args.debug:
        try:
            import debug
        except ImportError as e:
            print >> sys.stderr, e
    exposure = runIsr()
    if args.displayResult:
        ds9.mtv(exposure.getMaskedImage())
