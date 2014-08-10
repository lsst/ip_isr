from lsst.ip.isr import AssembleCcdTask
import lsst.afw.display.ds9 as ds9
import exampleUtils

def makeAssemblyInput(isPerAmp):
    '''!Make the input to pass to the assembly task
    \param[in] isPerAmp -- If True, return a dictionary of amp exposures keyed by amp name.
                           If False, return a single exposure with amps mosaiced preserving non-science pixels
                           (e.g. overscan)
    \return Either a dictionary of amp exposures or an exposure contining the mosaiced amps.
    '''

    #number of amps in x and y
    nAmpX = 3
    nAmpY = 2

    #number of science pixels in each amp in x and y
    nPixX = 512
    nPixY = 1024

    #number of prescan rows
    pre = 4

    #number of horizontal overscan columns
    hOscan = 10

    #number of vertical overscan rows
    vOscan = 15

    #number of pixels in the extended register
    ext = 1

    #First get the per amp input data and assemble if necessary
    detector = exampleUtils.createDetector(nAmpX, nAmpY, nPixX, nPixY, pre, hOscan, vOscan, ext, True)
    inputData = exampleUtils.makeAmpInput(detector)
    if not isPerAmp:
        noTrimAssembleConfig = AssembleCcdTask.ConfigClass()
        noTrimAssembleConfig.doTrim = False #Preserve non-science pixels
        noTrimAssembleTask = AssembleCcdTask(config=noTrimAssembleConfig)
        ccdAssemblyInput = noTrimAssembleTask.assembleCcd(inputData)
        #create a detector describing a mosaiced amp grid and set it on output data
        detector = exampleUtils.createDetector(nAmpX, nAmpY, nPixX, nPixY, pre, hOscan, vOscan, ext, isPerAmp)
        ccdAssemblyInput.setDetector(detector)
        return ccdAssemblyInput
    else:
        return inputData

def runAssembler():
    '''Run the task to assemble amps into a ccd'''

    #Create the assemble task with default config
    assembleConfig = AssembleCcdTask.ConfigClass()
    assembleTask = AssembleCcdTask(config=assembleConfig)
    frame = 0

    #The assemble task can operate in two ways:
    #1. On a dictionary of amp size images
    #2. On a single amp mosaic image
    #Both methods should produce the same output.
    for isPerAmp in (True, False):
        assemblyInput = makeAssemblyInput(isPerAmp)
        assembledExposure = assembleTask.assembleCcd(assemblyInput)
        ds9.mtv(assembledExposure.getMaskedImage(), frame=frame, title="Per amp input is %s"%(isPerAmp))
        frame += 1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of AssembleCcdTask")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    args = parser.parse_args()

    if args.debug:
        try:
            import debug
        except ImportError as e:
            print >> sys.stderr, e
    runAssembler()
