from lsst.ip.isr import AssembleCcdTask
import lsst.afw.display.ds9 as ds9
import exampleUtils


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
        assemblyInput = exampleUtils.makeAssemblyInput(isPerAmp)
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
