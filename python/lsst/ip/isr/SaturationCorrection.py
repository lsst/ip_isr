"""

@brief: Implementation of the stage, Saturation Correction, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
 contact: nms@astro.washington.edu
 file created: 

@file
"""

TODO: replace double gaussian psf
TODO: replace lsst.detection.defects
import lsst.detection as det
import lsst.detection.defects as defects

import os
import math
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as measAlg
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import isrLib

def saturationCorrection(chunkExposure, isrPolicy, lookupTable):

    """Saturation Correction for the Chunk Exposure

    Detect and mask pixels that are saturated in the A/D converter or
    have excessive non-linear response in the Chunk Exposure.  Uses
    lsst.afw.detection.DetectionSet to create footprints of
    saturated pixels above a threshold as given in the ISR Policy
    File, a lookup table or given in the metadata as the 'SATURATE'
    field.

    Lookup table is assumed to be in two column (comma delimited)
    format.  First column holds the chunkID (should match the metadata
    ID) and the second holds the threshold value for the chunk.

    Grow by additional pixels (as given in the ISR Policy) to mask
    charge spillover.  Set appropriate bits in the Mask and
    interpolate over masked pixels using lsst.meas.algorithms.interpolateOverDefects.
    A simple PSF estimation is used for now.

    @return chunkExposure with saturated pixels masked, grown, and interpolated
  
    @throw NotFound if any metadata or policy parameter can not be obtained
   
    TO DO (as of Tue 12/11/08):
    - delineate between A/D saturated pixels and other?
    - mask grown sat pixels with a different bitmask than for sat pix?
    - Calculate additional SDQA metrics as requested by SDQA team
    """

    stage = "lsst.ip.isr.saturationCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
        pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File.")
        satPolicy = isrPolicy.getPolicy("saturationPolicy")
        useLookup = satPolicy.getBool("satLookup")
        satTable = satPolicy.getBool("satTable")
        useDefSat = satPolicy.getBool("useDefSat")
        satGrow = satPolicy.getInt("grow")
        satThresh = satPolicy.getInt("threshold")
        nSatPixMin = satPolicy.getInt("nSatPixMin")
        psfFwhm = satPolicy.getInt("psfFWHM")
        chunkType = isrPolicy.getString("chunkType")
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Can not parse the ISR Policy File." % (e,)) 
        raise pexExcept.NotFound, "In %s: Can not obtain policy parameters from the ISR Policy File." % (stage,)

    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()

    try:
        pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the Metadata.")
        satField = chunkMetadata.findUnique("SATURATE")
        satLimit = satField.getValueInt()
        pexLog.Trace("%s" % (stage,), 4, "SATURATE: %s" % (satLimit,)) 
        print "satLimit: %s" % (satLimit,)
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Cannot get SATURATE from the Metadata." % (e,)) 
        raise pexExcept.NotFound, "In %s: Can not obtain SATURATE from the Metadata." % (stage,)  

    chunkMask = chunkMaskedImage.getMask()
    satMaskBit = chunkMask.getPlaneBitMask("SAT")
    growMaskBit = chunkMask.getPlaneBitMask("GROW")

   
    # Determine the threshold value for the chunk
    if useLookup == "true":
        pexLog.Trace("%s" % (stage,), 4, "Using Lookup Table: %s" % (satTable,))
        if chunkType == "amp":
            ampidField = chunkMetadata.findUnique("AMPLIST")
            if ampidField:
                id = ampidField.getValueString()
            else:
                raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)
        
        else:
            ccdidField = chunkMetadata.findUnique("CCDLIST")
            if ccdidField:
                id = ccdidField.getValueString()
            else:
                raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

        # assuming table has two comma-delimited columns, the first
        # holds the chunkID, the second holds the threshold value.
 
        lookupTable = open(satTable, "rU")  
        chunkRows = lookupTable.readlines()
        numChunks = len(chunkRows)
        print 'Number of chunks: ', numChunks
        for chunks in chunkRows:
            # strip trailing whitespace, returns, etc.
            chunks = chunks.strip()
            # ignore blank lines
            if not chunks:
                continue
            # ignore comment lines
            if chunks.startswith("#"):
                continue
            chunkList = chunks.split()
            if len(chunkList) < numChunks or len(chunkList) > numChunks:
                raise pexExcept.LengthError, "In %s: Could not parse the Lookup Table." %(stage,)
            
        for ids in chunkList:
            try:
                idField, thresh = ids.split(",")
            except Exception, e:
                pexLog.Trace("%s" % (stage,), 4, "Can not parse the Lookup Table." % (e,)) 
                raise pexExcept.LengthError, "In %s: Could not parse the Lookup Table." %(stage,)

        maxInd = numChunks - 1
        print "MaxInd: ", maxInd
        for ind in range 0, maxInd:
            if id = idField[ind]:
                thresholdVal = thresh[ind]
                break
            else:
                continue
            
    elif useDefSat == "true":
        thresholdVal = satThresh
    else:
        thresholdVal = satLimit
        
    # QQQ: Can we distinguish between pixels saturated in the A/D
    # converter and those just saurated on chip?  If so, we don't want
    # to grow around the A/D saturated pixels (in unbinned data).

    pexLog.Trace("%s" % (stage,), 4, "Finding saturated footprints.")
    satFootprintSet = afwDet.DetectionSetD(chunkMaskedImage, afwDet.Threshold(thresholdVal))
    satFootprintList = satFootprintSet.getFootprints()
    numSatFootprints = len(satFootprintList)
    pexLog.Trace("%s" % (stage,), 4, "Found %s saturated footprints." % (numSatFootprints,))
  
    numSatPix = 0
    for feet in satFootprintList:
        numSatPix = numSatPix + feet.getNpix()

    pexLog.Trace("%s" % (stage,), 4, "Found %s saturated pixels." % (numSatPix,))
    
    pexLog.Trace("%s" % (stage,), 4, "Growing around all saturated footprints.")
    grownSatFootprintSet = afwDet.DetectionSetD(satFootprintSet, satGrow)
    grownSatFootprintList = grownSatFootprintSet.getFootprints()
    numGrownSatFootprints = len(grownSatFootprintList)
    pexLog.Trace("%s" % (stage,), 4, "Found %s grown saturated footprints." % (numGrownSatFootprints,))

    numGrownSatPix = 0
    for bigFeet in grownSatFootprintList:
        numGrownSatPix = numGrownSatPix + bigFeet.getNpix()
    
    pexLog.Trace("%s" % (stage,), 4, "Found %s saturated pixels." % (numGrownSatPix,))
    
    # Mask all of the grown saturated pixel footprints.  Using "SAT"
    # bitmask for all pixels in the footprint.

    # QQQ: Do we want to distinguish between grown pixels and those
    # that were actually saturated?  Use "GROW" bitmask above for
    # grown pixels only?
    
    pexLog.Trace("%s" % (stage,), 4, "Masking all saturated pixels.")
    afwDet.setMaskFromFootprintList(chunkMask, grownSatFootprintList, satMaskBit)

    pexLog.Trace("%s" % (stage,), 4, "Interpolating over all saturated footprints.")  
    psf = det.dgPSF(psfFwhm/(2*math.sqrt(2*math.log(2)))) 
    chunkMask.addMaskPlane("INTERP")

    badPixList = defects.policyToBadRegionList(os.path.join(os.environ["DETECTION_DIR"], "pipeline/BadPixels.paf"))

    measAlg.interpolateOverDefects(chunkMaskedImage, psf, badPixList)
    #measAlg.interpolateOverDefects(chunkMaskedImage, psf, grownSatFootprintList)
    
    # Record the stage provenance to the Image Metadata (for now) this
    # will eventually go into an ISR_info object.

    pexLog.Trace("%s" % (stage,), 4, "Recording ISR provenance information.")
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_TH", satThresh))
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_PIX", numSatPix))    
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_FP", numSatFootprints))
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_END", "Completed Successfully")) 
    chunkMaskedImage.setMetadata(chunkMetadata)

    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    pexLog.Trace("%s" % (stage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - numSatFootprints
    - numSatPix
    - numGrownSatFootprints
    - numGrownSatPix                          
    """
    
    pexLog.Trace("%s" % (stage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (stage,))
                              
