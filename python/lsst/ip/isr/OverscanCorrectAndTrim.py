"""
@brief: Implementation of the stage, Overscan Correct and Trim, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
 contact: nms@astro.washington.edu
 file created: Thu Nov 20, 2008  

@file
"""
import eups
import os

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run this program.")
overscanExposureOutPath = os.path.join(dataDir, "overscanStripTestExposure_1")

def stringParser(dataString):

    """
    Simple code to parse data sections of CCD headers.
    dataString is any data section in the form '[colBegin:colEnd,rowBegin,rowEnd]'

    Return a vw::BBox2i, colBegin, colEnd, rowBegin, rowEnd
    """

    dsString = dataString[1:-1]
    xdsPart, ydsPart = dsString.split(",")
    colb, cole = xdsPart.split(":")
    print "ColBegin, ColEnd :", colb, cole
    rowb, rowe = ydsPart.split(":")
    print "RowBegin, RowEnd: ", rowb, rowe
    colEnd = int(cole)
    colBegin = int(colb)
    rowEnd = int(rowe)
    rowBegin = int(rowb)
    colSpan = colEnd - colBegin
    print "ColSpan: ", colSpan
    rowSpan = rowEnd - rowBegin
    print "RowSpan: ", rowSpan 
    dataBbox = afwImage.BBox2i(colBegin, rowBegin, colSpan, rowSpan)
    return dataBbox, colBegin, colEnd, rowBegin, rowEnd


def overscanCorrectAndTrim(chunkExposure, isrPolicy):

    """Overscan Correct and Trim
    
    Subtract a single constant (using all input pixels for the
    statistic) or one-dimensional function from the Chunk Exposure
    which varies along the overscan.  Use spline or polynomial fit
    along the overscan region (want a slowly varying function).  For
    the 1-D representation, Pan-STARRS determines the input values to
    the fit as representations on the coordinate along the overscan,
    with the statistic derived from the pixels in the perpendicular
    direction at each location.  Sigma-clipping on input data is a
    necessary option.
  
    The Image is also trimmed to remove the overscan region and any
    unexposed pixels.
  
    @return chunkExposure corrected for overscan and trimmed of all
    non-illuminated pixels
  
    @throw Runtime if sage has already been run on the image
    @throw Runtime if any Policy or metadata value can not be obtained
    @throw InvalidParameter if functional form for the overscan fit is invalid 
   
    TO DO (as of 11/20/08):
    - Implement sigma-clipping
    - Trim ramp-up on overscan (do not include in fits)
    - Implement additional constant modes and spline fit
    - add a smoothing option to the functional fit? 
    - add any additional SDQA statistics requested by SDQA team
    """

    stage = "lsst.ip.isr.overscanCorrectAndTrim"   
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Overscan Correct and Trim")

    # Check for an overscan region.  If there isn't one, trim the
    # image, issue a logging message, and end the function

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File.")
    try:
        overPolicy = isrPolicy.getPolicy("overscanPolicy")
        method = overPolicy.getString("method")
        smooth = overPolicy.getBool("smooth")
        sigClip = overPolicy.getBool("sigClip")
        datasecKey = overPolicy.getString("datasec")
        biassecKey = overPolicy.getString("biassec")
        trimsecKey = overPolicy.getString("trimsec")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the ISR Policy File: %s" % e   

    try:
        if method == "function":
            funcForm = overPolicy.getString("funcForm")
            funcOrder = overPolicy.getInt("funcOrder")
            stepSize = overPolicy.getDouble("stepSize")
        else: 
            constantMeth = overPolicy.getString("constantMeth")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the Overscan Policy File for 'method': %s" % e

    try:
        if smooth == "true":
            smoothOpt = overPolicy.getString("smoothOpt")
            if smoothOpt == "gaussian":
                smoothSigma = overPolicy.getInt("smoothSigma")
            else:
                smoothOrder = overPolicy.getInt("smoothOrder")
    except pexEx.LsstExceptionStack, e:
         print "Cannot parse the Overscan Policy File for 'smooth': %s" % e

    try:
        if sigClip == "true":
            sigClipVal = overPolicy.getInt("sigClipVal") 
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the Overscan Policy File for 'sigClip': %s" % e
    
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()

    try:
        pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from chunkMetadata.")
        datasecProp = chunkMetadata.findUnique(datasecKey)
        datasec = datasecProp.getValueString()
        print "DATASEC: ", datasec
        trimsecProp = chunkMetadata.findUnique(trimsecKey)
        if trimsecProp:
            trimsec = trimsecProp.getValueString()
            print "TRIMSEC: ", trimsec
        biassecProp = chunkMetadata.findUnique(biassecKey)
        biassec = biassecProp.getValueString()
        print "BIASSEC: ", biassec
    except pexEx.LsstExceptionStack, e:
        print "Cannot get the requested data sections from the chunkExposure: %s" % e      

    datasecBbox, dcb, dce, drb, dre = stringParser(datasec)
    biassecBbox, bcb, bce, brb, bre = stringParser(biassec)
    if trimsec:
        trimsecBbox, tcb, tce, trb, tre = stringParser(trimsec)

    else: 	
        tcb = dcb - bce
        tce = dce - bce
        trb = drb - brb
        tre - dre - bre
        trimsecBbox = afwImage.BBox2i(tcb, trb, (tce - tcb), (tre - trb))
    
    overscanExposure = afwImage.ExposureD()
    overscanMaskedImage = afwImage.MaskedImageD()
    overscanExposure = chunkExposure.getSubExposure(biassecBbox)
    overscanMaskedImage = chunkExposure.getMaskedImage()
    pexLog.Trace("%s" % (stage,), 4, "Writing Overscan Exposure to FitsStorage.")
    overscanExposure.writeFits(overscanExposureOutPath)

    trimmedExposure = afwImage.ExposureD()
    if trimsec:
        trimmedExposure = chunkExposure.getSubExposure(trimsecBbox)
    else:
        trimmedExposure = chunkExposure.getSubExposure(trimExp)

    # create a MaskedImage holding the overscan region to be
    # subtracted QQQ: NEED TO TRIM RAMP UP HERE in addition to
    # overscan region because we do not want to use the ramp in the
    # overscan subtraction. How do I determine appropriate ramp to
    # trim (need to create a Bbox)?
    
    trimmedMaskedImage = afwImage.MaskedImageD()
    trimmedMaskedImage = trimmedExposure.getMaskedImage()
    trimmedChunkMetadata = trimmedMaskedImage.getImage().getMetaData()
    trimmedVarianceMetadata = trimmedMaskedImage.getVariance().getMetaData()

    # Remove the datasec and biassec from the trimmed Chunk Exposure's
    # metadata...which currently lives in both the image and the
    # variance image.

    pexLog.Trace("%s" % (stage,), 4, "Removing DATASEC and BIASSEC keywords from metadata.")
    trimmedChunkMetadata.deleteAll(datasecKey)
    trimmedChunkMetadata.deleteAll(biassecKey)
    trimmedVarianceMetadata.deleteAll(datasecKey)
    trimmedVarianceMetadata.deleteAll(biassecKey)

    if trimsec:
        trimmedChunkMetadata.deleteAll(trimsecKey)
        trimmedVarianceMetadata.deleteAll(trimsecKey)
        
    trimmedMaskedImage.setMetadata(trimmedChunkMetadata)
    trimmedMaskedImage.setMetadata(trimmedVarianceMetadata)

    if method == "FUNCTION":
        if funcForm == "POLYNOMIAL":
            # Find the best fit function to the overscan region and apply
            # this function to the Chunk Exposure.  Best fit determined
            # via Chi^2 minimization.

            print "Performing Chi^2 minimization for functional fit for ", funcForm
            overscanFit = ipIsr.findBestFit(overscanMaskedImage, funcForm, funcOrder, stepSize)
        
            for i in range (0, overscanFit.parameterList.size()):
                parameters[i] = overscanFit.parameterList[i]      

                print "Minimization Fit Parameters :", overscanFit.parameterList
                order = overscanFit.parameterList.size() - 1
                polyFunction = afwMath.PolynomialFunction1D()
                polyFunction(order)
                polyFunction.setParameters(parameters)
        
                ipIsr.fitFunctionToImage(trimmedMaskedImage, polyFunction)

                # Record the polynomial provenance information

                pexLog.Trace("%s" % (stage,), 4, "Recording ISR provenance information for %s." % (funcForm))
                trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_FN", funcForm))
                trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_OR", order))
                trimmedMaskedImage.setMetadata(trimmedChunkMetadata)
        
        elif funcForm == "SPLINE":
            # need to add a spline function to afw/math/FunctionLibrary
            raise pexEx.NotFound, "Function SPLINE not implemented." 

        elif funcForm == "CHEBYCHEV":
            # need to add a spline function to afw/math/FunctionLibrary
            raise pexEx.NotFound, "Function CHEBYCHEV not implemented." 

        else:
            raise pexEx.NotFound, "Function not implemented. Use 'POLYNOMIAL', 'SPLINE', or 'CHEVYCHEV'." 
       
    else: 
        
        # Subtract a constant value.  For now, compute the mean for
        # the constant value to be subtracted. ADD THE OTHERS LATER...
        
        if constantMeth == "MEAN":

            pexLog.Trace("%s" % (stage,), 4, "Computing mean for overscan region.")
            mu = ipIsr.easyMean(overscanMaskedImage)
	    print "Mean Value: ", mu
            
            pexLog.Trace("%s" % (stage,), 4, "Subtracting mean from MaskedImmge.")
            trimmedMaskedImage -= mu

            pexLog.Trace("%s" % (stage,), 4, "Recording ISR provenance information for mean.")
            trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_FN", constantMeth)) 
            trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_MU", mu))    
            # trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_SD", sigma))
        trimmedMaskedImage.setMetadata(trimmedChunkMetadata)
            
        elif constantMeth == "MEDIAN":
            # not yet implemented
            raise pexExcept.InvalidParameter, "Median is not yet implemented."
        elif constantMeth == "MODE":
            # not yet implemented
            raise pexExcept.InvalidParameter, "Mode is not yet implemented."
        else: 
            raise pexExcept.InvalidParameter, "Invalid method requested for computing the constant overscan value."
        
    # Record final sub-stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final ISR stage provenance information.")
    trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_END", "Completed Successfully")) 
    trimmedMaskedImage.setMetadata(trimmedChunkMetadata);

    # Set the original chunkExposure to the new trimmed chunkExposure
    chunkMaskedImage = trimmedMaskedImage
    chunkExposure = trimmedExposure
    
    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    pexLog.Trace("%s" % (stage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - n
    - mu
    - sigma
    """                       
    pexLog.Trace("%s" % (stage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (stage,))
