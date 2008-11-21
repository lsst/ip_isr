
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

    overStage = "lsst.ip.isr.overscanCorrectAndTrim"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (overStage,))

    # Parse the Policy File
    try:
     overscanPolicy = isrPolicy.getPolicy("overscanPolicy")
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
            numParams = overPolicy.getDouble("numParams")
        else: 
            constantMeth = overPolicy.getString("constantMeth")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the Overscan Policy File for 'method': %s" % e

     try:
         if smooth = true:
             smoothOpt = overPolicy.getString("smoothOpt")
             if smoothOpt == "gaussian":
                 smoothSigma = overPolicy.getInt("smoothSigma")
             else:
                 smoothOrder = overPolicy.getInt("smoothOrder")
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the Overscan Policy File for 'smooth': %s" % e

     try:
         if sigClip = true:
             sigClipVal = overPolicy.getInt("sigClipVal") 
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the Overscan Policy File for 'sigClip': %s" % e
    
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()

    try:
        pexLog.Trace("%s" % (satStage,), 4, "Obtaining additional parameters from chunkMetadata.")
        datasecProp = chunkMetadata.findUnique(datasecKey)
        datasec = datasecProp.getValue()
        print "DATASEC: ", datasec
        trimsecProp = chunkMetadata.findUnique(trimsecKey)
        if trimsecProp:
            trimsec = trimsecProp.getValue()
            print "TRIMSEC: ", trimsec
        biassecProp = chunkMetadata.findUnique(biassecKey)
        biassec = biassecProp.getValue()
        print "BIASSEC: ", biassec
    except pexEx.LsstExceptionStack, e:
        print "Cannot get the requested data sections from the chunkExposure: %s" % e      

    datasecBbox = ipIsr.stringParse(datasec)
    biassecBbox = ipIsr.stringParse(biassec)
    if trimsec:
        trimsecBbox = ipIsr.stringParse(trimsec)

    trimExp = datasecBbox - biassecBbox

    overscanExposure = afwImage.ExposureD()
    overscanMaskedImage = afwImage.MaskedImageD()
    overscanExposure = chunkExposure.getSubExposure(biassecBbox)
    overscanMaskedImage = chunkExposure.getMaskedImage()
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
    trimmedVarianceMetadata = trimmedMaskedImage.getVariance()->getMetaData();

    # Remove the datasec and biassec from the trimmed Chunk Exposure's
    # metadata...which currently lives in both the image and the
    # variance image.

    trimmedChunkMetadata.deleteAll(datasecKey)
    trimmedChunkMetadata.deleteAll(biassecKey)
    trimmedVarianceMetadata.deleteAll(datasecKey)
    trimmedVarianceMetadata.deleteAll(biassecKey)

    if trimsec:
        trimmedChunkMetadata.deleteAll(trimsecKey)
        trimmedVarianceMetadata.deleteAll(trimsecKey)
        
    trimmedMaskedImage.setMetadata(trimmedChunkMetadata)
    trimmedMaskedImage.setMetadata(trimmedVarianceMetadata)

    if method == "function":

        # Find the best fit function to the overscan region and apply
        # this function to the Chunk Exposure.  Best fit determined
        # via Chi^2 minimization.

        std.fill(parameterList.begin(), parameterList.end(), numParams)
        std.fill(stepSizeList.begin(), stepSizeList.end(), stepSize)
            
        # collapse the overscan region down to a vector.  Get vectors
        # of measurements, variances, and x,y positions
        
        ipIsr.makeVectorFromRegion(overscanMeasurementList, overscanVarianceList, colPositionList, rowPositionList, overscanMaskedImage)

        # ideally we need a function factory to which we feed
        # different functional forms and get out a Function1 or
        # Function2...but for now, pick a few typical choices and code
        # them here (poly for now 10/22/08...ADD THE OTHERS LATER)

        if funcForm == "polynomial":

            polyFunc1 = afwMath.PolynomialFunction1D()
            polyFunc1(funcOrder)
            sigma = 1.0; // initial guess

            # find the best fit function
            overscanFit  = afwMath.minimize( 
                polyFunc1, 
                parameterList, 
                stepSizeList, 
                overscanMeasurementList, # overscan values 
                overscanVarianceList,    # variance for each value 
                colPositionList,         # x position  
                //rowPositionList,       # y position (if Function2 function)
                sigma 
                )  

            for i = 0, i < overscanFit.parameterList.size(), ++i:
                parameters[i] = overscanFit.parameterList[i]      
            
            order = overscanFit.parameterList.size() - 1
            polyFunction = afwMath.PolynomialFunction1D()
            polyFunction(order)
            polyFunction.setParameters(parameters)
        
            ipIsr.fitFunctionToImage(trimmedMaskedImage, polyFunction)
       
            # Record the stage provenance to the Image Metadata (for now) this
            # will eventually go into an ISR_info object.

            pexLog.Trace("%s" % (satStage,), 4, "Recording ISR provenance information for %s." % (funcForm))
            trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_FN", funcForm))
            trimmedMaskedImage.setMetadata(trimmedChunkMetadata)
            
        else if funcForm == "spline":
            # not yet implemented
            raise pexExcept.InvalidParameter, "Spline is not yet implemented."
        else:
            raise pexExcept.InvalidParameter, "Invalid functional form for overscan fit requested."
        
    else: 

        # Subtract a constant value.  For now, compute the mean for
        # the constant value to be subtracted. ADD THE OTHERS LATER...
        
        if constantMeth == "mean":

            ipIsr.easyMean(n, mu, sigma, overscanMaskedImage)
            trimmedMaskedImage -= mu
        else if constantMeth == "median":
            # not yet implemented
            raise pexExcept.InvalidParameter, "Median is not yet implemented."
        else if constantMeth == "mode":
            # not yet implemented
            raise pexExcept.InvalidParameter, "Mode is not yet implemented."
        else: 
            raise pexExcept.InvalidParameter, "Invalid method for computing the overscan value requested."
        trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_FN", constantMeth)) 
        trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_MU", mu))    
        trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_SD", sigma))
        trimmedMaskedImage.setMetadata(trimmedChunkMetadata)
        
    # Record final sub-stage provenance to the Image Metadata
    trimmedChunkMetadata.addProperty(dafBase.DataProperty("OVER_END", "Completed Successfully")) 
    trimmedMaskedImage.setMetadata(trimmedChunkMetadata);

    # Set the original chunkExposure to the new trimmed chunkExposure
    chunkMaskedImage = trimmedChunkMaskedImage
    chunkExposure = trimmedExposure
    
    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    # pexLog.Trace("%s" % (overStage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - n
    - mu
    - sigma
    """                       
    pexLog.Trace("%s" % (overStage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (overStage,))
                              
