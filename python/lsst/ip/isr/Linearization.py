"""
@brief Implementation of the stage, Linearization, for the nightly
 Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
         University of Washington
         nms@astro.washington.edu
         
 file created Mon Nov 24, 2008         

@file
"""

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.detection as det
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

def linearization(chunkExposure, isrPolicy, lookupTable):

    """Linearization
    
    @brief Correct for non-linear effects of mapping from electrons to
           ADU (amplifier roll-off). Must be done AFTER the Bias
           Subtraction stage.  Apply correction as a function of pixel
           value from either a lookup table for fitted function
           (polynomial or spline).
   
    @return chunkExposure with pixel values linearized
  
    @throw NotFound if any metadata or policy parameter can not be obtained
    @throw InvalidParameter if functional form for the lineaization fit is invalid 
    
    TO DO (as of Tue 11/24/08):
   
    - Calculate additional SDQA metrics as requested by SDQA team
    - inplement option for using satLimit from LookupTable??
    - sigma clipping?
    """
    stage = "lsst.ip.isr.linearization"   
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Linearization Stage")

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        linPolicy = isrPolicy.getPolicy("linearizePolicy")
        linearizeType = linPolicy.getString("linearizeType")
        lookupTableName = linPolicy.getString("lookupTableName")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the ISR Policy File: %s" % e   
        raise pexEx.NotFound, "Can not parse the Linearization policy file"
    
    # Get the MaskedImage and Metadata
    
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()

    try:
        pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from chunkMetadata.")

       
    except pexEx.LsstExceptionStack, e:
        print "Cannot get the ... from the chunkExposure: %s" % e      


    # currently accepts as "FUNCTION" either a polynomial or spline
    if linearizeType == "FUNCTION":

        # Get the functional form and coeffs for the polynomial from the policy
        try:
            funcForm = linPolicy.getString("funcForm")
            funcOrder = linPolicy.getInt("funcOrder")
            stepSize = linPolicy.getDouble("stepSize")
        except pexEx.LsstExceptionStack, e:
            pexLog.Trace("%s" % (stage,), 4, "Cannot get FUNCTION parameters from the Policy.")
            raise pexEx.NotFound, "Can not parse the Linearization Policy File"
            
        if funcForm == "POLYNOMIAL":
            pexLog.Trace("%s" % (stage,), 4, "Entering polynomial fitting block.")
            polyFunction = afwMath.PolynomialFunction1D(funcOrder)

            # Find the best fit function to the Exposure and apply
            # this function.  Best fit determined via Chi^2
            # minimization.

            pexLog.Trace("%s" % (stage,), 4, "Performing Chi^2 minimization for functional fit for %s." % (funcForm,))
            functionFit = ipIsr.findBestFit(chunkMaskedImage, funcForm, funcOrder, stepSize)

            for i in range (0, functionFit.parameterList.size()):
                parameters[i] = functionFit.parameterList[i]

            print "Minimization Fit Parameters :", functionFit.parameterList
            polyFunction.setParameters(parameters)

            pexLog.Trace("%s" % (stage,), 4, "Fitting minimized polynomial function.")
            ipIsr.fitFunctionToImage(chunkMaskedImage, polyFunction)
            
            pexLog.Trace("%s" % (stage,), 4, "Recording ISR functional fit provenance information.")
            chunkMetadata.addProperty(dafBase.DataProperty("LIN_FN", funcForm))    
            chunkMetadata.addProperty(dafBase.DataProperty("LIN_OR", funcOrder))
            chunkMaskedImage.setMetadata(chunkMetadata)
        
        elif funcForm == "SPLINE":
            # need to add a spline function to afw/math/FunctionLibrary
            raise pexEx.NotFound, "Function SPLINE not implemented." 

        elif funcForm == "CHEBYCHEV":
            # need to add a spline function to afw/math/FunctionLibrary
            raise pexEx.NotFound, "Function CHEBYCHEV not implemented." 

        else:
            raise pexEx.NotFound, "Function not implemented. Use 'POLYNOMIAL', 'SPLINE', or 'CHEVYCHEV'." 
            
    if linearizeType == "LOOKUP":

        pexLog.Trace("%s" % (stage,), 4, "Adding lookup table 'deltas' to MaskedImage")

        #QQQ: wil the lookup tabel deltas have uncertainties?  If so,
        #the variance will change.  It is currently not adjusted.
        
        ipIsr.iterateTable(chunkMaskedImage, lookupTable)
        
        pexLog.Trace("%s" % (stage,), 4, "Recording ISR lokup table provenance information.")
        chunkMetadata.addProperty(dafBase.DataProperty("LIN_LU", lookupTableName))
        chunkMaskedImage.setMetadata(chunkMetadata)
        
    # Record the stage provenance to the Image Metadata (for now) this
    # will eventually go into an ISR_info object.

    pexLog.Trace("%s" % (stage,), 4, "Recording final ISR provenance information.")
   
    chunkMetadata.addProperty(dafBase.DataProperty("LIN_END", "Completed Successfully")) 
    chunkMaskedImage.setMetadata(chunkMetadata)

    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    pexLog.Trace("%s" % (stage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - ?                     
    """
    
    pexLog.Trace("%s" % (stage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (stage,))
                              
