#!/usr/bin/env python

import re

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase

from lsst.pex.harness.Stage import Stage

propertySetTypeInfos = {}

def setTypeInfos():
    global propertySetTypeInfos
    p = dafBase.PropertySet()
    p.set("str",  "bar"); propertySetTypeInfos["string"] = p.typeOf("str")
    p.set("int",      3); propertySetTypeInfos["int"]    = p.typeOf("int")
    p.set("float",  3.1); propertySetTypeInfos["float"]  = p.typeOf("float")
    p.set("bool",  True); propertySetTypeInfos["bool"]   = p.typeOf("bool")
setTypeInfos()

def validateMetadata(metadata, metadataPolicy):
    paramNames = metadataPolicy.paramNames(1)
    for paramName in paramNames:
        if not metadata.exists(paramName):
            raise RuntimeError, 'Unable to find \'%s\' in metadata' % (paramName,)
        # TBD; VALIDATE AGAINST DICTIONARY FOR TYPE ETC
    return True

def transformMetadata(metadata, datatypePolicy, metadataPolicy, suffix):
    paramNames = metadataPolicy.paramNames(1)
    for paramName in paramNames:
        # If it already exists don't try and update it
        if metadata.exists(paramName):
            continue
        
        mappingKey = paramName + suffix
        if datatypePolicy.exists(mappingKey):
            keyword = datatypePolicy.getString(mappingKey)
            if metadata.typeOf(keyword) == propertySetTypeInfos["string"]:
                val = metadata.getString(keyword).strip()

                # some FITS files have extra words after the field name
                if paramName == "datasetId" and val.find(' ') > 0:
                    val = val[:val.index(' ')]

                metadata.set(paramName, val)
            else:
                metadata.copy(paramName, metadata, keyword)
#            metadata.copy(paramName, metadata, keyword)
            metadata.copy(keyword + "_original", metadata, keyword)
            metadata.remove(keyword)
    
    # Any additional operations on the input data?
    if datatypePolicy.exists('convertDateobsToTai'):
        if datatypePolicy.getBool('convertDateobsToTai'):
            dateObs  = metadata.getDouble('dateObs')
            dateTime = dafBase.DateTime(dateObs, dafBase.DateTime.UTC)
            dateObs  = dateTime.mjd(dafBase.DateTime.TAI)
            metadata.setDouble('dateObs', dateObs)

    if datatypePolicy.exists('convertDateobsToMidExposure'):
        if datatypePolicy.getBool('convertDateobsToMidExposure'):
            dateObs  = metadata.getDouble('dateObs')
            dateObs += metadata.getDouble('expTime') * 0.5 / 3600. / 24.
            metadata.setDouble('dateObs', dateObs)

    dateTime = dafBase.DateTime(metadata.getDouble('dateObs'))
    metadata.setDateTime('taiObs', dateTime)

    if datatypePolicy.exists('trimFilterName'):
        if datatypePolicy.getBool('trimFilterName'):
            filter = metadata.getString('filter')
            filter = re.sub(r'\..*', '', filter)
            metadata.setString('filter', filter)

    if datatypePolicy.exists('convertVisitIdToInt'):
        if datatypePolicy.getBool('convertVisitIdToInt'):
            visitId  = metadata.getString('visitId')
            metadata.setInt('visitId', int(visitId))


class ValidateMetadataStage(Stage):

    """Validates that every field in metadataPolicy exists in the
    input metadata, before sending the event down the pipeline.  This
    will evolve as the pipeline evolves and the metadata requirements
    of each stage evolves.

    For the input of external (non-LSST) data these data's metadata
    should be generally be run through the TransformMetadata stage,
    with a survey-specific policy file specifying this mapping.
    """

    def process(self):
        clipboard = self.inputQueue.getNextDataset()
        metadataPolicy = self._policy.getPolicy("metadata")
        imageMetadataKey = self._policy.get("imageMetadataKey")
        metadata = clipboard.get(imageMetadataKey)
        validateMetadata(metadata, metadataPolicy)
        self.outputQueue.addDataset(clipboard)
    
class TransformMetadataStage(Stage):

    """This stage takes an input set of metadata and transforms this
    to the LSST standard.  It will be input-dataset specific, and the
    mapping is described in the datatypePolicy.  The standard is to
    have a string in the datatypePolicy named metadataKeyword that
    represents the location of LSST metadata in the particular data
    set."""

    def process(self):
        clipboard = self.inputQueue.getNextDataset()
        metadataPolicy = self._policy.getPolicy("metadata")
        datatypePolicy = self._policy.getPolicy("datatype")
        imageKey = self._policy.get("imageKey")
        metadataKey = self._policy.get("metadataKey")
        exposureMetadataKey = self._policy.get("exposureMetadataKey")
        wcsKey = self._policy.get("wcsKey")
        decoratedImage = clipboard.get(imageKey)
        exposureMetadata = clipboard.get(exposureMetadataKey)
        metadata = decoratedImage.getMetadata()

        if self._policy.exists("suffix"):
            suffix = self._policy.get("suffix")
        else:
            suffix = "Keyword"

        if self._policy.exists("computeWcsGuess"):
            if self._policy.getBool("computeWcsGuess"):
                wcs = afwImage.Wcs(metadata)
                ampBBoxKey = self._policy.getString("ampBBoxKey")
                ampBBox = clipboard.get(ampBBoxKey)
                wcs.shiftReferencePixel(ampBBox.getX0(), ampBBox.getY0())
                clipboard.put(wcsKey, wcs)

        transformMetadata(metadata, datatypePolicy, metadataPolicy, suffix)

        metadata.setLongLong('ampExposureId',
                exposureMetadata.get('ampExposureId'))
        metadata.setLongLong('ccdExposureId',
                exposureMetadata.get('ccdExposureId'))
        metadata.setLongLong('fpaExposureId',
                exposureMetadata.get('fpaExposureId'))
        metadata.set('url', metadata.get('filename'))
        metadata.set('ampId', clipboard.get('ampId'))

        clipboard.put(metadataKey, metadata)
        clipboard.put(imageKey, decoratedImage.getImage())

        self.outputQueue.addDataset(clipboard)

class TransformExposureMetadataStage(Stage):

    """This stage takes a list of input Exposures and transforms the metadata
    of each one to the LSST standard.  It will be input-dataset specific, and
    the mapping is described in the datatypePolicy.  The standard is to have a
    string in the datatypePolicy named metadataKeyword that represents the
    location of LSST metadata in the particular data set."""

    def process(self):
        clipboard = self.inputQueue.getNextDataset()
        metadataPolicy = self._policy.getPolicy("metadata")
        datatypePolicy = self._policy.getPolicy("datatype")
        exposureKeys = self._policy.getStringArray("exposureKey")
        ampBBoxKey = self._policy.getString("ampBBoxKey")
        ampBBox = clipboard.get(ampBBoxKey)

        if self._policy.exists("suffix"):
            suffix = self._policy.get("suffix")
        else:
            suffix = "Keyword"

        for exposureKey in exposureKeys:
            exposure = clipboard.get(exposureKey)
            metadata = exposure.getMetadata()
            transformMetadata(metadata, datatypePolicy, metadataPolicy, suffix)
            exposure.getMaskedImage().setXY0(ampBBox.getLLC())

        self.outputQueue.addDataset(clipboard)
