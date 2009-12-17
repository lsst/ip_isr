#!/usr/bin/env python

import re

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
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


