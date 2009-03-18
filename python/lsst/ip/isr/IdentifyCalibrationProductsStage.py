import os

from lsst.pex.harness.Stage import Stage
from lsst.daf.base import DateTime, PropertySet
from lsst.daf.persistence import LogicalLocation

import lsst.ip.isr.calibDatabase as calibDatabase

class IdentifyCalibrationProductsStage(Stage):
    def __init__(self, stageId=-1, stagePolicy=None):
        Stage.__init__(self, stageId, stagePolicy)
        self.cdb = calibDatabase.CalibDB(self._policy.get("calibDbPath"))

    def process(self):
        self.activeClipboard = self.inputQueue.getNextDataset()

        eventName = self._policy.get("eventName")
        event = self.activeClipboard.get(eventName)
        when = DateTime(event.get("dateObs"))

        ccdId = self.activeClipboard.get("ccdId")
        ampId = self.activeClipboard.get("ampId")
        
        expTime = event.get("expTime")
        darkPolicy = self._policy.get("darkPolicy")
        darkCalibList = self.cdb.lookup(when, "dark", ccdId, ampId, all=True)
        darkTimeList = []
        for d in darkCalibList:
            darkTimeList.append(d.expTime)
        darkTimeList.sort()
        if darkPolicy == "min":
            darkExpTime = darkTimeList[0]
        elif darkPolicy == "max":
            darkExpTime = darkTimeList[-1]
        elif darkPolicy == "closest":
            minDist = abs(expTime - darkTimeList[0])
            minExpTime = darkTimeList[0]
            for i in xrange(1, len(darkTimeList)):
                dist = abs(expTime - darkTimeList[i])
                if dist < minDist:
                    minDist = dist
                    minExpTime = darkTimeList[i]
            darkExpTime = minExpTime
        else:
            raise RuntimeError, "Unrecognized darkPolicy: " + str(darkPolicy)


        biasPath = self.cdb.lookup(when, "bias", ccdId, ampId)
        darkPath = self.cdb.lookup(when, "dark", ccdId, ampId,
                expTime=darkExpTime)
        defectPath = self.cdb.lookup(when, "defect", ccdId, ampId)
        flatPath = self.cdb.lookup(when, "flat", ccdId, ampId,
                filter=event.get("filter"))
#         fringePath = self.cdb.lookup(when, "fringe", ccdId, ampId,
#                 filter=event.get("filter"))
        linearizePath = self.cdb.lookup(when, "linearize", ccdId, ampId)

        pathPrefix = LogicalLocation(self._policy.get("pathPrefix")).locString()
        calibData = PropertySet()
        calibData.set("biasPath", os.path.join(pathPrefix, biasPath))
        calibData.set("darkPath", os.path.join(pathPrefix, darkPath))
        calibData.set("defectPath", os.path.join(pathPrefix, defectPath))
        calibData.set("flatPath", os.path.join(pathPrefix, flatPath))
#         calibData.set("fringePath", os.path.join(pathPrefix, fringePath))
        calibData.set("linearizePath", os.path.join(pathPrefix, linearizePath))

        outputKey = self._policy.get("outputKey")
        self.activeClipboard.put(outputKey, calibData)

        self.outputQueue.addDataset(self.activeClipboard)
