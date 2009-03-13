from lsst.pex.harness.Stage import Stage
from lsst.daf.base import DateTime

import lsst.ip.isr.calibDatabase as cdb

class IdentifyCalibrationProductsStage(Stage):
    def process(self):
        self.activeClipboard = self.inputQueue.getNextDataset()

        event = self.activeClipboard.get("triggerImageprocEvent0")
        when = DateTime(event.get("dateObs"))
        
        if self.cdb is None:
            self.cdb = calibDatabase.CalibDB(self._policy.get("dbPath"))

        expTime = event.get("expTime")
        darkPolicy = self._policy.get("darkPolicy")
        darkCalibList = cdb.lookup(when, "dark", event.get("ccdId"),
                event.get("ampId"), all=True)
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


        biasPath = cdb.lookup(when, "bias", event.get("ccdId"),
                event.get("ampId"))
        darkPath = cdb.lookup(when, "dark", event.get("ccdId"),
                event.get("ampId"), expTime=darkExpTime)
        defectPath = cdb.lookup(when, "defect", event.get("ccdId"),
                event.get("ampId"))
        flatPath = cdb.lookup(when, "flat", event.get("ccdId"),
                event.get("ampId"), filter=event.get("filter"))
        fringePath = cdb.lookup(when, "fringe", event.get("ccdId"),
                event.get("ampId"), filter=event.get("filter"))
        linearizePath = cdb.lookup(when, "linearize", event.get("ccdId"),
                event.get("ampId"))

        calibData = PropertySet()
        calibData.set("bias", biasPath)
        calibData.set("dark", darkPath)
        calibData.set("defect", defectPath)
        calibData.set("flat", flatPath)
        calibData.set("fringe", fringePath)
        calibData.set("linearize", linearizePath)
        self.activeClipboard.put("calibData", calibData)

        self.outputQueue.addDataset(self.activeClipboard)
