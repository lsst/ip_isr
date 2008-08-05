import lsst.pex.harness.Stage
	import lsst.daf.base as dafBase
	import lsst.daf.persistence as dafPers
	import lsst.afw.image
	
	class MetadataManipStage(lsst.pex.harness.Stage.Stage):
	    def preprocess(self):
	        print 'Python MetadataManipStage preprocess'
	        self.activeClipboard = self.inputQueue.getNextDataset()
	
	        inputDP = self.activeClipboard.get("triggerVisitEvent")
	
	        visitDP = dafBase.DataProperty.createPropertyNode("visit2exposure")
	        visitDP.addProperty(dafBase.DataProperty("visitId",
	            inputDP.findUnique("visitId").getValueInt()))
	        visitDP.addProperty(dafBase.DataProperty.createInt64DataProperty(
	            "exposureId", inputDP.findUnique("exposureId").getValueInt() << 1))
	        self.activeClipboard.put("visit2exposure", visitDP)
	
	        outputDP = dafBase.DataProperty.createPropertyNode("rawFPAExposure")
	
	        outputDP.addProperty(dafBase.DataProperty.createInt64DataProperty(
	            "rawFPAExposureId",
	            inputDP.findUnique("exposureId").getValueInt() << 1))
	        outputDP.addProperty(dafBase.DataProperty("ra",
	            inputDP.findUnique("FOVRA").getValueDouble()))
	        outputDP.addProperty(dafBase.DataProperty("decl",
	            inputDP.findUnique("FOVDec").getValueDouble()))
	        outputDP.addProperty(dafBase.DataProperty.createFloatDataProperty("equinox",
	            inputDP.findUnique("equinox").getValueDouble()))
	        outputDP.addProperty(dafBase.DataProperty("filterId",
	            self._lookup(inputDP.findUnique("filterName").getValueString())))
	        outputDP.addProperty(dafBase.DataProperty.createDateTimeDataProperty(
	            "dateObs",
	            dafPers.DateTime(inputDP.findUnique("visitTime").getValueDouble())
	            ))
	        outputDP.addProperty(dafBase.DataProperty.createFloatDataProperty(
	            "expTime", inputDP.findUnique("exposureTime").getValueDouble()))
	        outputDP.addProperty(dafBase.DataProperty.createFloatDataProperty(
	            "airmass", inputDP.findUnique("airMass").getValueDouble()))
	
	        self.activeClipboard.put("rawFPAExposure", outputDP)
	
	        # Rely on postprocess() to put the activeClipboard on the outputQueue.
	       
	    def _lookup(self, filterName):
	        dbLocation = dafPers.LogicalLocation( \
	                "mysql://lsst10.ncsa.uiuc.edu:3306/test")
	        filterDB = lsst.afw.image.Filter(dbLocation, filterName)
	        filterId = filterDB.getId()
                return filterId
