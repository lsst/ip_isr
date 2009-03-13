"""
Manage a "data base" of calibration data

E.g.
import glob
import lsst.ip.isr.calibDatabase as calibDatabase
import lsst.daf.base

>>> calibDatabase.writeCalibValidityPaf(glob.glob("/lsst/images/repository/calib/03Am06/*-*-c00[89]-*[01].fits"),
                                        fd=open("/home/rhl/LSST/ip/isr/pipeline/calibDatabase.paf", "w"),
                                        "/lsst/images/repository/calib")

>>> when = lsst.daf.base.DateTime(2003, 07, 21, 0, 0, 0)
>>> cdb = calibDatabase.CalibDB("/home/rhl/LSST/ip/isr/pipeline/calibDatabase.paf")
>>> print cdb.lookup(when, "bias", "CCD009", 1)

(N.b. CCDs may be named as 11 or "CCD011", amplifiers as 1 or "Amplifier001")
"""

import datetime, os, re, sys
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexExcept

class _CalibData(object):
    """Contain what we know about calibration data"""

    def __init__(self, fileName, version, validFrom, validTo, expTime=0, filter=None):
        self.fileName = fileName
        self.version = version
        self.expTime = expTime
        self.validFrom = validFrom
        self.validTo = validTo
        self.filter = filter

def needExpTime(calibType):
    return calibType in ("bias", "dark")

def needFilter(calibType):
    return calibType in ("flat", "fringe")

def writeCalibValidityPaf(fileNameList, fd=sys.stdout, stripPrefix=None):

    if isinstance(fileNameList, str):
        fileNameList = [fileNameList]

    if stripPrefix:
        stripPrefix = re.sub(r"/*$", "", stripPrefix)

    calibInfo = {}
    for fileName in fileNameList:
        dirName = os.path.split(os.path.dirname(fileName))[-1]
        basename = os.path.basename(fileName)

        name, suffix = os.path.splitext(basename)
        if suffix == ".fits":
            calibType, extra, ccd, amp = name.split("-")
        elif suffix == ".paf":
            calibType, ccd, amp = name.split("-")
        else:
            print >> sys.stderr, "I don't recognize %s's suffix" % fileName
            continue

        expTime, filter = None, None
        if needExpTime(calibType):
            expTime = extra
        elif needFilter(calibType):
            filter = extra
        elif calibType in ("scatter"):
            continue;                   # ignoring scattered light frames
        elif calibType in ("defect"):
            pass
        else:
            print >> sys.stderr, "I don't know what to do with %s" % (fileName)
            continue

        if re.search(r"^c\d+$", ccd):
            ccd = "CCD%03d" % int(ccd[1:])
        else:
            print >> sys.stderr, "Expected cXXXX, saw %s in %s" % (ccd, fileName)
            continue

        if re.search(r"^a\d+$", amp):
            amp = "Amplifier%03d" % int(amp[1:])
        else:
            print >> sys.stderr, "Expected cXXXX, saw %s in %s" % (amp, fileName)
            continue

        if calibType == "defect":
            crunid = dirName            # we'll check that they are the same
            validFrom = "1970/01/01"
            validTo = validFrom
        else:
            meta = afwImage.readMetadata(fileName, 0)

            crunid = meta.get("CRUNID").strip()
            validFrom = meta.get("TVSTART")   # data validity start time
            validTo = meta.get("TVSTOP")     # data validity end time
        #
        # Reformat YYYY/MM/DD dates to ISO
        #
        validFrom = "%4d-%02d-%02dT00:00:00.00Z" % tuple([int(x) for x in validFrom.split("/")])
        validTo = "%4d-%02d-%02dT00:00:00.00Z" %  tuple([int(x) for x in validTo.split("/")])

        if crunid != dirName:
            print >> sys.stderr, "Expected CRUNID = %s, saw %s" % (crunid, dirName)

        if not calibInfo.has_key(ccd):
            calibInfo[ccd] = {}
        if not calibInfo[ccd].has_key(amp):
            calibInfo[ccd][amp] = {}
        if not calibInfo[ccd][amp].has_key(calibType):
            calibInfo[ccd][amp][calibType] = []

        if stripPrefix:
            if os.path.commonprefix([fileName, stripPrefix]) == stripPrefix:
                fileName = fileName[len(stripPrefix) + 1:]

        calibInfo[ccd][amp][calibType].append(_CalibData(fileName, crunid, validFrom, validTo,
                                                         expTime=expTime, filter=filter))
    #
    # Write that out
    #
    print >> fd, "calibrations: {"

    for ccd in sorted(calibInfo.keys()):
        print >> fd, "    %s: {" % ccd

        for amp in sorted(calibInfo[ccd].keys()):
            print >> fd, "        %s: {" % amp
            for calibType in sorted(calibInfo[ccd][amp].keys()):
                for calib in calibInfo[ccd][amp][calibType]:
                    print >> fd, """\
            %s: {
                fileName:  %-25s  # Calibration file
                version:   %-25s  # This calibration's version
                validFrom: %-25s  # Start time of this calibration's validity
                validTo:   %-25s  # End time of this calibration's validity
""" % (calibType, ('"%s"' % calib.fileName), ('"%s"' % calib.version),
       '"%s"' % calib.validFrom,  '"%s"' % calib.validTo),

                    if needExpTime(calibType):
                        print >> fd, "                expTime:   %4s                       # Exposure time" % calib.expTime

                    if needFilter(calibType):
                        print >> fd, "                filter:    \"%s\"                        # Filter" % calib.filter
                    
                    print >> fd, "            }"


            print >> fd, "        }"

        print >> fd, "    }"

    print >> fd, "}"
    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def DateTimeFromIsoStr(str, scale=dafBase.DateTime.TAI):
    """Convert a format of the form 2003-07-20T23:12:19.00Z (n.b. no fractional seconds)"""
    yearStr, timeStr = str.split("T")
    year, month, day = [int(x) for x in yearStr.split("-")]
    timeStr = re.sub(r"(\.0+)?Z$", "", timeStr)

    hr, min, sec = [int(x) for x in timeStr.split(":")]

    return dafBase.DateTime(year, month, day, hr, min, sec, scale)

class CalibDB(object):
    """A class to find the proper calibration files for a given type of calibration"""

    def __init__(self, calibDatabasePaf):
        """Read calibration file in calibDatabasePaf"""

        self.calibDatabasePaf = calibDatabasePaf

        try:
            self.calibPolicy = pexPolicy.Policy(self.calibDatabasePaf)
        except pexExcept.LsstCppException, e:
            raise "Failed to read %s: %s" % (self.calibDatabasePaf, e)

    def lookup(self, lsstDateTime, calibType, CCD="CCD009", amplifier=1, filter=None, expTime=None):
        """Find the  proper calibration given an lsst::daf::data::DateTime, a calib type, a CCD and an amplifier; if appropriate, a filter may also be specified

Calibrations are only valid for a range of times (special case:  if the times are equal, it is
assumed that the files are always valid)

Valid calibTypes are bias, dark, defect, flat, fringe, and linearize
"""

        if isinstance(CCD, int):
            CCD = "CCD%03d" % CCD
        if isinstance(amplifier, int):
            amplifier = "Amplifier%03d" % amplifier

        if calibType not in ("bias", "dark", "defect", "flat", "fringe", "linearize"):
            raise RuntimeError, ("Unknown calibration type: %s" % calibType)
        #
        # A placeholder
        #
        if calibType == "linearize":
            return "linearizationLookupTable.paf"

        if calibType == "bias":
            if expTime:
                raise RuntimeError, ("You may not specify an expTime for a bias: %s" % expTime)
            expTime = 0

        if needExpTime(calibType) and expTime is None:
            raise RuntimeError, ("Please specify an expTime for your %s" % (calibType))

        if needFilter(calibType) and not filter:
            raise RuntimeError, ("Please specify a filter for your %s" % (calibType))

        try:
            for calib in self.calibPolicy.getPolicy("calibrations").getPolicy(CCD).getPolicy(amplifier).getArray(calibType):
                validTo = DateTimeFromIsoStr(calib.get("validTo"))
                validFrom = DateTimeFromIsoStr(calib.get("validFrom"))

                if validFrom.nsecs() == validTo.nsecs() or validFrom.nsecs() <= lsstDateTime.nsecs() < validTo.nsecs():
                    if needExpTime(calibType):
                        if calib.get("expTime") != expTime:
                            continue

                    if needFilter(calibType):
                        if calib.get("filter") != filter:
                            continue

                    fileName = calib.get("fileName")

                    return fileName

        except IndexError, e:
            pass
        except TypeError, e:
            pass
        except pexExcept.LsstCppException, e:
            pass

        ctype = calibType
        if needExpTime(calibType):
            ctype += " %s" % expTime
        if needFilter(calibType):
            ctype += " %s" % filter

        raise RuntimeError, "Unable to locate %s for %s %s for %s" % (ctype, CCD, amplifier,
                                                                      datetime.datetime.fromtimestamp(int(lsstDateTime.nsecs()/1e9)).strftime("%Y-%m-%dT%H:%M:%SZ"))
