# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
Shutter motion profile storage class
"""

__all__ = ["ShutterMotionProfile", "ShutterMotionProfileFull"]

from astropy.table import Table
from scipy.optimize import newton
import numpy as np

from lsst.ip.isr import IsrCalib


class ShutterMotionProfile(IsrCalib):
    """Shutter motion profile measurements.

    Parameters
    ----------
    log : `logging.Logger`, optional
        Log to write messages to. If `None` a default logger will be used.
    **kwargs :
        Additional parameters.
    """

    _OBSTYPE = "shutterMotionProfile"
    _SCHEMA = "ShutterMotionProfile"
    _VERSION = 1.0

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Quantities that come from `encodeSamples`
        self.time_tai = []
        self.time_mjd = []
        self.position = []
        self.hall_time_tai = []
        self.hall_time_mjd = []
        self.hall_position = []
        self.hall_sensorId = []
        self.hall_isOn = []
        self.fit_name = []
        self.fit_start_time = []
        self.fit_pivot1 = []
        self.fit_pivot2 = []
        self.fit_jerk0 = []
        self.fit_jerk1 = []
        self.fit_jerk2 = []

        self.requiredAttributes.update(["time_tai", "time_mjd", "position",
                                        "hall_time_tai", "hall_time_mjd", "hall_position",
                                        "hall_sensorId", "hall_isOn",
                                        "fit_name", "fit_start_time", "fit_pivot1",
                                        "fit_pivot2", "fit_jerk0", "fit_jerk1", "fit_jerk2",
                                        ])

    def calculateMidpoint(self, modelName="hallSensorFit", skipPosition=False):
        """Calculate time of midpoint of travel for this profile.

        Derived from Shuang Liang's CTN-002 (https://ctn-002.lsst.io).
        Equation numbers listed are from this document.  As the fits
        have already been done, we can ignore the raw position/Hall
        sensor data.

        Parameters
        ----------
        modelName : `str`
            Fit model to use to calculate the midpoint.

        Returns
        -------
        tm_accel : `float`
            The time of the midpoint from the start of motion in
            seconds, as derived from the point where the acceleration
            on the shutter is zero.
        tm_position : `float`
            The time of the midpoint from the start of motion in
            seconds, as derived from the point where the shutter
            position is midway between its starting and ending
            locations.

        Raises
        ------
        RuntimeError
            Raised if the requested ``modelName`` is not found in the
            calibration.
        """
        modelIndex = -1
        for idx, name in enumerate(self.fit_name):
            if name == modelName:
                modelIndex = idx
        if modelIndex == -1:
            raise RuntimeError(f"Unknown model {modelName} requested.")

        # Alias to follow technote
        t0 = self.fit_start_time[modelIndex]
        t1 = self.fit_pivot1[modelIndex]
        t2 = self.fit_pivot2[modelIndex]

        # Equation (3.1)
        j0 = self.fit_jerk0[modelIndex]
        j1 = self.fit_jerk1[modelIndex]

        # Equation (3.2)
        a1 = j0*t1

        # Equation (3.4)
        A1 = a1 - j1*t1

        # First estimate of midpoint, where acceleration is zero.
        # a = 0 = A1 + j1*t  (Equation 5.1)
        def acc(t):
            return A1 + j1 * t

        try:
            tm_accel = newton(acc, 0.5*(t2 + t1))
        except Exception as e:
            self.log.warn(f"Midpoint calculation (from acceleration) failed to converge: {e}")
            tm_accel = np.nan

        if skipPosition:
            tm_position = np.nan
        else:
            # Second estimate of midpoint, when s is halfway betweeen
            # start and final position.  Equation (5.2).
            V1 = t1**2 * (j0 - j1)/2. - t1*A1
            S1 = t1**3 * (j0 - j1)/6. - t1**2 * A1/2. - t1*V1
            Smid = 0.5*(self.metadata["startPosition"] + self.metadata["endPosition"])

            def pos(t):
                return j1*(t**3)/6. + A1*(t**2)/2. + V1*t + S1 - Smid

            try:
                tm_position = newton(pos, tm_accel)
            except Exception as e:
                self.log.warn(f"Midpoint calculation (from position) failed to converge: {e}")
                tm_position = np.nan

        # Restore t0 so these can be compared to raw timestamps.
        return tm_accel + t0, tm_position + t0

    @classmethod
    def fromExposure(cls, exposure, direction="open"):
        """Construct a ShutterMotionProfile from an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposuref`
            Exposure to read header information from.
        direction : `str`, optional
            Direction of shutter to construcxt.  Should be one of
            "open" or "close".

        Returns
        -------
        calib : `lsst.ip.isr.ShutterMotionProfile`
            Constructed profile.

        """
        if direction not in ('open', 'close'):
            raise ValueError(f"Unknown shutter direction {direction} suppled.")

        keywords = [f'SHUTTER {direction.upper()} STARTTIME TAI ISOT',
                    f'SHUTTER {direction.upper()} STARTTIME TAI MJD',
                    f'SHUTTER {direction.upper()} SIDE',
                    f'SHUTTER {direction.upper()} MODEL',
                    f'SHUTTER {direction.upper()} HALLSENSORFIT MODELSTARTTIME',
                    f'SHUTTER {direction.upper()} HALLSENSORFIT PIVOTPOINT1',
                    f'SHUTTER {direction.upper()} HALLSENSORFIT PIVOTPOINT2',
                    f'SHUTTER {direction.upper()} HALLSENSORFIT JERK0',
                    f'SHUTTER {direction.upper()} HALLSENSORFIT JERK1',
                    f'SHUTTER {direction.upper()} HALLSENSORFIT JERK2']

        for kw in keywords:
            if kw not in exposure.metadata:
                raise RuntimeError(f"Expected header keyword not found: {kw}.")

        calib = cls()
        calib.time_tai.append(exposure.metadata[f'SHUTTER {direction.upper()} STARTTIME TAI ISOT'])
        calib.time_mjd.append(exposure.metadata[f'SHUTTER {direction.upper()} STARTTIME TAI MJD'])
        calib.metadata['SIDE'] = exposure.metadata[f'SHUTTER {direction.upper()} SIDE']

        calib.fit_name.append("hallSensorFit")
        # calib.metadata['']
        calib.fit_start_time.append(
            exposure.metadata[f'SHUTTER {direction.upper()} HALLSENSORFIT MODELSTARTTIME']
        )
        calib.fit_pivot1.append(exposure.metadata[f'SHUTTER {direction.upper()} HALLSENSORFIT PIVOTPOINT1'])
        calib.fit_pivot2.append(exposure.metadata[f'SHUTTER {direction.upper()} HALLSENSORFIT PIVOTPOINT2'])
        calib.fit_jerk0.append(exposure.metadata[f'SHUTTER {direction.upper()} HALLSENSORFIT JERK0'])
        calib.fit_jerk1.append(exposure.metadata[f'SHUTTER {direction.upper()} HALLSENSORFIT JERK1'])
        calib.fit_jerk2.append(exposure.metadata[f'SHUTTER {direction.upper()} HALLSENSORFIT JERK2'])

        return calib

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a ShutterMotionProfile from a dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.ShutterMotionProfile
            Constructed calibration.

        Raises
        ------
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration type.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary["fileType"]:
            raise RuntimeError(f"Incorrect calibration supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['OBSTYPE']}")
        motionProfile = dictionary.pop("motionProfile")

        encodeSamples = motionProfile.pop("encodeSamples")
        hallTransitions = motionProfile.pop("hallTransitions")
        fitResults = motionProfile.pop("fitResults")

        if "metadata" in dictionary:
            metadata = dictionary.pop("metadata")
            for key, value in metadata.items():
                dictionary[key] = value
        calib.setMetadata(dictionary)

        formatVersion = calib.metadata["version"]

        startTime = motionProfile.pop("startTime")
        if formatVersion == 1.0:
            # Original format.
            motionProfile["startTime_tai"] = startTime["tai"]
            motionProfile["startTime_mjd"] = startTime["mjd"]
        else:
            # Update to clarify all times are in the TAI system.
            motionProfile["startTime_tai"] = startTime["tai"]["isot"]
            motionProfile["startTime_mjd"] = startTime["tai"]["mjd"]

        calib.readEncodeSamples(encodeSamples, formatVersion)
        calib.readHallTransitions(hallTransitions, formatVersion)
        calib.readFitResults(fitResults)

        calib.updateMetadata(**motionProfile)
        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.

        The dictionary should be able to be round-tripped through
        `fromDict`.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()
        formatVersion = self.metadata["version"]

        if formatVersion == 1.0:
            outDict = {
                "fileName": self.metadata["fileName"],
                "fileType": self.metadata["fileType"],
                "metadata": {
                    "CALIBCLS": "lsst.ip.isr.ShutterMotionProfile",
                    "OBSTYPE": self._OBSTYPE,
                },
                "obsId": self.metadata["obsId"],
                "version": self.metadata.get("version", -1),
                "motionProfile": {
                    "startTime": {
                        "tai": self.metadata["startTime_tai"],
                        "mjd": self.metadata["startTime_mjd"],
                    },
                    "startPosition": self.metadata["startPosition"],
                    "targetPosition": self.metadata["targetPosition"],
                    "endPosition": self.metadata["endPosition"],
                    "targetDuration": self.metadata["targetDuration"],
                    "actionDuration": self.metadata["actionDuration"],
                    "side": self.metadata["side"],
                    "isOpen": self.metadata["isOpen"],
                    "encodeSamples": self.writeEncodeSamples(),
                    "hallTransitions": self.writeHallTransitions(),
                    "fitResults": self.writeFitResults(),
                },
            }
        elif formatVersion == 2.0:
            outDict = {
                "fileName": self.metadata["fileName"],
                "fileType": self.metadata["fileType"],
                "metadata": {
                    "CALIBCLS": "lsst.ip.isr.ShutterMotionProfile",
                    "OBSTYPE": self._OBSTYPE,
                },
                "obsId": self.metadata["obsId"],
                "version": self.metadata.get("version", -1),
                "motionProfile": {
                    "startTime": {
                        "tai": {
                            "isot": self.metadata["startTime_tai"],
                            "mjd": self.metadata["startTime_mjd"],
                        },
                    },
                    "startPosition": self.metadata["startPosition"],
                    "targetPosition": self.metadata["targetPosition"],
                    "endPosition": self.metadata["endPosition"],
                    "targetDuration": self.metadata["targetDuration"],
                    "actionDuration": self.metadata["actionDuration"],
                    "side": self.metadata["side"],
                    "isOpen": self.metadata["isOpen"],
                    "encodeSamples": self.writeEncodeSamples(),
                    "hallTransitions": self.writeHallTransitions(),
                    "fitResults": self.writeFitResults(),
                },
            }
        else:
            raise RuntimeError(f"Unknown file version: {formatVersion}")
        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        This method uses the `fromDict` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables to use to construct the crosstalk
            calibration.  For shutter motion profiles, the first table
            contains the samples, the second the Hall transition data,
            and the third the model fits.

        Returns
        -------
        calib : `lsst.ip.isr.ShutterMotionProfile`
            The calibration defined in the tables.
        """
        samples = tableList[0]
        transitions = tableList[1]
        modelFits = tableList[2]

        metadata = samples.meta

        calib = cls()
        calib.time_tai = np.squeeze(samples["TIME_TAI"].data).tolist()
        if hasattr(calib.time_tai[0], "decode"):
            calib.time_tai = [time.decode("utf-8") for time in calib.time_tai]
        calib.time_mjd = np.squeeze(samples["TIME_MJD"].data).tolist()
        calib.position = np.squeeze(samples["POSITION"].data).tolist()

        calib.hall_time_tai = np.squeeze(transitions["HALL_TIME_TAI"].data).tolist()
        if hasattr(calib.hall_time_tai[0], "decode"):
            calib.hall_time_tai = [time.decode("utf-8") for time in calib.hall_time_tai]
        calib.hall_time_mjd = np.squeeze(transitions["HALL_TIME_MJD"].data).tolist()
        calib.hall_position = np.squeeze(transitions["HALL_POSITION"].data).tolist()
        calib.hall_sensorId = np.squeeze(transitions["HALL_SENSORID"].data).tolist()
        calib.hall_isOn = np.squeeze(transitions["HALL_ISON"].data).tolist()

        calib.fit_model = modelFits.meta["FIT_MODEL"]

        calib.fit_name = np.squeeze(modelFits["FIT_NAME"].data).tolist()
        if hasattr(calib.fit_name[0], "decode"):
            calib.fit_name = [fit.decode("utf-8") for fit in calib.fit_name]
        calib.fit_start_time = np.squeeze(modelFits["FIT_START_TIME"].data).tolist()
        calib.fit_pivot1 = np.squeeze(modelFits["FIT_PIVOT1"].data).tolist()
        calib.fit_pivot2 = np.squeeze(modelFits["FIT_PIVOT2"].data).tolist()
        calib.fit_jerk0 = np.squeeze(modelFits["FIT_JERK0"].data).tolist()
        calib.fit_jerk1 = np.squeeze(modelFits["FIT_JERK1"].data).tolist()
        calib.fit_jerk2 = np.squeeze(modelFits["FIT_JERK2"].data).tolist()

        if "OBSTYPE" not in metadata:
            metadata["OBSTYPE"] = cls._OBSTYPE

        # This translation is needed to support correct
        # round-tripping.  It's not an elegant solution.
        for key in ("fileName", "fileType", "obsId", "version", "side", "isOpen"):
            if key.upper() in metadata:
                value = metadata.pop(key.upper())
                metadata[key] = value
        for key in ("CALIB_ID", "DETECTOR", "DET_NAME", "DET_SER", "FILTER", "INSTRUME",
                    "RAFTNAME", "SEQCKSUM", "SEQFILE", "SEQNAME", "SLOTNAME"):
            if key in metadata:
                if metadata[key] == "":
                    metadata[key] = None

        calib.updateMetadata(**metadata)
        return calib

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables containing the shutter motion profile
            information.
        """
        self.updateMetadata()

        samples = Table(
            {"TIME_TAI": np.array(self.time_tai).tolist(),
             "TIME_MJD": np.array(self.time_mjd).tolist(),
             "POSITION": np.array(self.position).tolist()},
            names=("TIME_TAI", "TIME_MJD", "POSITION"),
            dtype=("U32", "f8", "f8")
        )
        transitions = Table(
            {"HALL_TIME_TAI": np.array(self.hall_time_tai).tolist(),
             "HALL_TIME_MJD": np.array(self.hall_time_mjd).tolist(),
             "HALL_POSITION": np.array(self.hall_position).tolist(),
             "HALL_SENSORID": np.array(self.hall_sensorId).tolist(),
             "HALL_ISON": np.array(self.hall_isOn).tolist()},
            names=("HALL_TIME_TAI", "HALL_TIME_MJD", "HALL_POSITION",
                   "HALL_SENSORID", "HALL_ISON"),
            dtype=("U32", "f8", "f8", "i4", "?")
        )
        modelFits = Table(
            {"FIT_NAME": np.array(self.fit_name).tolist(),
             "FIT_START_TIME": np.array(self.fit_start_time).tolist(),
             "FIT_PIVOT1": np.array(self.fit_pivot1).tolist(),
             "FIT_PIVOT2": np.array(self.fit_pivot2).tolist(),
             "FIT_JERK0": np.array(self.fit_jerk0).tolist(),
             "FIT_JERK1": np.array(self.fit_jerk1).tolist(),
             "FIT_JERK2": np.array(self.fit_jerk2).tolist()},
            names=("FIT_NAME", "FIT_START_TIME", "FIT_PIVOT1", "FIT_PIVOT2",
                   "FIT_JERK0", "FIT_JERK1", "FIT_JERK2"),
            dtype=("U32", "f8", "f8", "f8", "f8", "f8", "f8")
        )
        modelFits.meta["FIT_MODEL"] = self.fit_model

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        samples.meta = outMeta

        return [samples, transitions, modelFits]

    def readEncodeSamples(self, inputSamples, formatVersion):
        """Read a list of input samples into the calibration.

        Parameters
        ----------
        inputSamples : `list` [`dict` [`str` `str`]]
            List of dictionaries of samples.
        formatVersion : `float`
            Version of the file format to read.

        Raises
        ------
        RuntimeError
            Raised if the calibration has already read samples, or if
            the format is not known.
        """
        if len(self.time_tai) != 0:
            raise RuntimeError("Cannot re-read already-read calibration.")

        if formatVersion == 1.0:
            for sample in inputSamples:
                self.time_tai.append(sample["time"]["tai"])
                self.time_mjd.append(sample["time"]["mjd"])
                self.position.append(sample["position"])
        elif formatVersion == 2.0:
            for sample in inputSamples:
                self.time_tai.append(sample["tai"]["isot"])
                self.time_mjd.append(sample["tai"]["mjd"])
                self.position.append(sample["position"])
        else:
            raise RuntimeError(f"Unknown file version: {formatVersion}")

    def writeEncodeSamples(self):
        """Return list of samples as dictionaries.

        Returns
        -------
        inputSamples : `list` [`dict` [`str` `str`]]
            List of dictionaries of samples.

        Raises
        ------
        RuntimeError
            Raised if the calibration has not read samples.
        """
        if len(self.time_tai) == 0:
            raise RuntimeError("Cannot export empty calibration.")

        formatVersion = self.metadata["version"]

        samples = []
        if formatVersion == 1.0:
            for tai, mjd, position in zip(self.time_tai, self.time_mjd, self.position):
                sample = {"time": {"tai": tai, "mjd": mjd},
                          "position": position}
                samples.append(sample)
        elif formatVersion == 2.0:
            for tai, mjd, position in zip(self.time_tai, self.time_mjd, self.position):
                sample = {"tai": {"isot": tai, "mjd": mjd},
                          "position": position}
                samples.append(sample)
        else:
            raise RuntimeError(f"Unknown file version: {formatVersion}")

        return samples

    def readHallTransitions(self, inputTransitions, formatVersion):
        """Read a list of input samples into the calibration.

        Parameters
        ----------
        inputTransitions : `list` [`dict` [`str` `str`]]
            List of dictionaries of transitions.
        formatVersion : `float`
            Version of the file format to read.

        Raises
        ------
        RuntimeError
            Raised if the calibration has already read samples, or if
            the format is not known.
        """
        if len(self.hall_time_tai) != 0:
            raise RuntimeError("Cannot re-read alreday-read calibration.")

        if formatVersion == 1.0:
            for transition in inputTransitions:
                self.hall_time_tai.append(transition["time"]["tai"])
                self.hall_time_mjd.append(transition["time"]["mjd"])
                self.hall_position.append(transition["position"])
                self.hall_sensorId.append(transition["sensorId"])
                self.hall_isOn.append(bool(transition["isOn"]))
        elif formatVersion == 2.0:
            for transition in inputTransitions:
                self.hall_time_tai.append(transition["tai"]["isot"])
                self.hall_time_mjd.append(transition["tai"]["mjd"])
                self.hall_position.append(transition["position"])
                self.hall_sensorId.append(transition["sensorId"])
                self.hall_isOn.append(bool(transition["isOn"]))
        else:
            raise RuntimeError(f"Unknown file version: {formatVersion}")

    def writeHallTransitions(self):
        """Return list of samples as dictionaries.

        Returns
        -------
        inputTransitions : `list` [`dict` [`str` `str`]]
            List of dictionaries of Hall transitions

        Raises
        ------
        RuntimeError
            Raised if the calibration has not read Hall
            transitions.
        """
        if len(self.hall_time_tai) == 0:
            raise RuntimeError("Cannot export empty calibration.")

        formatVersion = self.metadata["version"]
        if formatVersion not in (1.0, 2.0):
            raise RuntimeError(f"Unknown file version: {formatVersion}")
        transitions = []

        for tai, mjd, position, sensorId, isOn in zip(
                self.hall_time_tai,
                self.hall_time_mjd,
                self.hall_position,
                self.hall_sensorId,
                self.hall_isOn):
            if formatVersion == 1.0:
                transition = {"time": {"tai": tai, "mjd": mjd},
                              "position": position,
                              "sensorId": sensorId,
                              "isOn": isOn}
            elif formatVersion == 2.0:
                transition = {"tai": {"isot": tai, "mjd": mjd},
                              "position": position,
                              "sensorId": sensorId,
                              "isOn": isOn}
            transitions.append(transition)
        return transitions

    def readFitResults(self, fitResults):
        """Read a list of fit results into the calibration.

        Parameters
        ----------
        inputTransitions : `list` [`dict` [`str` `str`]]
            List of dictionaries of fit results.

        Raises
        ------
        RuntimeError
            Raised if the calibration has already read fit results.
        """
        if len(self.fit_name) != 0:
            raise RuntimeError("Cannot re-read already-read fit results.")
        self.fit_model = fitResults.pop("Model")

        for fitName, fitModel in fitResults.items():
            if hasattr(fitName, "decode"):
                fitName = fitName.decode("utf-8")
            self.fit_name.append(fitName)
            self.fit_start_time.append(fitModel["ModelStartTime"])
            self.fit_pivot1.append(fitModel["PivotPoint1"])
            self.fit_pivot2.append(fitModel["PivotPoint2"])
            self.fit_jerk0.append(fitModel["Jerk0"])
            self.fit_jerk1.append(fitModel["Jerk1"])
            self.fit_jerk2.append(fitModel["Jerk2"])

    def writeFitResults(self):
        """Return list of samples as dictionaries.

        Returns
        -------
        inputTransitions : `list` [`dict` [`str` `str`]]
            List of dictionaries of Hall transitions

        Raises
        ------
        RuntimeError
            Raised if the calibration has not read Hall
            transitions.
        """
        if len(self.fit_name) == 0:
            raise RuntimeError("Cannot export empty calibration.")

        fitResults = {"Model": self.fit_model}
        for fitName, startTime, pivot1, pivot2, jerk0, jerk1, jerk2 in zip(
                self.fit_name, self.fit_start_time,
                self.fit_pivot1, self.fit_pivot2,
                self.fit_jerk0, self.fit_jerk1, self.fit_jerk2):
            fitResults[fitName] = {"ModelStartTime": startTime,
                                   "PivotPoint1": pivot1,
                                   "PivotPoint2": pivot2,
                                   "Jerk0": jerk0,
                                   "Jerk1": jerk1,
                                   "Jerk2": jerk2}
        return fitResults


class ShutterMotionProfileFull(IsrCalib):
    """Class to hold both open and close profiles, as stored in the
    exposure headers.

    Parameters
    ----------
    log : `logging.Logger`, optional
        Log to write messages to. If `None` a default logger will be used.
    **kwargs :
        Additional parameters.
    """
    _OBSTYPE = "shutterMotionProfileFull"
    _SCHEMA = "ShutterMotionProfileFull"
    _VERSION = 1.0

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.profile_open = None
        self.profile_close = None

        self.requiredAttributes.update(["profile_open", "profile_close"])

    @classmethod
    def fromExposure(cls, exposure):
        """Construct a ShutterMotionProfileFull from an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposuref`
            Exposure to read header information from.
        direction : `str`, optional
            Direction of shutter to construcxt.  Should be one of
            "open" or "close".

        Returns
        -------
        calib : `lsst.ip.isr.ShutterMotionProfile`
            Constructed profile.
        """
        calib = cls()
        calib.profile_open = ShutterMotionProfile.fromExposure(exposure, direction="open")
        calib.profile_close = ShutterMotionProfile.fromExposure(exposure, direction="close")

        return calib

    def calculateMidpoints(self):
        # This is at least a start for downstream calculations.
        midpoint_open, _ = self.profile_open.calculateMidpoint(skipPosition=True)
        midpoint_close, _ = self.profile_close.calculateMidpoint(skipPosition=True)

        return midpoint_open, midpoint_close
