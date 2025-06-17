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

__all__ = ["ShutterMotionProfile"]

from astropy.table import Table
from scipy.optimize import newton

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

    _OBSTYPE = 'shutterMotionProfile'
    _SCHEMA = 'ShutterMotionProfile'
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

        self.requiredAttributes.update(['time_tai', 'time_mjd', 'position',
                                        'hall_time_tai', 'hall_time_mjd', 'hall_position',
                                        'hall_sensorId', 'hall_isOn',
                                        'fit_name', 'fit_start_time', 'fit_pivot1',
                                        'fit_pivot2', 'fit_jerk0', 'fit_jerk1', 'fit_jerk2',
                                        ])

    def calculateMidpoint(self, modelName="hallSensorFit"):
        """Calculate time of midpoint of travel for this profile.

        Derived from Shuang Liang's CTN-002.  Equation numbers listed
        are from this document.  As the fits have already been done,
        we can ignore the raw position/Hall sensor data.

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
        # I think these are all "six-parameter" models.
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

        tm_accel = newton(acc, 0.5*(t2 + t1))

        # Second estimate of midpoint, when s is halfway betweeen
        # start and final position.  Equation (5.2).
        V1 = t1**2 * (j0 - j1)/2. - t1*A1
        S1 = t1**3 * (j0 - j1)/6. - t1**2 * A1/2. - t1*V1
        Smid = 0.5*(self.metadata['startPosition'] + self.metadata['endPosition'])

        def pos(t):
            return j1*(t**3)/6. + A1*(t**2)/2. + V1*t + S1 - Smid

        tm_position = newton(pos, tm_accel)

        # Restore t0 so these can be compared to raw timestamps.
        return tm_accel + t0, tm_position + t0

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

        if calib._OBSTYPE != dictionary['fileType']:
            raise RuntimeError(f"Incorrect photodiode supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['OBSTYPE']}")
        motionProfile = dictionary.pop('motionProfile')

        encodeSamples = motionProfile.pop("encodeSamples")
        calib.readEncodeSamples(encodeSamples)

        hallTransitions = motionProfile.pop("hallTransitions")
        calib.readHallTransitions(hallTransitions)

        fitResults = motionProfile.pop("fitResults")
        calib.readFitResults(fitResults)

        if 'metadata' in dictionary:
            metadata = dictionary.pop('metadata')
            for key, value in metadata.items():
                dictionary[key] = value
        calib.setMetadata(dictionary)

        startTime = motionProfile.pop("startTime")
        motionProfile["startTime_tai"] = startTime["tai"]
        motionProfile["startTime_mjd"] = startTime["mjd"]

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

        outDict = {"fileName": self.metadata["fileName"],
                   "fileType": self.metadata["fileType"],
                   "metadata": {
                       "CALIBCLS": "lsst.ip.isr.ShutterMotionProfile",
                       "OBSTYPE": self._OBSTYPE, },
                   "obsId": self.metadata["obsId"],
                   "version": self.metadata.get("version", -1),
                   "motionProfile": {
                       "startTime": {
                           "tai": self.metadata["startTime_tai"],
                           "mjd": self.metadata["startTime_mjd"], },
                       "startPosition": self.metadata["startPosition"],
                       "targetPosition": self.metadata["targetPosition"],
                       "endPosition": self.metadata["endPosition"],
                       "targetDuration": self.metadata["targetDuration"],
                       "actionDuration": self.metadata["actionDuration"],
                       "side": self.metadata["side"],
                       "isOpen": self.metadata["isOpen"],
                       "encodeSamples": self.writeEncodeSamples(),
                       "hallTransitions": self.writeHallTransitions(),
                       "fitResults": self.writeFitResults(), }, }

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
        calib.time_tai = samples["TIME_TAI"]
        calib.time_mjd = samples["TIME_MJD"]
        calib.position = samples["POSITION"]

        calib.hall_time_tai = transitions["HALL_TIME_TAI"]
        calib.hall_time_mjd = transitions["HALL_TIME_MJD"]
        calib.hall_position = transitions["HALL_POSITION"]
        calib.hall_sensorId = transitions["HALL_SENSORID"]
        calib.hall_isOn = transitions["HALL_ISON"]

        calib.fit_model = modelFits.meta["FIT_MODEL"]
        calib.fit_name = modelFits["FIT_NAME"]
        calib.fit_start_time = modelFits["FIT_START_TIME"]
        calib.fit_pivot1 = modelFits["FIT_PIVOT1"]
        calib.fit_pivot2 = modelFits["FIT_PIVOT2"]
        calib.fit_jerk0 = modelFits["FIT_JERK0"]
        calib.fit_jerk1 = modelFits["FIT_JERK1"]
        calib.fit_jerk2 = modelFits["FIT_JERK2"]

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

        samples = Table([{"TIME_TAI": self.time_tai,
                          "TIME_MJD": self.time_mjd,
                          "POSITION": self.position}])
        transitions = Table([{"HALL_TIME_TAI": self.hall_time_tai,
                              "HALL_TIME_MJD": self.hall_time_mjd,
                              "HALL_POSITION": self.hall_position,
                              "HALL_SENSORID": self.hall_sensorId,
                              "HALL_ISON": self.hall_isOn}])
        modelFits = Table([{"FIT_NAME": self.fit_name,
                            "FIT_START_TIME": self.fit_start_time,
                            "FIT_PIVOT1": self.fit_pivot1,
                            "FIT_PIVOT2": self.fit_pivot2,
                            "FIT_JERK0": self.fit_jerk0,
                            "FIT_JERK1": self.fit_jerk1,
                            "FIT_JERK2": self.fit_jerk2}])
        modelFits.meta["FIT_MODEL"] = self.fit_model

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        samples.meta = outMeta

        return [samples, transitions, modelFits]

    def readEncodeSamples(self, inputSamples):
        """Read a list of input samples into the calibration.

        Parameters
        ----------
        inputSamples : `list` [`dict` [`str` `str`]]
            List of dictionaries of samples.

        Raises
        ------
        RuntimeError
            Raised if the calibration has already read samples.
        """
        if len(self.time_tai) != 0:
            raise RuntimeError("Cannot re-read already-read calibration.")

        for sample in inputSamples:
            self.time_tai.append(sample["time"]["tai"])
            self.time_mjd.append(sample["time"]["mjd"])
            self.position.append(sample["position"])

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

        samples = []
        for tai, mjd, position in zip(self.time_tai, self.time_mjd, self.position):
            sample = {"time": {"tai": tai, "mjd": mjd},
                      "position": position}
            samples.append(sample)
        return samples

    def readHallTransitions(self, inputTransitions):
        """Read a list of input samples into the calibration.

        Parameters
        ----------
        inputTransitions : `list` [`dict` [`str` `str`]]
            List of dictionaries of transitions.

        Raises
        ------
        RuntimeError
            Raised if the calibration has already read samples.
        """
        if len(self.hall_time_tai) != 0:
            raise RuntimeError("Cannot re-read alreday-read calibration.")

        for transition in inputTransitions:
            self.hall_time_tai.append(transition["time"]["tai"])
            self.hall_time_mjd.append(transition["time"]["mjd"])
            self.hall_position.append(transition["position"])
            self.hall_sensorId.append(transition["sensorId"])
            self.hall_isOn.append(bool(transition["isOn"]))

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

        transitions = []
        for tai, mjd, position, sensorId, isOn in zip(
                self.hall_time_tai,
                self.hall_time_mjd,
                self.hall_position,
                self.hall_sensorId,
                self.hall_isOn):
            transition = {"time": {"tai": tai, "mjd": mjd},
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
