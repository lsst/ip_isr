.. lsst-task-topic:: lsst.ip.isr.ampOffset.AmpOffsetTask

#############
AmpOffsetTask
#############

``AmpOffsetTask`` is designed to measure background levels at amplifier boundaries, using those values to calculate the optimal offsets needed to align the flux levels to remove remnant shifts between amplifiers that are not corrected by bias, dark, and flat correction.

.. _lsst.ip.isr.ampOffset.AmpOffsetTask-processing-summary:

Processing summary
==================

``AmpOffsetTask`` is a stub task that should be subclassed by camera-specific implementations that can correctly map the measurements into amp-to-amp shifts.

.. _lsst.ip.isr.ampOffset.AmpOffsetTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.ampOffset.AmpOffsetTask

.. _lsst.ip.isr.ampOffset.AmpOffsetTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.ampOffset.AmpOffsetTask

.. _lsst.ip.isr.ampOffset.AmpOffsetTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.ampOffset.AmpOffsetTask

