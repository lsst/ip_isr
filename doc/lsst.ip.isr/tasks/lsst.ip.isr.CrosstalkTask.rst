.. lsst-task-topic:: lsst.ip.isr.CrosstalkTask

#############
CrosstalkTask
#############

``CrosstalkTask`` corrects the effects of crosstalk by subtracting scaled copies of the source amplifier images from the target amplifiers.  This task uses the `lsst.ip.isr.CrosstalkCalib` to hold the crosstalk coefficients, as well as a number of the utility functions to ensure that the source and template amplifier images are aligned such that the readout corner is in the same location for both.  Both intra- and inter- detector crosstalk signals can be corrected, although the latter that requires all the exposures for all potential source detectors be supplied.

.. _lsst.ip.isr.CrosstalkTask-processing-summary:

Processing summary
==================

``CrosstalkTask`` runs these operations:

#. Iterate over all potential source amplifier and target amplifiers, identifying which have significant crosstalk coefficients.
#. Flip amplifier images to match readout corners, and subtract the scaled version from the target amplifier.
#. Repeat this procedure for any inter-detector crosstalk signals that have supplied exposures.

.. _lsst.ip.isr.CrosstalkTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.CrosstalkTask

.. _lsst.ip.isr.CrosstalkTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.CrosstalkTask

.. _lsst.ip.isr.CrosstalkTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.CrosstalkTask

.. _lsst.ip.isr.CrosstalkTask-debug:
