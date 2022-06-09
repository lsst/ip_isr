.. lsst-task-topic:: lsst.ip.isr.OverscanCorrectionTask

######################
OverscanCorrectionTask
######################

``OverscanCorrectionTask`` calculates and subtracts the overscan from the input exposure.

.. _lsst.ip.isr.OverscanCorrectionTask-processing-summary:

Processing summary
==================

``OverscanCorrectionTask`` runs these operations:

#. Fits the overscan, either as a single constant value (for the 'MEAN', 'MEANCLIP', and 'MEDIAN' methods) or as a per-row vector (for the 'MEDIAN_PER_ROW', polynomial, and spline methods).
#. Subtracts that model from both the image and overscan arrays,
#. Optionally masks suspect pixels.

.. _lsst.ip.isr.OverscanCorrectionTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.OverscanCorrectionTask

.. _lsst.ip.isr.OverscanCorrectionTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.OverscanCorrectionTask

.. _lsst.ip.isr.OverscanCorrectionTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.OverscanCorrectionTask

.. _lsst.ip.isr.OverscanCorrectionTask-debug:

Debugging
=========

If ``debug.display`` is true for ``lsst.ip.isr.overscan``, then per-amplifier plots of the data values and the overscan model are shown for evaluation.
