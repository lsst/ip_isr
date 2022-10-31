.. lsst-task-topic:: lsst.ip.isr.deferredCharge.DeferredChargeTask

##################
DeferredChargeTask
##################

``DeferredChargeTask`` corrects the input exposure using a per-amplifier model of the charge-transfer inefficiency.

.. _lsst.ip.isr.deferredCharge.DeferredChargeTask-processing-summary:

Processing summary
==================

``DeferredChargeTask`` runs these operations:

#. Transforms the input amplifier data to have a common orientation via a series of flips.
#. Corrects the image data for local offsets.
#. Corrects the image data for serial traps.

.. _lsst.ip.isr.deferredCharge.DeferredChargeTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.deferredCharge.DeferredChargeTask

.. _lsst.ip.isr.deferredCharge.DeferredChargeTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.deferredCharge.DeferredChargeTask

.. _lsst.ip.isr.deferredCharge.DeferredChargeTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.deferredCharge.DeferredChargeTask

.. _lsst.ip.isr.deferredCharge.DeferredChargeTask-debug:

Debugging
=========

Debugging hooks exist to display the exposure and the locations of the random positions, as well as to display the fit between the science and fringe measurements.  These are disabled in the code by default.
