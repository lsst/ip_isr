.. lsst-task-topic:: lsst.ip.isr.vignette.VignetteTask

############
VignetteTask
############

``VignetteTask`` determines the intersection of an exposure with a circular vignetted model, optionally masking the pixels that fall outside this model.

.. _lsst.ip.isr.vignette.VignetteTask-processing-summary:

Processing summary
==================

``VignetteTask`` runs these operations:

#. Constructs a circular vignette model from configuration parameters.
#. Sets the exposure valid polygon using this model.
#. Optionally masks the pixels in the exposure that in the vignetted region.

.. _lsst.ip.isr.vignette.VignetteTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.vignette.VignetteTask

.. _lsst.ip.isr.vignette.VignetteTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.vignette.VignetteTask

.. _lsst.ip.isr.vignette.VignetteTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.vignette.VignetteTask
