.. lsst-task-topic:: lsst.ip.isr.assembleCcdTask.AssembleCcdTask

###############
AssembleCcdTask
###############

``AssembleCcdTask`` constructs a full detector image from individual segments.  The end result can either be untrimmed (all overscan and prescan regions are retained), or trimmed (these sections are removed, producing an image of only the imaging region).

.. _lsst.ip.isr.assembleCcdTask.AssembleCcdTask-processing-summary:

Processing summary
==================

``AssembleCcdTask`` runs these operations:

#. The exposure detector is identified, and a loop over each amplifier trims and shifts each segment into the correct location of the assembled exposure.
#. The exposure metadata is filtered to remove keywords that have context in the raw amplifier data, but not in the assembled image.
#. The WCS, filter, and visit info are also transferred from the input exposure to the output.


.. _lsst.ip.isr.assembleCcdTask.AssembleCcdTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.assembleCcdTask.AssembleCcdTask

.. _lsst.ip.isr.assembleCcdTask.AssembleCcdTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.assembleCcdTask.AssembleCcdTask

.. _lsst.ip.isr.assembleCcdTask.AssembleCcdTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.assembleCcdTask.AssembleCcdTask

.. _lsst.ip.isr.assembleCcdTask.AssembleCcdTask-debug:

Debugging
=========

If the a debug frame for ``assembledExposure`` is ``True`` for ``lsst.ip.isr.assembleCcdTask``, then the final assembled exposure is displayed.
