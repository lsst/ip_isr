.. lsst-task-topic:: lsst.ip.isr.AssembleCcdTask

###############
AssembleCcdTask
###############

``AssembleCcdTask`` assembles sections of an image into a larger mosaic.
The sub-sections are typically amplifier sections and are to be assembled into a detector size pixel grid.
The assembly is driven by the entries in the raw amp information.
The task can be configured to return a detector image with non-data (e.g. overscan) pixels included.
The task can also renormalize the pixel values to a nominal gain of 1.
The task also removes exposure metadata that has context in raw amps, but not in trimmed detectors (e.g. 'BIASSEC').

.. _lsst.ip.isr.AssembleCcdTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``AssembleCcdTask`` does not have a ``run`` method.
Instead, its main methods are `~lsst.ip.isr.AssembleCcdTask.assembleCcd` and `~lsst.ip.isr.AssembleCcdTask.postprocessExposure`.

.. ``ExampleTask`` runs this sequence of operations:

.. #. Runs this thing. (FIXME)

.. #. Processes processes that intermediate result. (FIXME)

.. #. Stores those results in this last step. (FIXME)

.. _lsst.ip.isr.AssembleCcdTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.AssembleCcdTask

.. _lsst.ip.isr.AssembleCcdTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.AssembleCcdTask

.. _lsst.ip.isr.AssembleCcdTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.AssembleCcdTask

.. _lsst.ip.isr.AssembleCcdTask-debug:

Debugging
=========

The available debug variables in AssembleCcdTask are:

``display``
    A dictionary containing debug point names as keys with frame number as value. Valid keys are:

    assembledExposure
        display assembled exposure
