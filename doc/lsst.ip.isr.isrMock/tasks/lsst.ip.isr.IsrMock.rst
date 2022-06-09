.. lsst-task-topic:: lsst.ip.isr.IsrMock

#######
IsrMock
#######

``IsrMock`` creates a variety of simulated images and calibration products, for use in testing.

.. _lsst.ip.isr.IsrMock-processing-summary:

Processing summary
==================

``IsrMock`` can generate simulated "raw" images (and the associated calibration products needed to correct them) that may contain any or all of the following:

#. A simulated sky with Poissonian noise.
#. A single Gaussian source.
#. A fringe signal.
#. Flat field correction.
#. Dark current and noise.
#. Bias offset and noise.
#. An overscan with a simple gradient functional form.
#. Crosstalk between amps.

The output exposure can be returned as a trimmed (overscan and prescan removed) single image, an untrimmed single image, or a dictionary containing each amplifier as a separate image indexed by the amplifier name.  This allows testing of the various methods of assembly.

.. _lsst.ip.isr.IsrMock-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.IsrMock

.. _lsst.ip.isr.IsrMock-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.IsrMock

.. _lsst.ip.isr.IsrMock-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.IsrMock

.. _lsst.ip.isr.IsrMock-debug:
