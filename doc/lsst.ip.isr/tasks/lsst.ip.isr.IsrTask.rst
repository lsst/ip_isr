.. lsst-task-topic:: lsst.ip.isr.IsrTask

#######
IsrTask
#######

``IsrTask`` performs instrument signal removal (ISR) from raw exposure data.  There are many configuration options available to enable and disable the many stages of reduction.  This allows for the processing of raw data into full ``postISRCCD`` files for further processing, as well as subsets of that processing to assist in the construction of calibration products.

.. _lsst.ip.isr.IsrTask-processing-summary:

Processing summary
==================

``IsrTask`` runs these operations:

#. Overscan correction using `lsst.ip.isr.OverscanTask`.
#. CCD assembly using `lsst.ip.isr.AssembleCcdTask`.
#. Bias correction using `lsst.ip.isr.isrFunctions.biasCorrection`.
#. Variance plane construction.
#. Linearity correction.
#. Crosstalk correction using `lsst.ip.isr.CrosstalkTask`.
#. Defect masking, NaN masking, saturation trail widening, and any camera specific masking.
#. Brighter-fatter correction using `lsst.ip.isr.isrFunction.brighterFatterCorrection`.
#. Dark correction using `lsst.ip.isr.isrFunctions.darkCorrection`.
#. Fringe correction using `lsst.ip.isr.FringeTask`.
#. Straylight correction using `lsst.ip.isr.StrayLightTask`.
#. Flat correction using `lsst.ip.isr.isrFunctions.flatCorrection`, or gain scaling using `lsst.ip.isr.isrFunctions.applyGains`.
#. Vignette polygon construction and masking using `lsst.ip.isr.VignetteTask`.
#. Attaching transmission curves using `lsst.ip.isr.isrFunctions.attachTransmissionCurve`.
#. Illumination correction using `lsst.ip.isr.isrFunctions.illuminationCorrection`.
#. Interpolation over masked pixels using `lsst.ip.isr.isrFunctions.interpolateFromMask`.
#. Amp-to-amp offset correction using `lsst.ip.isr.AmpOffsetTask`.

.. _lsst.ip.isr.IsrTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.IsrTask

.. _lsst.ip.isr.IsrTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.IsrTask

.. _lsst.ip.isr.IsrTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.IsrTask

.. _lsst.ip.isr.IsrTask-debug:

Debugging
=========

Debug break points exist after a number of the major steps, allowing the exposure processed to that level to be displayed.  Each of these break points are named after the most recent step completed, and include ``doBias``, ``doCrosstalk``, ``doAssembleCcd``, ``doBrighterFatter``, ``doDark``, ``doFringe``, ``doStrayLight``, ``doFlat``, and ``postISRCCD``.  
