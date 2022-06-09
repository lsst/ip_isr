.. lsst-task-topic:: lsst.ip.isr.FringeTask

##########
FringeTask
##########

``FringeTask`` removes the fringing signal that is evident in long wavelength data.

.. _lsst.ip.isr.FringeTask-processing-summary:

Processing summary
==================

``FringeTask`` runs these operations:

#. A background "pedestal" level is optionally subtracted from the input fringes.
#. A large number of random locations on the exposure are identified.
#. At each of those random positions, the local sky level is estimated by taking a heavily (3.0 sigma) clipped median in a box.
#. Those positions are measured similarly on the fringe frame.
#. An iterative clipping least squares process is run to find the best fitting scale factor to minimize

.. math::
y = science - scale * fringe

#. This scaled fringe signal is subtracted from the science image.
#. This process can be repeated if multiple fringe frames exist.


.. _lsst.ip.isr.FringeTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.FringeTask

.. _lsst.ip.isr.FringeTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.FringeTask

.. _lsst.ip.isr.FringeTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.FringeTask

.. _lsst.ip.isr.FringeTask-debug:

Debugging
=========

Debugging hooks exist to display the exposure and the locations of the random positions, as well as to display the fit between the science and fringe measurements.  These are disabled in the code by default.
