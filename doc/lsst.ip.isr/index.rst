.. py:currentmodule:: lsst.ip.isr

.. _lsst.ip.isr:

###########
lsst.ip.isr
###########

The ``lsst.ip.isr`` module provides instrument signature removal (ISR) related tasks.

ISR includes steps such as combining multiple amplifiers into one full CCD image, corrections for overscans, crosstalk, bias and dark frames, and the creation of variance and mask planes.

.. _lsst.ip.isr-using:

Using lsst.ip.isr
=================

.. toctree::
   :maxdepth: 1

``lsst.ip.isr`` is generally used as the initial step of pipeline processing, performing the initial image processing for ``lsst.drp.pipe``, ``lsst.ap.pipe``, and ``lsst.cp.pipe``.  This module also contains the definitions for calibration products that are not just a simple image.  The final major component is a set of tasks designed to simulate raw data and calibration products, to allow for the functions and methods to be properly tested with known inputs.

   calibration-types
   isrmock-reference

.. _lsst.ip.isr-contributing:

Contributing
============

``lsst.ip.isr`` is developed at https://github.com/lsst/ip_isr.
You can find Jira issues for this module under the `ip_isr <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20ip_isr>`_ component.

.. If there are topics related to developing this module (rather than using it), link to this from a toctree placed here.

.. .. toctree::
..    :maxdepth: 1

Task reference
==============

.. _lsst.ip.isr-pipeline-tasks:

Pipeline tasks
--------------

.. lsst-pipelinetasks::
   :root: lsst.ip.isr

.. _lsst.ip.isr-command-line-tasks:

Command-line tasks
------------------

.. lsst-cmdlinetasks::
   :root: lsst.ip.isr

.. _lsst.ip.isr-tasks:

Tasks
-----

.. lsst-tasks::
   :root: lsst.ip.isr
   :toctree: tasks

.. _lsst.ip.isr-pyapi:

Python API reference
====================

.. automodapi:: lsst.ip.isr
   :no-main-docstr:
   :no-inheritance-diagram:

.. automodapi:: lsst.ip.isr.ampOffset
   :no-main-docstr:
   :no-inheritance-diagram:

.. automodapi:: lsst.ip.isr.vignette
   :no-main-docstr:
   :no-inheritance-diagram:
