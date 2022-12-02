.. lsst-task-topic:: lsst.ip.isr.isrStatistics.IsrStatisticsTask

#################
IsrStatisticsTask
#################

``IsrStatisticsTask`` provides a way to measure some set of properties on the final ISR corrected exposure.  This is intended to allow these properties to be calculated without the need to add multiple additional tasks for each new measurements.

.. _lsst.ip.isr.isrStatistics.IsrStatisticsTask-processing-summary:

Processing summary
==================

``IsrStatisticsTask`` has a run method that can be simply extended for future statistics.  Currently, only the measurement of charge-transfer ineffiency metrics is being calculated.

.. _lsst.ip.isr.isrStatistics.IsrStatisticsTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.isr.isrStatistics.IsrStatisticsTask

.. _lsst.ip.isr.isrStatistics.IsrStatisticsTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.isr.isrStatistics.IsrStatisticsTask

.. _lsst.ip.isr.isrStatistics.IsrStatisticsTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.isr.isrStatistics.IsrStatisticsTask
