.. lsst-task-topic:: lsst.meas.algorithms.MeasureApCorrTask

#################
MeasureApCorrTask
#################

``MeasureApCorrTask`` measures aperture correction for the flux fields returned by `lsst.meas.base.getApCorrNameSet`.

.. _lsst.meas.algorithms.MeasureApCorrTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.MeasureApCorrTask

.. _lsst.meas.algorithms.MeasureApCorrTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.MeasureApCorrTask

.. _lsst.meas.algorithms.MeasureApCorrTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.MeasureApCorrTask

.. _lsst.meas.algorithms.MeasureApCorrTask-debug:

Debugging
=========

MeasureApCorrTask has a debug dictionary containing two boolean keys:

display
  display = lsstDebug.Info(__name__).display

doPause
  doPause = lsstDebug.Info(__name__).doPause
