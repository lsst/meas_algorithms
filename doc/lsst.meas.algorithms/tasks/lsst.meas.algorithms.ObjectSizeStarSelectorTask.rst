.. lsst-task-topic:: lsst.meas.algorithms.ObjectSizeStarSelectorTask

##########################
ObjectSizeStarSelectorTask
##########################

``ObjectSizeStarSelectorTask`` looks for a cluster of small objects in a size-magnitude plot.

.. _lsst.meas.algorithms.ObjectSizeStarSelectorTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.ObjectSizeStarSelectorTask

.. _lsst.meas.algorithms.ObjectSizeStarSelectorTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.ObjectSizeStarSelectorTask

.. _lsst.meas.algorithms.ObjectSizeStarSelectorTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.ObjectSizeStarSelectorTask

.. _lsst.meas.algorithms.ObjectSizeStarSelectorTask-debug:

Debugging
=========

ObjectSizeStarSelectorTask has a debug dictionary with the following boolean keys:

display
  display = lsstDebug.Info(__name__).display

displayExpsoure
  displayExpsoure = lsstDebug.Info(__name__).displayExposure

plotMagSize
  plotMagSize = lsstDebug.Info(__name__).plotMagSize

dumpData
  dumpData = lsstDebug.Info(__name__).dumpData


More information can be found in `lsstDebug`.
