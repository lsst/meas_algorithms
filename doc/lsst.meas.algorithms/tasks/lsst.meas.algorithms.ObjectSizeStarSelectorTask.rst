.. lsst-task-topic:: lsst.meas.algorithms.ObjectSizeStarSelectorTask

##########################
ObjectSizeStarSelectorTask
##########################

``ObjectSizeStarSelectorTask`` selects likely stars by looking for a cluster of small objects in a size-magnitude plot.

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
 bool; if True display debug information

displayExpsoure
  bool; if True display the exposure and spatial cells

plotMagSize
  bool; if True display the magnitude-size relation using matplotlib

dumpData
  bool; if True dump data to a pickle file


More information can be found in `lsstDebug`.
