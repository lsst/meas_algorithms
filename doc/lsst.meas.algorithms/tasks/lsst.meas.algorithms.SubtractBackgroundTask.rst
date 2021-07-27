.. lsst-task-topic:: lsst.meas.algorithms.SubtractBackgroundTask

######################
SubtractBackgroundTask
######################

``SubtractBackgroundTask`` fits a model of the background of an exposure and subtracts it.

.. _lsst.meas.algorithms.SubtractBackgroundTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.SubtractBackgroundTask

.. _lsst.meas.algorithms.SubtractBackgroundTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.SubtractBackgroundTask

.. _lsst.meas.algorithms.SubtractBackgroundTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.SubtractBackgroundTask

.. _lsst.meas.algorithms.SubtractBackgroundTask-indepth:

In Depth
========

Quantities set in exposure Metadata
-----------------------------------

The `run` method will optionally set the following items of exposure metadata;
the names may be overridden; the defaults are shown:

    - ``BGMEAN``: Mean value of the background
    - ``BGVAR``: Standard deviation of background

.. _lsst.meas.algorithms.SubtractBackgroundTask-debug:

Debugging
=========

SubtractBackgroundTask has a debug dictionary containing three integer keys:

unsubtracted
  ``int``; If >0: `fitBackground` displays the unsubtracted masked image overlaid with the grid of cells used to fit the background in the specified frame.

subtracted
  ``int``; If >0: `run` displays the background-subtracted exposure in the specified frame.

background
  ``int``; If >0: `run` displays the background image in the specified frame.
