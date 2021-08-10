.. lsst-task-topic:: lsst.meas.algorithms.psfDeterminer.BasePsfDeterminerTask

#####################
BasePsfDeterminerTask
#####################

``BasePsfDeterminerTask`` is the base class for PSF determiners: override ``determinePsf`` to implement your algorithm, and register your new Task with ``psfDeterminerRegistry`` to allow it to be used in a pipeline.


.. _lsst.meas.algorithms.BasePsfDeterminerTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.BasePsfDeterminerTask

.. _lsst.meas.algorithms.BasePsfDeterminerTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.BasePsfDeterminerTask

.. _lsst.meas.algorithms.BasePsfDeterminerTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.BasePsfDeterminerTask
