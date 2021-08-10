.. lsst-task-topic:: lsst.meas.algorithms.dynamicDetection.DynamicDetectionTask

####################
DynamicDetectionTask
####################

``DynamicDetectionTask`` detects sources on an image with a dynamic threshold.

.. _lsst.meas.algorithms.DynamicDetectionTask-summary:

Processing summary
==================

``DynamicDetectionTask`` runs this sequence of operations:

#. Detects sources using a lower threshold than normal in order to identify good sky regions.

#. Performs forced PSF photometry on identified sky regions.

#. Sets threshold so that the standard deviation of measurements matches median estimated error using PSF flux measurements and estimated errors.

.. _lsst.meas.algorithms.DynamicDetectionTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.DynamicDetectionTask

.. _lsst.meas.algorithms.DynamicDetectionTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.DynamicDetectionTask

.. _lsst.meas.algorithms.DynamicDetectionTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.DynamicDetectionTask

