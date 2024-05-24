.. lsst-task-topic:: lsst.meas.algorithms.normalizedCalibrationFlux.NormalizedCalibrationFluxTask

#############################
NormalizedCalibrationFluxTask
#############################

``NormalizedCalibrationFluxTask`` takes an unnormalized calibration flux (e.g. ``base_CompensatedTophatFlux_12``) and normalizes it to an overall calibration flux (e.g. ``base_CircularApertureFlux_12_0``).
This task is essentially a specialized wrapper around `lsst.meas.algorithms.MeasureApCorrTask` that can measure and/or apply an aperture correction map into a new field name given by the ``normalized_calibflux_name`` config field.

.. _lsst.meas.algorithms.NormalizedCalibrationFluxTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.MeasureApCorrTask

.. _lsst.meas.algorithms.NormalizedCalibrationFluxTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.NormalizedCalibrationFluxTask

.. _lsst.meas.algorithms.NormalizedCalibrationFluxTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.NormalizedCalibrationFluxTask

.. _lsst.meas.algorithms.NormalizedCalibrationFluxTask-debug:
