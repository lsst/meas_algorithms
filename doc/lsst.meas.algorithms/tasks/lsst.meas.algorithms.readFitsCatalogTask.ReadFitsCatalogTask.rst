.. lsst-task-topic:: lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask

###################
ReadFitsCatalogTask
###################

``ReadFitsCatalogTask`` reads an object catalog from a FITS file into a numpy array object suitable for use with `lsst.meas.algorithms.IngestIndexReferenceTask`.

.. _lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask

.. _lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask

.. _lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask

.. _lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask-examples:

Examples
========

A complete example of using ReadFitsCatalogTask:

.. code-block:: python

   from lsst.meas.algorithms.readFitsCatalogTask import ReadFitsCatalogTask
   filePath = "tests/data/testReadFitsCatalog.fits"
   task = ReadFitsCatalogTask()
   catalogArray = task.run(filePath)

The resulting `catalogArray` is a numpy structured array containing fields such as "name", "ra" and "dec"
and a few rows of data. For more complicated cases config parameters allow you to rename columns
and choose which HDU to read.
