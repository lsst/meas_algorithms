.. lsst-task-topic:: lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask

###################
ReadTextCatalogTask
###################

``ReadTextCatalogTask`` reads an object catalog from a UTF-8 encoded text file into a numpy array object suitable for use with `lsst.meas.algorithms.IngestIndexReferenceTask`.

.. _lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask

.. _lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask

.. _lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask

.. _lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask-examples:

Examples
========

Given a file named `table.csv` containing the following:
        ra,     dec,    flux
        5.5,    -45.2,  12453
        19.6,   34.2,   32123

you can read this file with the following code:

.. code-block:: py

   from lsst.meas.algorithms.readTextCatalogTask import ReadTextCatalogTask
   task = ReadTextCatalogTask()
   catalogArray = task.run("table.csv")
    
The resulting `catalogArray` is a numpy structured array containing three fields
("ra", "dec" and "flux") and two rows of data. For more complex cases,
config parameters allow you to specify the names of the columns (instead of using automatic discovery)
and set the number of rows to skip.
