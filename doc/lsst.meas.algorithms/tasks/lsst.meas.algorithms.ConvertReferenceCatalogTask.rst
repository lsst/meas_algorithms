.. lsst-task-topic:: lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask

###########################
ConvertReferenceCatalogTask
###########################

``ConvertReferenceCatalogTask`` converts an external catalog for use as an LSST Science Pipelines reference catalog, using a Hierarchical Triangular Mesh (HTM) indexing scheme. The format and layout of the input data is configurable. The output data is a collection of `lsst.afw.table.SimpleCatalog` files identified by their HTM pixel. This task is not available as a command-line task: see :ref:`creating-a-reference-catalog` for how to run the task.

.. _lsst.meas.algorithms.ConvertReferenceCatalogTask-summary:

Processing summary
==================

``ConvertReferenceCatalogTask`` uses Python `multiprocessing` to ingest multiple files in parallel, configured by :lsst-config-field:`~lsst.meas.algorithms.convertReferenceCatalog.IngestIndexedReferenceConfig.n_processes`.
Once it has generated the necessary multiprocessing file locks (one per output file: ~130,000 files for HTM ``depth=7``), it performs the following steps for each input file:

#. Reads the file using the configured :lsst-config-field:`~lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogConfig.file_reader` subtask (default: :lsst-task:`~lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogTask`).

#. Indexes the coordinates in the input data to determine which mesh pixel they go with, and thus which output file they will be written to.

#. Loops over the output pixels in this input file (where N is the number of sources in this pixel):

 #. Acquires the lock for this output file.

 #. Reads an existing output file and appends N new empty rows, or generates a new empty catalog with N rows.

 #. Fills in the empty rows of the catalog with the converted values from the input data.

 #. Writes the output file and releases the file lock.

.. lsst.meas.algorithms.ConvertReferenceCatalogTask-cli:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask

.. _lsst.meas.algorithms.ConvertReferenceCatalogTask-butler:

Butler datasets
===============

``ConvertReferenceCatalogTask`` does not behave in the same manner as most LSST Tasks.
When run directly through the `~lsst.meas.algorithms.ConvertReferenceCatalogTask.run` method, ``ConvertReferenceCatalogTask`` reads input from a collection of non-LSST files, and writes to a specified directory.
Note that configurations for ``ConvertReferenceCatalogTask``, and its subtasks, affect what the output dataset content is.

.. _lsst.meas.algorithms.ConvertReferenceCatalogTask-butler-outputs:

Output datasets
---------------

``ref_cat``
    An LSST-style reference catalog, consisting of one `lsst.afw.table.SimpleCatalog` per HTM pixel.


.. _lsst.meas.algorithms.ConvertReferenceCatalogTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask

.. _lsst.meas.algorithms.ConvertReferenceCatalogTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask

.. _lsst.meas.algorithms.ConvertReferenceCatalogTask-examples:

Examples
========

See :ref:`creating-a-reference-catalog` for a description of how to run the task to ingest the Gaia DR2 catalog.
