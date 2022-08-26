.. lsst-task-topic:: lsst.meas.algorithms.detection.SourceDetectionTask

###################
SourceDetectionTask
###################

``SourceDetectionTask`` Detect positive and negative sources on an exposure
and return a new table.SourceCatalog

.. _lsst.meas.algorithms.SourceDetectionTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.SourceDetectionTask

.. _lsst.meas.algorithms.SourceDetectionTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.SourceDetectionTask

.. _lsst.meas.algorithms.SourceDetectionTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.SourceDetectionTask

.. _lsst.meas.algorithms.SourceDetectionTask-examples:

Examples
========

This code is in `measAlgTasks.py` in the examples directory, and can be run as *e.g.*

.. code-block:: bash

   examples/measAlgTasks.py --doDisplay

The example also runs the SingleFrameMeasurementTask; see meas_algorithms_measurement_Example
for more explanation.

Import the task (there are some other standard imports; read the file if you're confused)

.. code-block:: py

   from lsst.meas.algorithms.detection import SourceDetectionTask

We need to create our task before processing any data as the task constructor can add an
extra column to the schema, but first we need an almost-empty Schema

.. code-block:: py

   schema = afwTable.SourceTable.makeMinimalSchema()

after which we can call the constructor:

.. code-block:: py

  config = SourceDetectionTask.ConfigClass()
  config.thresholdPolarity = "both"
  config.background.isNanSafe = True
  config.thresholdValue = 3
  detectionTask = SourceDetectionTask(config=config, schema=schema)

We're now ready to process the data (we could loop over multiple exposure/catalogues using
the same task objects). First create the output table:

.. code-block:: py

   table = afwTable.SourceTable.make(schema)

And process the image

.. code-block:: py

   result = detectionTask.run(table, exposure)

(You may not be happy that the threshold was set in the config before creating the Task
rather than being set separately for each exposure. You **can** reset it just before calling
the run method if you must, but we should really implement a better solution).

We can then unpack the results:

.. code-block:: py

   sources = result.sources
   print("Found %d sources (%d +ve, %d -ve)" % (len(sources), result.fpSets.numPos,
   result.fpSets.numNeg))

.. _lsst.meas.algorithms.SourceDetectionTask-debug:

Debugging
=========

The :ref:`pipetask run <pipetask-command>` command-line interface
supports a flag ``--debug`` to to import `debug.py` from your ``PYTHONPATH``; see
`lsstDebug` for more about `debug.py` files.

The available variables in `SourceDetectionTask` are:

`display`

        #. If True, display the exposure of afwDisplay.Display's frame 0. Positive detections
           in blue, negative detections in cyan.
        #. If display > 1, display the convolved exposure on frame 1
