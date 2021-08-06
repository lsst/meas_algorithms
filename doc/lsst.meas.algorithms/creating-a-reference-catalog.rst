.. _creating-a-reference-catalog:

#########################################
How to generate an LSST reference catalog
#########################################

The LSST Data Management Science Pipeline uses external reference catalogs to fit astrometric and photometric calibration models.
In order to use these external catalogs with our software, we have to convert them into a common format.
This page describes how to "convert" and "ingest" an external catalog for use as a reference catalog for LSST.

The process for generating an LSST-style HTM indexed reference catalog is similar to that of running other LSST Tasks: write an appropriate ``Config`` and run :lsst-task:`~lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask`
The differences are in how you prepare the input data, how you go about running the ``Task``, and what you do with the final output.

Ingesting a large reference catalog can be a slow process.
Ingesting all of Gaia DR2 took a weekend on a high-performance workstation running with 8 parallel processes, for example.

This page uses `Gaia DR2`_ as an example.

.. note::

    If you have an already existing converted reference catalog on disk (for example, the PS1 or Gaia DR2 catalogs that were used in gen2 butlers), you can skip the first few steps, and just :ref:`register and ingest the files <lsst.meas.algorithms-refcat-ingest>` directly, after creating a suitable filaname->htm7-index astropy table file.

.. _Gaia DR2: https://www.cosmos.esa.int/web/gaia/dr2

1. Gathering data
=================

:lsst-task:`~lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask` reads reference catalog data from one or more text or FITS files representing an external catalog (e.g. :file:`GaiaSource*.csv.gz`).
In order to convert these files, you must have a copy of them on a local disk.
Network storage (such as NFS and GPFS) are not recommended for this work, due to performance issues involving tens of thousands of small files.
Ensure that you have sufficient storage capacity.
For example, the GaiaSource DR2 files take 550 GB of space, and the converted LSST reference catalog takes another 200 GB.

To write the config for the converter, you will need to have a document that describes the columns in the input data, including their units and any caveats that may apply (for example, Gaia DR2 does not supply magnitude errors).
In this example, we used the `Gaia Source Catalog data model <https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html>`_ document.

If the files are text files of some sort, check that you can read one of them with `astropy.io.ascii.read`, which is what the converter uses to read text files. For example:

.. code-block:: python

    import astropy.io.ascii
    data = astropy.io.ascii.read('gaia_source/GaiaSource_1000172165251650944_1000424567594791808.csv.gz', format='csv')
    print(data)

The default Config assumes that the files are readable with ``format="csv"``; you can change that to a different ``format`` if necessary (see `lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogConfig` for how to configure that file reader).

.. _lsst.meas.algorithms-refcat-config:

2. Write a Config for the conversion
====================================

`~lsst.meas.algorithms.ConvertReferenceCatalogConfig` specifies what fields in the input files get translated to the output data, and how they are converted along the way.
See the :ref:`ConvertReferenceCatalogTask configuration documentation <lsst.meas.algorithms.ConvertReferenceCatalogTask-configs>` docs for the different available options.

This is an example configuration that was used to convert the Gaia DR2 catalog (saved as ``gaia_dr2_config.py`` in the current directory and provided to the next step as a ``configFile``):

.. code-block:: python

    # The name of the output reference catalog dataset.
    config.dataset_config.ref_dataset_name = "gaia_dr2"

    # Gaia has a specialized convert manager class.
    from lsst.meas.algorithms import convertRefcatManager
    config.manager.retarget(convertRefcatManager.ConvertGaiaManager)

    # Ingest the data in parallel with this many processes.
    config.n_processes = 8

    # These define the names of the fields from the gaia_source data model:
    # https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html

    config.id_name = "source_id"
    config.ra_name = "ra"
    config.dec_name = "dec"
    config.ra_err_name = "ra_error"
    config.dec_err_name = "dec_error"

    config.parallax_name = "parallax"
    config.parallax_err_name = "parallax_error"
    config.coord_err_unit = "milliarcsecond"

    config.pm_ra_name = "pmra"
    config.pm_ra_err_name = "pmra_error"
    config.pm_dec_name = "pmdec"
    config.pm_dec_err_name = "pmdec_error"

    config.epoch_name = "ref_epoch"
    config.epoch_format = "jyear"
    config.epoch_scale = "tcb"

    # NOTE: these names have `_flux` appended to them when the output Schema is created,
    # while the Gaia-specific class handles the errors.
    config.mag_column_list = ["phot_g_mean", "phot_bp_mean", "phot_rp_mean"]

    # These fields are brought along unmodified.
    config.extra_col_names = ["astrometric_excess_noise", "phot_variable_flag"]

In order to deal with the way that Gaia released their photometric data, we have subclassed the conversion manager as `~lsst.meas.algorithms.convertRefcatManager.ConvertGaiaManager`.
This class special-cases the calculation of the flux and flux errors from the values in the Gaia DR2 catalog, which cannot be handled via the simple Config system used above.

.. _lsst.meas.algorithms-refcat-convert:

3. Convert the files to the LSST format
=======================================

:doc:`scripts/convertReferenceCatalog` takes three parameters: output path, configuration file, and quoted input file glob.
See the commandline reference for more details about these parameters.

External catalogs may be split across tens of thousands of files: attempting to specify the full list on the command line is likely to be impossible due to limits imposed by the underlying operating system and shell.
You must specify the input file list as a quoted glob expression; the converter will expand it before processing.
In this example, the output will be written to ``gaia-refcat/`` in the current directory.

You must first run ``setup meas_algorithms`` to use the ``convertReferenceCatalog`` script.

.. prompt:: bash

    convertReferenceCatalog gaia-refcat/ gaia_dr2_config.py "/project/shared/data/gaia_dr2/gaia_source/csv/GaiaSource*" &> convert-gaia.log

To test the conversion without processing the full catalog (which can take many hours), specify a glob pattern that only matches a few files.
For example, ``GaiaSource_970*.csv.gz`` will only process 6 of the GaiaSource files.

Monitor the log file in a new terminal with:

.. prompt:: bash

    tail -f convert-gaia.log

Check the log ouput after several hours: ``ConvertReferenceCatalogTask`` reports progress in 1% intervals.

.. _lsst.meas.algorithms-refcat-ingest:

4. Ingest the files into the butler
===================================

When ``convertReferenceCatalog`` has finished, it will print the two commands you need to run to register the new refcat dataset type and ingest your converted output into it.
For the example we are using here, these commands would be:

.. prompt:: bash

    butler register-dataset-type REPO gaia_dr2 SimpleCatalog htm7
    butler ingest-files -t direct REPO gaia_dr2 refcats gaia/filename_to_htm.ecsv

where REPO is the path to the butler repository that you are ingesting the data into.
We use the ``direct`` transfer mode here to keep the files in the directory we converted them to. See ``butler ingest-files -h`` for other options.

These commands should finish in a short amount of time, logging a message about how many files were ingested.
You can query the ``refcats`` collection to see whether your htm shards appear:

.. prompt:: bash

    butler query-datasets --collections refcats REPO

For LSST staff using ``lsst-dev``, see the `Reference catalogs policy <https://developer.lsst.io/services/datasets.html#reference-catalogs>`_ in the Developer Guide for additional policy about adding reference catalogs to the common repo.
