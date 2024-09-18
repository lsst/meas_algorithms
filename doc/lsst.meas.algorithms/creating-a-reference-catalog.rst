.. _creating-a-reference-catalog:

#########################################
How to generate an LSST reference catalog
#########################################

The LSST Data Management Science Pipeline uses external reference catalogs to fit astrometric and photometric calibration models.
In order to use these external catalogs with our software, we have to convert them into a common format.
This page describes how to "convert" and "ingest" an external catalog for use as a reference catalog for LSST.

The process for generating an LSST-style HTM indexed reference catalog is similar to that of running other LSST Tasks: write an appropriate ``Config`` and run :lsst-task:`~lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask`
The differences are in how you prepare the input data, how you go about running the ``Task``, and what you do with the final output.

Converting a large reference catalog can be a slow process.
Converting all of Gaia DR2 took a weekend on a high-performance workstation running with 8 parallel processes, for example.

This page uses `Gaia DR2`_ as an example.

.. note::

    If you have an already existing converted reference catalog on disk (for example, the PS1 or Gaia DR2 catalogs that were used in gen2 butlers), you can skip the first few steps, and just :ref:`register and ingest the files <lsst.meas.algorithms-refcat-ingest>` directly, after :ref:`creating a suitable filename to htm7-index astropy lookup table <lsst.meas.algorithms-refcat-existing>`.

.. _Gaia DR2: https://www.cosmos.esa.int/web/gaia/dr2

1. Gathering data
=================

:lsst-task:`~lsst.meas.algorithms.convertReferenceCatalog.ConvertReferenceCatalogTask` reads reference catalog data from one or more text or FITS files representing an external catalog (e.g. :file:`GaiaSource*.csv.gz`).
In order to convert these files, you must have a copy of them on a local disk.
Network storage (such as NFS and GPFS) are not recommended for this work, due to performance issues involving tens of thousands of small files.
Ensure that you have sufficient storage capacity.
For example, the GaiaSource DR2 files take 550 GB of space, and the converted LSST reference catalog takes another 200 GB.

To write the config for the converter, you will need to have a document that describes the columns in the input data, including their units and any caveats that may apply (for example, Gaia DR2 supplies magnitudes in the Vega system and no magnitude errors, so we convert the native fluxes and errors in e-/s to AB magnitudes with a customized class).
In this example, we used the `Gaia Source Catalog data model <https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html>`_ document.

If the files are text files of some sort, check that you can read one of them with `astropy.io.ascii.read`, which is what the converter uses to read text files. For example:

.. code-block:: python

    import astropy.io.ascii
    data = astropy.io.ascii.read('gaia_source/GaiaSource_1000172165251650944_1000424567594791808.csv.gz', format='csv')
    print(data)

The default Config assumes that the files are readable with ``format="csv"``; you can change that to a different ``format`` if necessary (see `~lsst.meas.algorithms.readTextCatalogTask.ReadTextCatalogConfig` for how to configure that file reader, and `~lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogConfig`).

.. _lsst.meas.algorithms-refcat-config:

2. Write a README for the reference catalog
===========================================
All reference catalogs should have a README file that describes the catalog, its provenance, and how it was converted.
This README should be stored in the same directory as the reference catalog files.
An example of a README can be found in the `Common Dataset Organization and Policy <https://developer.lsst.io/usdf/datasets.html#reference-catalogs>` page.

3. Write a Config for the conversion
====================================

`~lsst.meas.algorithms.ConvertReferenceCatalogConfig` specifies what fields in the input files get translated to the output data, and how they are converted along the way.
See the :ref:`ConvertReferenceCatalogTask configuration documentation <lsst.meas.algorithms.ConvertReferenceCatalogTask-configs>` docs for the different available options.

This is an example configuration that was used to convert the Gaia DR2 catalog (saved as ``gaia_dr2_config.py`` in the current directory and provided to the next step as a ``configFile``):

.. code-block:: python

    # The name of the output reference catalog dataset.
    config.dataset_config.ref_dataset_name = "gaia_dr2"

    # Gaia has a specialized convert manager class to handle the flux conversion.
    from lsst.meas.algorithms import convertRefcatManager
    config.manager.retarget(convertRefcatManager.ConvertGaiaManager)

    # Ingest the data in parallel with this many processes; this is sized to
    # fill a single node on lsst-devl.
    config.n_processes = 48

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

    # This is a required config field, and is used to populate the output schema:
    # we append `_flux` and `_fluxErr` to them in the output schema.
    # The Gaia-specific convert manager class handles the flux/flux error math,
    # using the flux fields (that are in e-/s units).
    config.mag_column_list = ["phot_g_mean", "phot_bp_mean", "phot_rp_mean"]

    # These fields are brought along unmodified.
    config.extra_col_names = ["astrometric_excess_noise", "phot_variable_flag"]

In order to deal with the way that Gaia released their photometric data, we have subclassed the conversion manager as `~lsst.meas.algorithms.convertRefcatManager.ConvertGaiaManager`.
This class special-cases the calculation of the flux and flux errors from the values in the Gaia DR2 catalog, which cannot be handled via the simple Config system used above.

.. _lsst.meas.algorithms-refcat-convert:

3. Convert the files to the LSST format
=======================================

:doc:`scripts/convertReferenceCatalog` takes three parameters: output path, configuration file, and quoted input file glob.
See the command line reference for more details about these parameters.

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

Check the log output after several hours: ``ConvertReferenceCatalogTask`` reports progress in 1% intervals.

.. _lsst.meas.algorithms-refcat-ingest:

4. Ingest the files into the butler
===================================

When ``convertReferenceCatalog`` has finished, a new directory (named ``gaia-refcat/`` in the example above) will now exist containing the HTM-indexed files for the input catalog in the LSST format.
Our convention for reference catalogs is described in `DMTN-167 <https://dmtn-167.lsst.io/#reference-catalogs>`_ and is used in the commands listed below.

Three final steps are now required to register the new refcat dataset type, ingest your converted output into a RUN collection, and CHAIN that collection to the primary ``refcats`` collection.
If using ``convertReferenceCatalog``, these commands will have been printed on the command line when it finished.
For the example we are using here, these commands would be:

.. prompt:: bash

    butler register-dataset-type REPO gaia_dr2_20200414 SimpleCatalog htm7
    butler ingest-files -t direct REPO gaia_dr2_20200414 refcats/DM-NNNNN gaia/filename_to_htm.ecsv
    butler collection-chain REPO --mode extend refcats refcats/DM-NNNNN

where ``REPO`` is the path to the butler repository that you are ingesting the data into, and ``DM-NNNNN`` is the ticket you are tracking this refcat ingest on.
The name of the reference catalog is used as the dataset type and should include a date string (``_2020041411`` in this example) to distinguish between different versions of the same reference catalog.
We use the ``direct`` transfer mode here to leave the files in the directory they were converted into: ``gaia_dr2/``.
See ``butler ingest-files -h`` for other options, including ``copy``, ``move`` and ``link`` transfer modes.

These commands should finish in a short amount of time, logging a message about how many files were ingested.
You can query the ``refcats`` collection to see whether your htm shards appear:

.. prompt:: bash

    butler query-datasets --collections refcats REPO

For LSST staff using ``lsst-devl``, see the `Reference catalogs policy <https://developer.lsst.io/services/datasets.html#reference-catalogs>`_ in the Developer Guide for additional policy about adding reference catalogs to the common repo.

.. _lsst.meas.algorithms-refcat-existing:

5. Ingesting pre-existing reference catalogs
============================================

.. note::

    The ``.ecsv`` files described here have already been created for the reference catalogs used in ``/repo/main`` at USDF, in the ``/sdf/group/rubin/datasets/refcats/htm/v1`` directory.
    If you wish to use these refcats in your own butler repo, just run the three commands at the end of this section; there is no need to generate the htm7 index lookup table file.

Already existing reference catalogs (for example, the PS1 or Gaia DR2 catalogs that were used in gen2 butlers) can be directly ingested into a gen3 repo as they are already in the LSST format.
To ingest already existing HTM-indexed files in the LSST format, first create a suitable filename to htm7-index astropy lookup table, and then follow the steps above to :ref:`ingest the files into the butler <lsst.meas.algorithms-refcat-ingest>`.

This is an example script that creates an ``.ecsv`` lookup table for the ``butler ingest-files`` command, from all of the HTM indexed files in a given directory (`refcat_dir` here). We use the existing Gaia DR2 catalog on lsst-devl in this example:

.. code-block:: python

    """Generate an astropy-readable .ecsv lookup file for `butler ingest-files`, to ingest an existing gen2 refcat.
    """

    import os
    import glob
    import astropy.table

    refcat_dir = "/datasets/refcats/htm/v1/gaia_dr2_20200414"
    out_dir = "/path/to/my/output/directory"

    out_file = f"{out_dir}/{os.path.basename(refcat_dir)}.ecsv"

    table = astropy.table.Table(names=("filename", "htm7"), dtype=("str", "int"))
    files = glob.glob(f"{refcat_dir}/[0-9]*.fits")

    for i, file in enumerate(files):
        # running status, overwriting each print statement as it proceeds
        print(f"{i}/{len(files)} ({100*i/len(files):0.1f}%)", end="\r")

        # extract file index; add row to table
        file_index = int(os.path.basename(os.path.splitext(file)[0]))
        table.add_row((file, file_index))

    table.write(out_file)
    print(f"Output written to: {out_file}")

Once this script is complete, finish reference catalog ingestion by following the :ref:`file ingestion instructions above <lsst.meas.algorithms-refcat-ingest>`.
In particular, you need to change the name of the registered dataset type to "gaia_dr2_20200414" (the reference catalog used in the Python code block above), and the filename to the generated .ecsv file ("gaia_dr2_20200414.ecsv" in the Python code block above):

.. prompt:: bash

    butler register-dataset-type REPO gaia_dr2_20200414 SimpleCatalog htm7
    butler ingest-files -t direct REPO gaia_dr2_20200414 refcats/DM-NNNNN gaia_dr2_20200414.ecsv
    butler collection-chain REPO --mode extend refcats refcats/DM-NNNNN
