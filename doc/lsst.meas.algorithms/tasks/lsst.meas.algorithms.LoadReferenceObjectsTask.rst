.. lsst-task-topic:: lsst.meas.algorithms.LoadReferenceObjectsTask

########################
LoadReferenceObjectsTask
########################

Implementations must subclass this class, override the loadSkyCircle method, and will typically override the value of ConfigClass with a task-specific config class.

``LoadReferenceObjectsTask`` acts as an abstract base class for tasks that load objects from a reference catalog in a particluar region of sky.

.. _lsst.meas.algorithms.LoadReferenceObjectsTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.LoadReferenceObjectsTask

.. _lsst.meas.algorithms.LoadReferenceObjectsTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.algorithms.LoadReferenceObjectsTask

.. _lsst.meas.algorithms.LoadReferenceObjectsTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.algorithms.LoadReferenceObjectsTask

.. _lsst.meas.algorithms.LoadReferenceObjectsTask-indepth:

In Depth
========

Catalog Schema
--------------

    - ``coord``: ICRS position of star on sky (an ``lsst.geom.SpherePoint``)
    - ``centroid``: position of star on an exposure, if relevant (an ``lsst.afw.Point2D``)
    - ``hasCentroid``: is centroid usable? (a Flag)
    - ``<referenceFilterName>_flux``: brightness in the specified reference catalog filter (nJy)
        Note: you can use astropy.units to convert from AB Magnitude to nJy: `u.Magnitude(value, u.ABmag).to_value(u.nJy)`
    - ``<referenceFilterName>_fluxErr`` (optional): brightness standard deviation (nJy); omitted if no data is available; possibly nan if data is available for some objects but not others
    - ``camFlux``: brightness in default camera filter (nJy); omitted if `defaultFilter` not specified
    - ``camFluxErr``: brightness standard deviation for default camera filter; omitted if `defaultFilter` not specified or standard deviation not available that filter
    - ``<cameraFilterName>_camFlux``: brightness in specified camera filter (nJy)
    - ``<cameraFilterName>_camFluxErr`` (optional): brightness standard deviation in specified camera filter (nJy); omitted if no data is available;
        possibly nan if data is available for some objects but not others
    - ``photometric`` (optional): is the object usable for photometric calibration? (a Flag)
    - ``resolved`` (optional): is the object spatially resolved? (a Flag)
    - ``variable`` (optional): does the object have variable brightness? (a Flag)
    - ``coord_raErr``: uncertainty in `coord` along the direction of right ascension (radian, an Angle) = uncertainty in ra * cos(dec); nan if unknown
    - ``coord_decErr``: uncertainty in `coord` along the direction of declination (radian, an Angle);
        nan if unknown

The following are optional; fields should only be present if the information is available for at least some objects.
    Numeric values are `nan` if unknown:

    - ``epoch``: date of observation as TAI MJD (day)
    - ``pm_ra``: proper motion along the direction of right ascension (rad/year, an Angle) = dra/dt * cos(dec)
    - ``pm_dec``: proper motion along the direction of declination (rad/year, and Angle)
    - ``pm_raErr``: uncertainty in `pm_ra` (rad/year)
    - ``pm_decErr``: uncertainty in `pm_dec` (rad/year)
    - ``pm_ra_dec_Cov``: covariance between pm_ra and pm_dec (rad2/year2)
    - ``pm_flag``: set if proper motion, error or covariance is bad
    - ``parallax``: parallax (rad, an Angle)
    - ``parallaxErr``: uncertainty in `parallax` (rad)
    - ``parallax_flag``: set if parallax value or parallaxErr is bad
    - ``coord_ra_pm_ra_Cov``: covariance between coord_ra and pm_ra (rad2/year)
    - ``coord_ra_pm_dec_Cov``: covariance between coord_ra and pm_dec (rad2/year)
    - ``coord_ra_parallax_Cov``: covariance between coord_ra and parallax (rad2/year)
    - ``coord_dec_pm_ra_Cov``: covariance between coord_dec and pm_ra (rad2/year)
    - ``coord_dec_pm_dec_Cov``: covariance between coord_dec and pm_dec (rad2/year)
    - ``coord_dec_parallax_Cov``: covariance between coord_dec and parallax (rad2/year)
    - ``pm_ra_parallax_Cov``: covariance between pm_ra and parallax (rad2/year)
    - ``pm_dec_parallax_Cov``: covariance between pm_dec and parallax (rad2/year)
