.. _reference-catalog-schema:

#########################################
LSST reference catalog schema description
#########################################

The reference catalogs for use with LSST Data Management Science Pipelines have the following schema:

Catalog Schema
--------------

    - ``coord_ra``: ICRS position of star on sky (an ``lsst.geom.Angle``)
    - ``coord_dec``: ICRS position of star on sky (an ``lsst.geom.Angle``)
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

    - ``coord_ra_coord_dec_Cov``: covariance between ``coord_ra`` and ``coord_dec`` (rad^2)
    - ``epoch``: date of observation as TAI MJD (day)
    - ``pm_ra``: proper motion along the direction of right ascension (rad/year, an Angle) = dra/dt * cos(dec)
    - ``pm_dec``: proper motion along the direction of declination (rad/year, and Angle)
    - ``pm_raErr``: uncertainty in `pm_ra` (rad/year)
    - ``pm_decErr``: uncertainty in `pm_dec` (rad/year)
    - ``pm_ra_pm_dec_Cov``: covariance between pm_ra and pm_dec (rad^2/year^2)
    - ``pm_flag``: set if proper motion, error or covariance is bad
    - ``parallax``: parallax (rad, an Angle)
    - ``parallaxErr``: uncertainty in `parallax` (rad)
    - ``parallax_flag``: set if parallax value or parallaxErr is bad
    - ``coord_ra_pm_ra_Cov``: covariance between coord_ra and pm_ra (rad^2/year)
    - ``coord_ra_pm_dec_Cov``: covariance between coord_ra and pm_dec (rad^2/year)
    - ``coord_ra_parallax_Cov``: covariance between coord_ra and parallax (rad^2)
    - ``coord_dec_pm_ra_Cov``: covariance between coord_dec and pm_ra (rad^2/year)
    - ``coord_dec_pm_dec_Cov``: covariance between coord_dec and pm_dec (rad^2/year)
    - ``coord_dec_parallax_Cov``: covariance between coord_dec and parallax (rad^2)
    - ``pm_ra_parallax_Cov``: covariance between pm_ra and parallax (rad^2/year)
    - ``pm_dec_parallax_Cov``: covariance between pm_dec and parallax (rad^2/year)
