.. lsst-task-topic:: lsst.meas.algorithms.installGaussianPsf.InstallGaussianPsfTask

######################
InstallGaussianPsfTask
######################

``InstallGaussianPsfTask`` installs a Gaussian PSF model in an exposure.

.. _lsst.meas.algorithms.installGaussianPsf.InstallGaussianPsfTask-summary:

Processing summary
==================

``InstallGaussianPsfTask`` installs a Gaussian PSF model in an exposure, creating a new PSF with the same sigma and width if a PSF already exists for the exposure. If there is not a PSF model for the exposure, the sigma and width are taken from the config.

.. _lsst.meas.algorithms.installGaussianPsf.InstallGaussianPsfTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.algorithms.installGaussianPsf.InstallGaussianPsfTask

.. _lsst.meas.algorithms.installGaussianPsf.InstallGaussianPsfTask-examples:

Examples
========

.. code-block:: python

    from lsst.afw.image import ExposureF
    from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask, FwhmPerSigma

    exposure = ExposureF(100, 100)
    task = InstallGaussianPsfTask()
    task.run(exposure=exposure)

     # This particular exposure had no PSF model to begin with, so the new PSF model
     # uses the config's FWHM. However, measured FWHM is based on the truncated
     # PSF image, so it does not exactly match the input
     measFwhm = exposure.getPsf().computeShape().getDeterminantRadius() * FwhmPerSigma
     assert abs(measFwhm - task.config.fwhm) < 1e-3
