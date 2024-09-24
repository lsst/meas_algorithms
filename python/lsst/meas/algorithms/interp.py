from . import legacyInterpolateOverDefects
from . import InterpolateOverDefectGaussianProcess

__all__ = ["interpolateOverDefects"]


def interpolateOverDefects(
    image,
    psf,
    badList,
    fallbackValue=0.0,
    useFallbackValueAtEdge=False,
    fwhm=1.0,
    useLegacyInterp=True,
    maskNameList=None,
    **kwargs
):
    if useLegacyInterp:
        return legacyInterpolateOverDefects(
            image, psf, badList, fallbackValue, useFallbackValueAtEdge
        )
    else:
        gp = InterpolateOverDefectGaussianProcess(image, fwhm=fwhm,
                                                  defects=maskNameList, **kwargs)
        return gp.run()
