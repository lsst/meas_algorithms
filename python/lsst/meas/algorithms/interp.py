from . import interpolateOverDefectsOld
from . import interpolateOverDefectsGP

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
        return interpolateOverDefectsOld(
            image, psf, badList, fallbackValue, useFallbackValueAtEdge
        )
    else:
        return interpolateOverDefectsGP(image, fwhm, maskNameList, **kwargs)
