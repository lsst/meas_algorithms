from . import interpolateOverDefectsOld

__all__ = ['interpolateOverDefects']


def interpolateOverDefects(image, psf, badList, fallbackValue=0.0,
                           useFallbackValueAtEdge=False, useNewAlgorithm=False, **kwargs):
    if useNewAlgorithm:
        raise RuntimeError("didn't write this.")
    else:
        return interpolateOverDefectsOld(image, psf, badList, fallbackValue, useFallbackValueAtEdge)
