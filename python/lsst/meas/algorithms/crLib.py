__all__ = ["findCosmicRays"]

import warnings

import lsst.pex.policy as pexPolicy
from .cr import _findCosmicRays


def findCosmicRays(image, psf, bkgd, ps, keep=False):
    """Find cosmic rays in an image."""
    if isinstance(ps, pexPolicy.Policy):
        warnings.warn("findCosmicRays acceptance of Policy is deprecated. Use PropertySet "
                      "(support will be removed after v19.0)",
                      category=FutureWarning, stacklevel=2)
        ps = ps.asPropertySet()
    return _findCosmicRays(image, psf, bkgd, ps, keep)


findCosmicRays.__doc__ = _findCosmicRays.__doc__
