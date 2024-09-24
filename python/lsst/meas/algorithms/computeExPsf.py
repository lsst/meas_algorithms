# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


from lsst.meas.algorithms.treecorrUtils import TreecorrConfig
from lsst.pipe.base import Task
import lsst.pipe.base as pipeBase
import treecorr
import copy
import numpy.typing as npt


__all__ = ("ComputeExPsfTask")

class ComputeExPsfTask(Task):
    """Compute Ex for PSF.

    Compute scalar correlation function from
    PSF ellipticity residuals to compute TEx
    metrics.

    Parameters
    ----------
    de1: `np.ndarray`
        PSF ellipticity residuals component 1.
    de2: `np.ndarray`
        PSF ellipticity residuals component 2.
    ra: `np.ndarray`
        Right ascension coordinate.
    dec: `np.ndarray`
        Declination coordinate.
    units: `str`
        In which units are ra and dec. units supported
        are the same as the one in treecorr.

    Returns
    -------
    struct : `lsst.pipe.base.Struct`
            The struct contains the following data:
        ``E1``: `float`
            <de1 de1> scalar correlation function, compute
            in an angular bin define in TreecorrConfig.
        ``E2``: `float`
            <de2 de2> scalar correlation function, compute
            in an angular bin define in TreecorrConfig.
        ``Ex``: `float`
            <de1 de2> scalar cross-correlation function, compute
            in an angular bin define in TreecorrConfig.
    """

    ConfigClass = TreecorrConfig
    _DefaultName = "computeExPsf"

    def run(
        self,
        de1: npt.NDArray,
        de2: npt.NDArray,
        ra: npt.NDArray,
        dec: npt.NDArray,
        units: str = "arcmin",
    ) -> pipeBase.Struct:

        kwargs_cat = {
            "ra": ra,
            "dec": dec,
            "ra_units": units,
            "dec_units": units,
        }

        cat1 = treecorr.Catalog(k=de1, **kwargs_cat)
        cat2 = treecorr.Catalog(k=de2, **kwargs_cat)

        config_kk = self.config.toDict()

        kk = treecorr.KKCorrelation(config_kk)

        kk.process(cat1)
        kk_E1 = copy.deepcopy(kk.xi[0])
        kk.process(cat2)
        kk_E2 = copy.deepcopy(kk.xi[0])
        kk.process(cat1, cat2)
        kk_Ex = copy.deepcopy(kk.xi[0])

        return pipeBase.Struct(E1=kk_E1, E2=kk_E2, Ex=kk_Ex)
