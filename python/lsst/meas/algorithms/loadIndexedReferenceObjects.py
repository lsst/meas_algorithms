#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ["LoadIndexedReferenceObjectsConfig", "LoadIndexedReferenceObjectsTask"]

from lsst.meas.algorithms import LoadReferenceObjectsTask, LoadReferenceObjectsConfig
import lsst.pex.config as pexConfig
from deprecated.sphinx import deprecated


@deprecated(reason=("This config is no longer used; it will be removed after v25. "
                    "Please use LoadReferenceObjectsConfig instead."),
            version="v25.0", category=FutureWarning)
class LoadIndexedReferenceObjectsConfig(LoadReferenceObjectsConfig):
    ref_dataset_name = pexConfig.Field(
        dtype=str,
        default='cal_ref_cat',
        doc='Name of the ingested reference dataset',
        deprecated='This field is no longer used. It will be removed after v25.',
    )


@deprecated(reason=("This task is used in gen2 only; it will be removed after v25. "
                    "See DM-35671 for details on updating code to avoid this warning."),
            version="v25.0", category=FutureWarning)
class LoadIndexedReferenceObjectsTask(LoadReferenceObjectsTask):
    """Stub of the LoadIndexedReferenceObjectsTask to allow retargeting before removal."""
    ConfigClass = LoadIndexedReferenceObjectsConfig
    _DefaultName = 'LoadIndexedReferenceObjectsTask'
