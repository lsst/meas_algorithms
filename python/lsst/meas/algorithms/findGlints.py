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

__all__ = ["FindGlintsConfig", "FindGlintsTask"]

import collections
import math
import scipy.stats

import lsst.afw.table
import lsst.pex.config
import lsst.pipe.base


class FindGlintsConfig(lsst.pex.config.Config):
    min_points = lsst.pex.config.Field(
        doc="Minimum number of points in a given radius.",
        dtype=int,
        default=5,
        check=lambda x: x >= 3,
    )
    radius = lsst.pex.config.Field(
        doc="Radius to search for glint trails from each source (pixels).",
        dtype=float,
        default=500,
    )
    maxgcr = lsst.pex.config.Field(
        doc="Maximum RMS deviation from a straight line (pixels).",
        dtype=float,
        default=1.0,
    )


class FindGlintsTask(lsst.pipe.base.Task):
    ConfigClass = FindGlintsConfig
    _DefaultName = "findGlints"

    def run(self, catalog):
        matches = lsst.afw.table.matchXy(catalog, self.config.radius)
        per_id = collections.defaultdict(list)
        for match in matches:
            per_id[match.first["id"]].append(match)
        counts = {id: len(value) for id, value in per_id.items()}

        trails = []
        trailed_ids = set()
        for id in counts:
            # Don't search this point if it was already included in a trail.
            if counts[id] < self.config.min_points or id in trailed_ids:
                continue

            self.log.info("id=%d has %d matches", id, counts[id])
            if (trail := self._search_one(per_id[id], catalog)) is not None:
                self.log.info("Found glint trail with %d points.", len(trail))
                trails.append(trail)
                trailed_ids.update(trail["id"])

        # self.deduplicate(trails)
        return trails

    def _search_one(self, matches, catalog):
        """Search one set of matches for a possible trail."""
        components = collections.defaultdict(list)
        # Use normalized distances from the first record to all the others.
        xy_deltas = {pair.second["id"]: (pair.second.getX() - pair.first.getX(),
                                         pair.second.getY() - pair.first.getY()) for pair in matches}

        for i, (id1, pair1) in enumerate(xy_deltas.items()):
            distance = math.sqrt(pair1[0]**2 + pair1[1]**2)
            for j, (id2, pair2) in enumerate(xy_deltas.items()):
                if i == j:
                    continue
                delta = abs(pair1[0] * pair2[1] - pair1[1] * pair2[0])
                if delta / distance < 2 * self.config.maxgcr:
                    components[i].append(j)

        longest = -1
        length = 0
        # TODO: I think there's a itertools thing to do this?
        for id, matched in components.items():
            if len(matched) > length:
                length = len(matched)
                longest = id
        length += 2  # to account for the base source and the first pair
        if length < self.config.min_points:
            return None

        x = [matches[i].second.getX() for i in components[longest]]
        x.append(matches[0].first.getX())
        y = [matches[i].second.getY() for i in components[longest]]
        y.append(matches[0].first.getY())
        result = scipy.stats.linregress(x, y)
        print(result)

        points = self._other_points(result.slope, result.intercept, catalog)
        trail = catalog[points]

        # TODO: outlier rejection here!
        if result.stderr <= self.config.maxgcr and length >= self.config.min_points:
            # TODO: recenter on the middle?
            return trail
        else:
            return None

    def _other_points(self, m, b, catalog):
        """Find all catalog records that could lie on this line."""
        dist = abs(b + m * catalog.getX() - catalog.getY()) / math.sqrt(1 + m**2)
        # 2x max error to search more broadly: outlier rejection may change
        # the line parameters some.
        return dist < 2 * self.config.maxgcr

    def deduplicate(self, glints):
        """Remove duplicate glint trails."""
        pass
