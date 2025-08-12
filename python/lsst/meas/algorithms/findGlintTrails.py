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

__all__ = ["FindGlintTrailsConfig", "FindGlintTrailsTask", "GlintTrailParameters"]

import collections
import dataclasses
import math

import numpy as np
import scipy.spatial
import sklearn.linear_model

import lsst.afw.table
import lsst.pex.config
import lsst.pipe.base


class FindGlintTrailsConfig(lsst.pex.config.Config):
    radius = lsst.pex.config.Field(
        doc="Radius to search for glint trail candidates from each source (pixels).",
        dtype=float,
        default=500,
    )
    min_points = lsst.pex.config.Field(
        doc="Minimum number of points to be considered a possible glint trail.",
        dtype=int,
        default=5,
        check=lambda x: x >= 3,
    )
    threshold = lsst.pex.config.Field(
        doc="Maximum root mean squared deviation from a straight line (pixels).",
        dtype=float,
        default=15.0,
    )
    seed = lsst.pex.config.Field(
        doc="Random seed for RANSAC fitter, to ensure stable fitting.",
        dtype=int,
        default=42,
    )
    bad_flags = lsst.pex.config.ListField[str](
        doc="Do not fit sources that have these flags set.",
        default=["ip_diffim_DipoleFit_classification",
                 "is_negative",
                 ],
    )


@dataclasses.dataclass(frozen=True, kw_only=True)
class GlintTrailParameters:
    """Holds values from the line fit to a single glint trail."""
    slope: float
    intercept: float
    stderr: float
    length: float  # pixels
    angle: float  # radians, from +X axis


class FindGlintTrailsTask(lsst.pipe.base.Task):
    """Find glint trails in a catalog by searching for sources that lie in a
    line.

    Notes
    -----
    For each source ("anchor") in the input catalog that was not included in
    an an earlier iteration as part of a trail:
     * Find all sources within a given radius.
     * For each pair of anchor and match, identify the other sources that
       could lie on the same line(s).
     * Take the longest set of such pairs as a candidate trail.
     * Fit a line to the identified pairs with the RANSAC algorithm.
     * Find all sources in the catalog that could lie on that line.
     * Refit a line to all of the matched sources.
     * If the error is below the threshold and the number of sources on the
       line is greater than the minimum, return the sources that were
       considered inliers during the fit, and the fit parameters.
    """

    ConfigClass = FindGlintTrailsConfig
    _DefaultName = "findGlintTrails"

    def run(self, catalog):
        """Find glint trails in a catalog.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            Catalog to search for glint trails.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results as a struct with attributes:

            ``trails``
                Catalog subsets containing sources in each trail that was found.
                (`list` [`lsst.afw.table.SourceCatalog`])
            ``trailed_ids``
                Ids of all the sources that were included in any fit trail.
                (`set` [`int`])
            ``parameters``
                Parameters of all the trails that were found.
                (`list` [`GlintTrailParameters`])
        """
        good_catalog = self._select_good_sources(catalog)

        matches = lsst.afw.table.matchXy(good_catalog, self.config.radius)
        per_id = collections.defaultdict(list)
        for match in matches:
            per_id[match.first["id"]].append(match)
        counts = {id: len(value) for id, value in per_id.items()}

        trails = []
        parameters = []
        trailed_ids = set()
        # Search starting with the source with the largest number of matches.
        for id in dict(sorted(counts.items(), key=lambda item: item[1], reverse=True)):
            # Don't search this point if it was already included in a trail.
            if counts[id] < self.config.min_points or id in trailed_ids:
                continue

            self.log.debug("id=%d at %.1f,%.1f has %d matches within %d pixels.",
                           id,
                           per_id[id][0].first.getX(),
                           per_id[id][0].first.getY(),
                           counts[id],
                           self.config.radius)
            if (trail := self._search_one(per_id[id], good_catalog)) is not None:
                trail, result = trail
                # Check that we didn't already find this trail.
                n_new = len(set(trail["id"]).difference(trailed_ids))
                if n_new > 0:
                    self.log.info("Found %.1f pixel length trail with %d points, "
                                  "%d not in any other trail (slope=%.4f, intercept=%.2f)",
                                  result.length, len(trail), n_new, result.slope, result.intercept)
                    trails.append(trail)
                    trailed_ids.update(trail["id"])
                    parameters.append(result)

        self.log.info("Found %d glint trails containing %d total sources.",
                      len(trails), len(trailed_ids))
        return lsst.pipe.base.Struct(trails=trails,
                                     trailed_ids=trailed_ids,
                                     parameters=parameters)

    def _select_good_sources(self, catalog):
        """Return sources that could possibly be in a glint trail, i.e. ones
        that do not have bad flags set.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            Original catalog to be selected from.

        Returns
        -------
        good_catalog : `lsst.afw.table.SourceCatalog`
            Catalog that has had bad sources removed.
        """
        bad = np.zeros(len(catalog), dtype=bool)
        for flag in self.config.bad_flags:
            bad |= catalog[flag]
        return catalog[~bad]

    def _search_one(self, matches, catalog):
        """Search one set of matches for a possible trail.

        Parameters
        ----------
        matches : `list` [`lsst.afw.table.Match`]
            Matches for one anchor source to search for lines.
        catalog : `lsst.afw.SourceCatalog`
            Catalog of all sources, to refit lines to.

        Returns
        -------
        trail, result : `tuple` or None
            If the no trails matching the criteria are found, return None,
            otherwise return a tuple of the sources in the trail and the
            trail parameters.
        """
        components = collections.defaultdict(list)
        # Normalized distances from the first record to all the others.
        xy_deltas = {pair.second["id"]: (pair.second.getX() - pair.first.getX(),
                                         pair.second.getY() - pair.first.getY()) for pair in matches}

        # Find all sets of pairs from this anchor that could lie on a line.
        for i, (id1, pair1) in enumerate(xy_deltas.items()):
            distance = math.sqrt(pair1[0]**2 + pair1[1]**2)
            for j, (id2, pair2) in enumerate(xy_deltas.items()):
                if i == j:
                    continue
                delta = abs(pair1[0] * pair2[1] - pair1[1] * pair2[0])
                # 2x threshold to search more broadly; will be refined later.
                if delta / distance < 2 * self.config.threshold:
                    components[i].append(j)

        # There are no lines with at least 3 components.
        if len(components) == 0:
            return None

        longest, value = max(components.items(), key=lambda x: len(x[1]))
        n_points = len(value)
        n_points += 2  # to account for the base source and the first pair
        if n_points < self.config.min_points:
            return None

        candidate = [longest] + components[longest]
        trail, result = self._other_points(n_points, candidate, matches, catalog)

        if trail is None or len(trail) < self.config.min_points:
            return None
        if result.stderr > self.config.threshold:
            self.log.info("Candidate trail with %d sources rejected with stderr %.6f > %.3f",
                          len(trail), result.stderr, self.config.threshold)
            return None
        else:
            return trail, result

    def _other_points(self, n_points, indexes, matches, catalog):
        """Find all catalog records that could lie on this line.

        Parameters
        ----------
        n_points : `int`
            Number of sources in this candidate trail.
        indexes : `list` [`int`]
            Indexes into matches on this candidate trail.
        matches : `list` [`lsst.afw.table.Match`]
            Matches for one anchor sources to search for lines.
        catalog : `lsst.afw.SourceCatalog`
            Catalog of all sources, to refit lines to.

        Returns
        -------
        trail : `lsst.afw.table.SourceCatalog`
            Sources that are in the fitted trail.
        result : `GlintTrailParameters`
            Parameters of the fitted trail.
        """

        def extract(fitter, x, y, prefix=""):
            """Extract values from the fit and log and return them."""
            x = x[fitter.inlier_mask_]
            y = y[fitter.inlier_mask_]
            predicted = fitter.predict(x).flatten()
            stderr = math.sqrt(((predicted - y.flatten())**2).sum())
            m, b = fitter.estimator_.coef_[0][0], fitter.estimator_.intercept_[0]
            self.log.debug("%s fit: score=%.6f, stderr=%.6f, inliers/total=%d/%d",
                           prefix, fitter.score(x, y), stderr, sum(fitter.inlier_mask_), len(x))
            # Simple O(N^2) search for longest distance; there will never be
            # enough points in a trail a for "faster" approach to be worth it.
            length = max(scipy.spatial.distance.pdist(np.hstack((x, y))), default=0)
            angle = math.atan(m)
            return GlintTrailParameters(slope=m, intercept=b, stderr=stderr, length=length, angle=angle)

        # min_samples=2 is necessary here for some sets of only 5 matches,
        # otherwise we sometimes get "UndefinedMetricWarning: R^2 score is not
        # well-defined with less than two samples" from RANSAC.
        fitter = sklearn.linear_model.RANSACRegressor(residual_threshold=self.config.threshold,
                                                      loss="squared_error",
                                                      random_state=self.config.seed,
                                                      min_samples=2)

        # The (-1,1) shape is to keep sklearn happy.
        x = np.empty(n_points).reshape(-1, 1)
        x[0] = matches[0].first.getX()
        x[1:, 0] = [matches[i].second.getX() for i in indexes]
        y = np.empty(n_points).reshape(-1, 1)
        y[0] = matches[0].first.getY()
        y[1:, 0] = [matches[i].second.getY() for i in indexes]

        fitter.fit(x, y)
        result = extract(fitter, x, y, prefix="preliminary")
        # Reject trails that have too many outliers after the first fit.
        if (n_inliers := sum(fitter.inlier_mask_)) < self.config.min_points:
            self.log.debug("Candidate trail rejected with %d < %d points.", n_inliers, self.config.min_points)
            return None, None

        # Find all points that are close to this line and refit with them.
        x = catalog["slot_Centroid_x"]
        y = catalog["slot_Centroid_y"]
        dist = abs(result.intercept + result.slope * x - y) / math.sqrt(1 + result.slope**2)
        # 2x threshold to search more broadly: outlier rejection may change
        # the line parameters some and we want to grab all candidates here.
        candidates = (dist < 2 * self.config.threshold).flatten()
        # min_samples>2 should make the fit more stable.
        fitter = sklearn.linear_model.RANSACRegressor(residual_threshold=self.config.threshold,
                                                      loss="squared_error",
                                                      random_state=self.config.seed,
                                                      min_samples=3)
        # The (-1,1) shape is to keep sklearn happy.
        x = x[candidates].reshape(-1, 1)
        y = y[candidates].reshape(-1, 1)
        fitter.fit(x, y)
        result = extract(fitter, x, y, prefix="final")

        return catalog[candidates][fitter.inlier_mask_], result
