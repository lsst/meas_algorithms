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

__all__ = ["FindGlintTrailsConfig", "FindGlintTrailsTask", "FittedLine", "GlintTrailParameters", "fit_line"]

import collections
import dataclasses
import math

import numpy as np
import scipy.spatial

import lsst.afw.table
import lsst.pex.config
import lsst.pipe.base


class FittedLine:
    """A straight line through ``center`` along the unit vector
    ``direction``.

    Parameters
    ----------
    center : `numpy.ndarray`
        ``(x, y)`` of a point on the line: the mean of the points it was fit
        to.
    direction : `numpy.ndarray`
        Unit vector along the line; its sign is arbitrary.
    """

    def __init__(self, center, direction):
        self.center = center
        self.direction = direction

    def distance(self, x, y):
        """Compute the perpendicular distance of points from the line.

        Parameters
        ----------
        x, y : `numpy.ndarray` or `float`
            Coordinates of the points.

        Returns
        -------
        distance : `numpy.ndarray` or `float`
            Distance of each point from the line, in the units of ``x`` and
            ``y``.
        """
        return np.abs((x - self.center[0])*self.direction[1] - (y - self.center[1])*self.direction[0])


def fit_line(x, y):
    """Fit the line that minimizes the summed squared perpendicular distances
    of points from it.

    Unlike a regression of ``y`` on ``x``, this fits lines at any angle,
    including vertical ones, equally well.

    Parameters
    ----------
    x, y : `numpy.ndarray`
        Coordinates of at least two points.

    Returns
    -------
    line : `FittedLine`
        The fitted line.
    """
    points = np.column_stack((x, y))
    center = points.mean(axis=0)
    direction = np.linalg.svd(points - center)[2][0]
    return FittedLine(center=center, direction=direction)


def _ransac_line(x, y, threshold, min_samples, seed, max_trials=100):
    """Fit a line to points that may include outliers, using RANSAC.

    Each trial fits a line to ``min_samples`` randomly chosen points.  The
    trial with the most points within ``threshold`` of its line wins, with
    ties going to the smaller summed squared distance.  The line is iteratively
    refit with new selections of inliers until the results stop changing, so
    that every inlier is within ``threshold`` pixels of the returned line.

    Parameters
    ----------
    x, y : `numpy.ndarray`
        Coordinates of the points.
    threshold : `float`
        Maximum perpendicular distance of an inlier from the line.
    min_samples : `int`
        Number of points each trial line is fit to.
    seed : `int`
        Seed for the random choice of points, so that fits are repeatable.
    max_trials : `int`, optional
        Maximum number of trial lines.

    Returns
    -------
    line : `FittedLine` or `None`
        The line fit to the inliers, or `None` if fewer than ``min_samples``
        points are inliers.
    inliers : `numpy.ndarray` [`bool`] or `None`
        Which of the points lie within ``threshold`` of ``line``.
    """
    if len(x) < min_samples:
        return None, None
    rng = np.random.default_rng(seed)
    best = None
    best_count = min_samples - 1
    best_sum = np.inf
    for _ in range(max_trials):
        sample = rng.choice(len(x), size=min_samples, replace=False)
        distance = fit_line(x[sample], y[sample]).distance(x, y)
        inliers = distance <= threshold
        count = np.sum(inliers)
        total = np.sum(distance[inliers]**2)
        if count > best_count or (count == best_count and total < best_sum):
            best, best_count, best_sum = inliers, count, total
            if count == len(x):
                break
    if best is None:
        return None, None
    inliers = best
    for _ in range(max_trials):
        line = fit_line(x[inliers], y[inliers])
        refit = line.distance(x, y) <= threshold
        if np.sum(refit) < min_samples:
            return None, None
        if np.array_equal(refit, inliers):
            break
        inliers = refit
    return line, refit


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
        doc="Maximum perpendicular distance of a source from the line fit to its trail (pixels)."
            " Real glints are highly collinear, with a typical deviation from a line of < 1 pixel.",
        dtype=float,
        default=1.5,
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
    """Holds values from the line fit to a single glint trail.

    ``slope`` and ``intercept`` are very large or non-finite for a
    near-vertical trail; ``angle`` with ``center_x`` and ``center_y`` describe
    the line at any orientation.
    """
    slope: float
    intercept: float
    stderr: float  # root mean squared perpendicular distance, pixels
    length: float  # pixels
    angle: float  # radians from the +X axis, in (-pi/2, pi/2]
    center_x: float  # pixels; the line passes through (center_x, center_y)
    center_y: float  # pixels


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
     * Fit a line to the identified pairs with the RANSAC algorithm,
       measuring distances perpendicular to the line so that trails at any
       angle are found.
     * Find all sources in the catalog that could lie on that line.
     * Refit a line to all of the matched sources.
     * If the line holds at least the minimum number of sources, return the
       sources that were considered inliers during the fit, and the fit
       parameters.
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
                                  "%d not in any other trail (angle=%.2f deg, center=%.1f,%.1f)",
                                  result.length, len(trail), n_new, math.degrees(result.angle),
                                  result.center_x, result.center_y)
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

        def extract(line, inliers, x, y, prefix=""):
            """Extract values from the fit and log and return them."""
            x = x[inliers]
            y = y[inliers]
            stderr = math.sqrt(np.mean(line.distance(x, y)**2))
            self.log.debug("%s fit: stderr=%.6f, inliers/total=%d/%d",
                           prefix, stderr, len(x), len(inliers))
            # Simple O(N^2) search for longest distance; there will never be
            # enough points in a trail a for "faster" approach to be worth it.
            length = max(scipy.spatial.distance.pdist(np.column_stack((x, y))), default=0)
            dx, dy = line.direction
            # The slope and angle do not depend on the direction's sign.
            if dx == 0:
                m, b = math.inf, math.nan
            else:
                m = dy/dx
                b = line.center[1] - m*line.center[0]
            angle = math.atan(m)
            return GlintTrailParameters(slope=m, intercept=b, stderr=stderr, length=length, angle=angle,
                                        center_x=line.center[0], center_y=line.center[1])

        x = np.empty(n_points)
        x[0] = matches[0].first.getX()
        x[1:] = [matches[i].second.getX() for i in indexes]
        y = np.empty(n_points)
        y[0] = matches[0].first.getY()
        y[1:] = [matches[i].second.getY() for i in indexes]

        line, inliers = _ransac_line(x, y, self.config.threshold, min_samples=2, seed=self.config.seed)
        if line is None:
            self.log.info("Glint trail interpolation could not find a valid fit.")
            return None, None
        extract(line, inliers, x, y, prefix="preliminary")
        # Reject trails that have too many outliers after the first fit.
        if (n_inliers := np.sum(inliers)) < self.config.min_points:
            self.log.debug("Candidate trail rejected with %d < %d points.",
                           n_inliers, self.config.min_points)
            return None, None

        # Find all points that are close to this line and refit with them.
        x = catalog["slot_Centroid_x"]
        y = catalog["slot_Centroid_y"]
        # 2x threshold to search more broadly: outlier rejection may change
        # the line parameters some and we want to grab all candidates here.
        candidates = line.distance(x, y) < 2*self.config.threshold
        x = x[candidates]
        y = y[candidates]
        # min_samples>2 should make the fit more stable.
        line, inliers = _ransac_line(x, y, self.config.threshold, min_samples=3, seed=self.config.seed)
        if line is None:
            self.log.info("Glint trail interpolation could not find a valid fit.")
            return None, None
        result = extract(line, inliers, x, y, prefix="final")

        return catalog[candidates][inliers], result
