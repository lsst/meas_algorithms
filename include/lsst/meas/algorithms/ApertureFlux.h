// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#if !defined(LSST_MEAS_ALGORITHMS_DETAIL_APERTUREPHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_DETAIL_APERTUREPHOTOMETRY_H 1
/**
 * @brief Compute a set of apertures using a naive algorithm
 * @ingroup meas/algorithms
 */
namespace lsst {
namespace meas {
namespace algorithms {

/**
 * Implement "Aperture" photometry.
 * @brief A class that knows how to calculate fluxes as a simple sum over a Footprint
 */
class ApertureFlux : public Algorithm {
public:
    typedef std::vector<double> VectorD;

    ApertureFlux(ApertureFluxControl const & ctrl, afw::table::Schema & schema) :
        Algorithm(ctrl),
        _fluxKey(schema.addField< afw::table::Array<double> >(
                ctrl.name, "sum of pixels in apertures", "dn", ctrl.radii.size())),
        _errKey(schema.addField< afw::table::Array<double> >(
                ctrl.name + ".err", "uncertainty for " + ctrl.name, "dn", ctrl.radii.size())),
        _nProfileKey(schema.addField<int>(ctrl.name + ".nProfile", "pixels",
                                          "Number of points in radial profile")),
        _flagKey(schema.addField<afw::table::Flag>(ctrl.name + ".flag", "success flag for " + ctrl.name))
    {}
    virtual ~ApertureFlux() {}

protected:
    afw::table::Key< afw::table::Array<double> > _fluxKey;
    afw::table::Key< afw::table::Array<double> > _errKey;
    afw::table::Key<int> _nProfileKey;
    afw::table::Key< afw::table::Flag > _flagKey;

private:
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(ApertureFlux);
};


}}}
#endif
