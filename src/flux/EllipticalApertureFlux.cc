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
 
#if 0 && defined(__ICC)
#pragma warning (push)
#pragma warning (disable: 21)           // type qualifiers are meaningless in this declaration
#pragma warning disable: 68)            // integer conversion resulted in a change of sign
#pragma warning (disable: 279)          // controlling expression is constant
#pragma warning (disable: 304)          // access control not specified ("public" by default)
#pragma warning (disable: 444)          // destructor for base class ... is not virtual
//#pragma warning (pop)
#endif

#include <cmath>
#include <limits>
#include <numeric>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/FluxControl.h"
#include "lsst/meas/algorithms/ApertureFlux.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * Implement "EllipticalAperture" photometry.
 * @brief A class that knows how to calculate fluxes as a simple sum over a Footprint
 */
class EllipticalApertureFlux : public ApertureFlux {
public:
    EllipticalApertureFlux(EllipticalApertureFluxControl const & ctrl, afw::table::Schema & schema) :
        ApertureFlux(ctrl, schema)
        {}

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(EllipticalApertureFlux);
};

/************************************************************************************************************/
/**
 * @brief Given an image and a source position, calculate a set of fluxes in elliptical apertures
 */
template <typename PixelT>
void EllipticalApertureFlux::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_flagKey, true);         // say we've failed so that's the result if we throw
    source.set(_nProfileKey, 0);        // no points measured

    if (source.getShapeFlag()) {        // the shape's bad; give up now
        return;
    }
    afw::geom::ellipses::Axes const shape = source.getShape();
    double const ellip = 1.0 - shape.getB()/shape.getA();
    double const theta = shape.getTheta();
  
    VectorD const & radii = static_cast<EllipticalApertureFluxControl const &>(getControl()).radii;
    int const nradii = radii.size();

    typename afw::image::Exposure<PixelT>::MaskedImageT const& mimage = exposure.getMaskedImage();
    double const xcen = center.getX();   // object's column position
    double const ycen = center.getY();   // object's row position

    double oradius = 0.0;                // old size of shape
    double const fac = ::sqrt(shape.getA()/shape.getB()); // calculateSincApertureFlux expects the major axis
    for (int i = 0; i != nradii; ++i) {
        double const radius = fac*radii[i];

        std::pair<double, double> flux =
            algorithms::photometry::calculateSincApertureFlux(mimage, xcen, ycen,
                                                              oradius, radius, theta, ellip);
        oradius = radius;

        source.set(_fluxKey[i], flux.first);
        source.set(_errKey[i],  flux.second);
    }
    source.set(_nProfileKey, nradii);
    source.set(_flagKey, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(EllipticalApertureFlux);

PTR(AlgorithmControl) EllipticalApertureFluxControl::_clone() const {
    return boost::make_shared<EllipticalApertureFluxControl>(*this);
}

PTR(Algorithm) EllipticalApertureFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    if (metadata) {
        std::string key = this->name + ".radii";
        std::replace(key.begin(), key.end(), '.', '_');
        metadata->add(key, radii, "Radii for aperture flux measurement");
    }
    return boost::make_shared<EllipticalApertureFlux>(*this, boost::ref(schema));
}

}}}
