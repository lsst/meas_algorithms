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
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/PhotometryControl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {
/**
 * Implement "Naive" photometry.
 * @brief A class that knows how to calculate photometrys as a simple sum over a Footprint
 */
template<typename ExposureT>
class NaivePhotometer : public Algorithm<afwDetection::Photometry, ExposureT>
{
public:
    typedef Algorithm<afwDetection::Photometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<NaivePhotometer> Ptr;
    typedef boost::shared_ptr<NaivePhotometer const> ConstPtr;

    explicit NaivePhotometer(NaivePhotometryControl const & ctrl) : AlgorithmT(), _radius(ctrl.radius) {}

    /// Modifier
    void setRadius(double radius) { _radius = radius; }

    /// Accessor
    double getRadius() const { return _radius; }

    virtual std::string getName() const { return "NAIVE"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<NaivePhotometer<ExposureT> >(*this);
    }

    virtual PTR(afwDetection::Photometry) measureSingle(
        afwDetection::Source const&,
        afwDetection::Source const&,
        ExposurePatch<ExposureT> const&
        ) const;

private:
    double _radius;
};

template <typename MaskedImageT>
class FootprintFlux : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0.0), _sumVar(0.0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = _sumVar = 0.0;
    }
    void reset(afwDetection::Footprint const&) {}        

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int,                                   ///< column-position of pixel
                    int                                    ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = loc.image(0, 0);
        typename MaskedImageT::Variance::Pixel vval = loc.variance(0, 0);
        _sum += ival;
        _sumVar += vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

    /// Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:
    double _sum;
    double _sumVar;
};

template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage,          ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0.0), _sumVar(0.0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset(afwDetection::Footprint const& foot) {
        _sum = _sumVar = 0.0;
        
        afwGeom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size for "
                                             "%d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    void reset() {}
    
    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        typename MaskedImageT::Variance::Pixel vval = iloc.variance(0, 0);
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
        _sumVar += wval*wval*vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the variance in the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;                                   // The variance of our desired sum
    int _x0, _y0;                                     // the origin of the current Footprint
};

            
/*****************************************************************************************************/
/**
 * Accumulate sum(x) and sum(x**2)
 */
template<typename T>
struct getSum2 {
    getSum2() : sum(0.0), sum2(0.0) {}
    
    getSum2& operator+(T x) {
        sum += x;
        sum2 += x*x;
        return *this;
    }
    
    double sum;                         // \sum_i(x_i)
    double sum2;                        // \sum_i(x_i^2)
};


/************************************************************************************************************/
/**
 * @brief Given an image and a pixel position, return a Photometry
 */
template<typename ExposureT>
PTR(afwDetection::Photometry) NaivePhotometer<ExposureT>::measureSingle(
    afwDetection::Source const& target,
    afwDetection::Source const& source,
    ExposurePatch<ExposureT> const& patch
    ) const
{
    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;

    CONST_PTR(ExposureT) exposure = patch.getExposure();
    MaskedImageT const& mimage = exposure->getMaskedImage();

    double const xcen = patch.getCenter().getX();   ///< object's column position
    double const ycen = patch.getCenter().getY();   ///< object's row position

    int const ixcen = afwImage::positionToIndex(xcen);
    int const iycen = afwImage::positionToIndex(ycen);

    // BBox for data image
    afwGeom::BoxI imageBBox(mimage.getBBox(afwImage::PARENT));

    /* ******************************************************* */
    // Aperture photometry
    FootprintFlux<typename ExposureT::MaskedImageT> fluxFunctor(mimage);
    afwDetection::Footprint const foot(
        afwGeom::PointI(ixcen, iycen), 
        getRadius(), 
        imageBBox
        );
    fluxFunctor.apply(foot);

    double aperFlux = fluxFunctor.getSum();

    double aperFluxErr = ::sqrt(fluxFunctor.getSumVar());
    return boost::make_shared<afwDetection::Photometry>(aperFlux, aperFluxErr);
}

} // anonymous

LSST_ALGORITHM_CONTROL_PRIVATE_IMPL(NaivePhotometryControl, NaivePhotometer)

}}}
