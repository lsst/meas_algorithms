#ifndef LSST_MEAS_ALGORITHMS_SECOND_MOMENTS_H
#define LSST_MEAS_ALGORITHMS_SECOND_MOMENTS_H

/// \file

#include "lsst/afw/detection.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"

namespace lsst{
namespace meas{
namespace algorithms{

/************************************************************************************************************/
/**
 * \brief Calculate a detected source's second moments
 */
template <typename MaskedImageT>
class FootprintSecondMoments : public lsst::afw::detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintSecondMoments(MaskedImageT const& mimage,      ///< The image the source lives in
                        double x,                           ///< x position of source
                        double y)                           ///< y position of source
        : detection::FootprintFunctor<MaskedImageT>(mimage), 
            _x(x), _y(y), 
            _n(0), _sum(0), 
            _sumXx(0), _sumYy(0), _sumXy(0) {
    }

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc,  ///< locator pointing at the pixel
                    int x,                                  ///< column-position of pixel
                    int y                                   ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image();

        _n++;
        _sum += val;
        double xdev = _x - lsst::afw::image::indexToPosition(x);
        double ydev = _y - lsst::afw::image::indexToPosition(y);
        
        _sumXx += xdev*xdev*val;
        _sumYy += ydev*ydev*val;
        _sumXy += xdev*ydev*val;
    }

    /// Return the number of pixels
    int getN() const { return _n; }
    /// Return the Footprint's flux
    double getSum() const { return _sum; } 
    /// Get centroid's x value
    double getX()  const { return _x; }
    /// Get centroid's y value
    double getY()  const { return _y; }  
    /// Get the second moment of x about x
    double getXx() const { return _sumXx/_sum; }   
    /// Get the second moment of y about y
    double getYy() const { return _sumYy/_sum; }    
    /// Get the second moment of x about y
    double getXy() const { return _sumXy/_sum; }
    
private:
    int _n;
    double _x, _y;
    double _sum, _sumXx, _sumYy, sumXy
};

}}} //namespace lsst::meas::algorithms

#endif //LSST_MEAS_ALGORITHMS_SECOND_MOMENTS_H
