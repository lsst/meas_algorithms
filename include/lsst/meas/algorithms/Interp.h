#if !defined(LSST_DETECTION_INTERP_H)
#define LSST_DETECTION_INTERP_H
//!
// Interpolate over defects in a MaskedImage
//
#include <vector>
#include <lsst/afw/image/MaskedImage.h>
#include "lsst/detection/PSF.h"

namespace lsst { namespace detection { namespace interp {
    /**
     * LPC coefficients for sigma = 1, S/N = infty
     */
    double const lpc_1_c1 = 0.7737;
    double const lpc_1_c2 = -0.2737;
    /**
     * LPC coefficients for sigma = 1/sqrt(2), S/N = infty. These are the coeffs
     * to use when interpolating at 45degrees to the row/column
     */
    double const lpc_1s2_c1 = 0.7358;
    double const lpc_1s2_c2 = -0.2358;
    /*
     * Used to debias min(x, y)
     */
    double const min_2Gaussian_bias = -0.5641895835;	///< The mean value of the minimum of two N(0,1) variates
    
    template <typename MaskedImageT>
    typename MaskedImageT::Image::Pixel singlePixel(int x, int y, MaskedImageT const& image,
                                                    bool horizontal, double minval);
}}}

namespace lsst { namespace detection {

/**
 * \brief Encapsulate information about a bad portion of a detector
 */
class Defect {
public:
    typedef boost::shared_ptr<Defect> Ptr; //!< shared pointer to Defect
    
    enum DefectPosition {
        LEFT = 1,                       //!< defect is at left boundary
        NEAR_LEFT,			//!< defect is near left boundary
        WIDE_LEFT,			//!< defect is wide at left boundary
        MIDDLE,				//!< defect is in middle of frame
        WIDE_NEAR_LEFT,			//!< defect is near left, and wide
        WIDE,				//!< defect is in middle, and wide
        WIDE_NEAR_RIGHT,                //!< defect is near right, and wide
        NEAR_RIGHT,			//!< defect is near right boundary
        WIDE_RIGHT,			//!< defect is wide at right boundary
        RIGHT				//!< defect is at right boundary
    };

    enum { WIDE_DEFECT = 11 };          //!< minimum width of a WIDE defect
    
    explicit Defect(const lsst::afw::image::BBox& bbox = lsst::afw::image::BBox() //!< Region's bounding box
                      ) :
        _bbox(bbox), _pos(static_cast<DefectPosition>(0)), _type(0)
        {
    }
    ~Defect() {}
    

    void classify(DefectPosition pos,   //!< Position of defect in chip
                  unsigned int type     //!< Type of defect
                 ) {
        _pos = pos;    
        _type = type;
    }

    unsigned int getType() const { return _type; } //!< Return the defect's interpolation type
    DefectPosition getPos() const { return _pos; } //!< Return the position of the defect
    lsst::afw::image::BBox const & getBBox() const { return _bbox; } //!< Return the Defect's bounding box
    int const getX0() const { return _bbox.getX0(); } //!< Return the Defect's left column
    int const getX1() const { return _bbox.getX1(); } //!< Return the Defect's right column
    int const getY0() const { return _bbox.getY0(); } //!< Return the Defect's bottom row
    int const getY1() const { return _bbox.getY1(); } //!< Return the Defect's top row    

    /**
     * Offset a Defect by <tt>(dx, dy)</tt>
     */
    void shift(int dx,                 //!< How much to move defect in column direction
               int dy                  //!< How much to move in row direction
              ) {
        _bbox.shift(dx, dy);
    }
private:
    lsst::afw::image::BBox _bbox;                   //!< Bounding box for bad pixels
    DefectPosition _pos;                //!< Position of defect
    unsigned int _type;                 //!< Type of defect
};

class PSF;
    
template <typename MaskedImageT>
void interpolateOverDefects(MaskedImageT &image,
                            detection::PSF const &psf,
                            std::vector<Defect::Ptr> &badList
                           );

}}

#endif
