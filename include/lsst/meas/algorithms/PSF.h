#if !defined(LSST_DETECTION_PSF_H)
#define LSST_DETECTION_PSF_H
//!
// Describe an image's PSF
//
#include <boost/shared_ptr.hpp>
#include <lsst/daf/data.h>

namespace lsst { namespace detection {

/*!
 * \brief Represent an image's PSF
 */
class PSF : public lsst::daf::data::LsstBase {
public:
    typedef boost::shared_ptr<PSF> Ptr;

    explicit PSF();
    virtual ~PSF() {}

    virtual std::string toString() const = 0;

    void setA(double const A) { _A = A; } ///< Set the central amplitude
    double getA() const { return _A; }  ///< Get the central amplitude
    
    virtual double getValue(double const col = 0, double const row = 0) const = 0; ///< Evaluate the PSF at (col, row)
private:
    double _A;                          ///< Central amplitude of PSF
};

/*!
 * \brief Represent a PSF as a circularly symmetrical double Gaussian
 *
 * @note This is just a place holder
 */
class dgPSF : public PSF {
public:
    typedef boost::shared_ptr<dgPSF> Ptr;

    explicit dgPSF() : PSF() {}
    explicit dgPSF(
                   double sigma1,         ///< Width of inner Gaussian
                   double sigma2 = 0,    ///< Width of outer Gaussian
                   double b = 0           ///< Central amplitude of outer Gaussian (inner amplitude == 1)
                  );

    std::string toString() const;

    double getValue(double const col = 0, double const row = 0) const;
private:
    double _sigma1;                     ///< Width of inner Gaussian
    double _sigma2;                     ///< Width of outer Gaussian
    double _b;                          ///< Central amplitude of outer Gaussian (inner amplitude == 1)
};

}}

#endif
