// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

/**
 * @file
 *
 * @brief Detect cosmic rays in a MaskedImage
 *
 * @ingroup detection
 */
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <string>
#include <typeinfo>

#include <iostream>


#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/log/Log.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/math/Random.h"
#include "lsst/meas/algorithms/CR.h"
#include "lsst/meas/algorithms/Interp.h"

/**
 * @todo These should go into afw --- actually, there're already there, but
 * in an anon namespace
 */
namespace lsst {
namespace afw {
namespace detection {
/**
 * run-length code for part of object
 */
class IdSpan {
public:
    typedef std::shared_ptr<IdSpan> Ptr;
    typedef std::shared_ptr<const IdSpan> ConstPtr;

    explicit IdSpan(int id, int y, int x0, int x1) : id(id), y(y), x0(x0), x1(x1) {}
    int id;                         /* ID for object */
    int y;                          /* Row wherein IdSpan dwells */
    int x0, x1;                     /* inclusive range of columns */
};
/**
 * comparison functor; sort by ID, then by row (y), then by column range start (x0)
 */
struct IdSpanCompar : public std::binary_function<const IdSpan::ConstPtr, const IdSpan::ConstPtr, bool> {
    bool operator()(const IdSpan::ConstPtr a, const IdSpan::ConstPtr b) {
        if (a->id < b->id) {
            return true;
        } else if(a->id > b->id) {
            return false;
        } else {
            if (a->y < b->y) {
                return true;
            } else if (a->y > b->y) {
                return false;
            } else {
                return (a->x0 < b->x0) ? true : false;
            }
        }
    }
};
/**
 * Follow a chain of aliases, returning the final resolved value.
 */
int resolve_alias(const std::vector<int>& aliases, /* list of aliases */
                  int id) {         /* alias to look up */
    int resolved = id;              /* resolved alias */

    while (id != aliases[id]) {
        resolved = id = aliases[id];
    }

    return(resolved);
}
}}} // namespace lsst::afw::detection

namespace lsst {
namespace meas {
namespace algorithms {

namespace geom = lsst::afw::geom;
namespace math = lsst::afw::math;
namespace image = lsst::afw::image;
namespace detection = lsst::afw::detection;

namespace {

template<typename ImageT, typename MaskT>
void removeCR(image::MaskedImage<ImageT, MaskT> & mi, std::vector<std::shared_ptr<detection::Footprint>> & CRs,
              double const bkgd, MaskT const , MaskT const saturBit, MaskT const badMask,
              bool const debias, bool const grow);

template<typename ImageT>
bool condition_3(ImageT *estimate, double const peak,
                 double const mean_ns, double const mean_we, double const mean_swne, double const mean_nwse,
                 double const dpeak,
                 double const dmean_ns, double const dmean_we,double const dmean_swne,double const dmean_nwse,
                 double const thresH, double const thresV, double const thresD,
                 double const cond3Fac
                );

/************************************************************************************************************/
//
// A class to hold a detected pixel
template<typename ImageT>
struct CRPixel {
    typedef typename std::shared_ptr<CRPixel> Ptr;

    CRPixel(int _col, int _row, ImageT _val, int _id = -1) :
        id(_id), col(_col), row(_row), val(_val) {
        _i = ++i;
    }
    ~CRPixel() {}

    bool operator< (const CRPixel& a) const {
        return _i < a._i;
    }

    int get_i() const {
        return _i;
    }

    int id;                             // Unique ID for cosmic ray (not cosmic ray pixel)
    int col;                            // position
    int row;                            //    of pixel
    ImageT val;                         // initial value of pixel
private:
    static int i;                       // current value of running index
    int mutable _i;                     // running index
};

template<typename ImageT>
int CRPixel<ImageT>::i = 0;

template<typename ImageT>
struct Sort_CRPixel_by_id {
    bool operator() (CRPixel<ImageT> const & a, CRPixel<ImageT> const & b) const {
        return a.id < b.id;
    }
};

/*****************************************************************************/
/*
 * This is the code to see if a given pixel is bad
 *
 * Note that the pixel in question is at index 0, so its value is pt_0[0]
 */
template <typename MaskedImageT>
bool is_cr_pixel(typename MaskedImageT::Image::Pixel *corr,      // corrected value
                 typename MaskedImageT::xy_locator loc,          // locator for this pixel
                 double const minSigma, // minSigma, or -threshold if negative
                 double const thresH, double const thresV, double const thresD, // for condition #3
                 double const bkgd,     // unsubtracted background level
                 double const cond3Fac // fiddle factor for condition #3
                )
{
    typedef typename MaskedImageT::Image::Pixel ImagePixel;
    //
    // Unpack some values
    //
    ImagePixel const v_00 = loc.image(0, 0);

    if (v_00 < 0) {
        return false;
    }
    /*
     * condition #1 is not applied on a pixel-by-pixel basis
     */
    ;
    /*
     * condition #2
     */
    ImagePixel const mean_we =   (loc.image(-1,  0) + loc.image( 1,  0))/2; // avgs of surrounding 8 pixels
    ImagePixel const mean_ns =   (loc.image( 0,  1) + loc.image( 0, -1))/2;
    ImagePixel const mean_swne = (loc.image(-1, -1) + loc.image( 1,  1))/2;
    ImagePixel const mean_nwse = (loc.image(-1,  1) + loc.image( 1, -1))/2;

    if (minSigma < 0) {         /* |thres_sky_sigma| is threshold */
        if (v_00 < -minSigma) {
            return false;
        }
    } else {
        double const thres_sky_sigma = minSigma*sqrt(loc.variance(0, 0));

        if (v_00 < mean_ns   + thres_sky_sigma &&
            v_00 < mean_we   + thres_sky_sigma &&
            v_00 < mean_swne + thres_sky_sigma &&
            v_00 < mean_nwse + thres_sky_sigma) {
            return false;
        }
    }
/*
 * condition #3
 *
 * Note that this uses mean_ns etc. even if minSigma is negative
 */
    double const dv_00 =      sqrt(loc.variance( 0,  0));
    // standard deviation of means of surrounding pixels
    double const dmean_we =   sqrt(loc.variance(-1,  0) + loc.variance( 1,  0))/2;
    double const dmean_ns =   sqrt(loc.variance( 0,  1) + loc.variance( 0, -1))/2;
    double const dmean_swne = sqrt(loc.variance(-1, -1) + loc.variance( 1,  1))/2;
    double const dmean_nwse = sqrt(loc.variance(-1,  1) + loc.variance( 1, -1))/2;

    if (!condition_3(corr,
                     v_00 - bkgd, mean_ns - bkgd, mean_we - bkgd, mean_swne - bkgd, mean_nwse - bkgd,
                     dv_00,      dmean_ns,       dmean_we,       dmean_swne,       dmean_nwse,
                     thresH, thresV, thresD, cond3Fac)){
        return false;
    }
/*
 * OK, it's a contaminated pixel
 */
    *corr += static_cast<ImagePixel>(bkgd);

    return true;
}

/************************************************************************************************************/
//
// Worker routine to process the pixels adjacent to a span (including the points just
// to the left and just to the right)
//
template <typename MaskedImageT>
void checkSpanForCRs(detection::Footprint *extras, // Extra spans get added to this Footprint
                     std::vector<CRPixel<typename MaskedImageT::Image::Pixel> >& crpixels,
                                        // a list of pixels containing CRs
                     int const y,   // the row to process
                     int const x0, int const x1, // range of pixels in the span (inclusive)
                     MaskedImageT& image, ///< Image to search
                     double const minSigma, // minSigma
                     double const thresH, double const thresV, double const thresD, // for cond. #3
                     double const bkgd, // unsubtracted background level
                     double const cond3Fac, // fiddle factor for condition #3
                     bool const keep // if true, don't remove the CRs
                    )
{
    typedef typename MaskedImageT::Image::Pixel MImagePixel;
    typename MaskedImageT::xy_locator loc = image.xy_at(x0 - 1, y); // locator for data

    int const imageX0 = image.getX0();
    int const imageY0 = image.getY0();

    for (int x = x0 - 1; x <= x1 + 1; ++x) {
        MImagePixel corr = 0;                // new value for pixel
        if (is_cr_pixel<MaskedImageT>(&corr, loc, minSigma, thresH, thresV, thresD,
                                     bkgd, cond3Fac)) {
            if (keep) {
                crpixels.push_back(CRPixel<MImagePixel>(x + imageX0, y + imageY0, loc.image()));
            }
            loc.image() = corr;

            std::vector<geom::Span> spanList(extras->getSpans()->begin(),
                                              extras->getSpans()->end());
            spanList.push_back(geom::Span(y + imageY0, x + imageX0, x +imageX0));
            extras->setSpans(std::make_shared<geom::SpanSet>(std::move(spanList)));
        }
        ++loc.x();
    }
}

/************************************************************************************************************/
/*
 * Find the sum of the pixels in a Footprint
 */
template <typename ImageT>
class CountsInCR {
public:
    CountsInCR(double const bkgd) :
        _bkgd(bkgd),
        _sum(0.0) {}

    void operator()(geom::Point2I const & point, ImageT const & val){
        _sum += val - _bkgd;
    }

    virtual void reset(detection::Footprint const&) {}
    virtual void reset() {
        _sum = 0.0;
    }

    double getCounts() const { return _sum; }
private:
    double const _bkgd;           // the Image's background level
    ImageT _sum;         // the sum of all DN in the Footprint, corrected for bkgd
};
}

/*
 * Restore all the pixels in crpixels to their pristine state
 */
template <typename ImageT>
static void reinstateCrPixels(
        ImageT *image,                                                // the image in question
        std::vector<CRPixel<typename ImageT::Pixel> > const& crpixels // a list of pixels with CRs
                             )
{
    if (crpixels.empty()) return;

    typedef typename std::vector<CRPixel<typename ImageT::Pixel> >::const_iterator crpixel_iter;
    for (crpixel_iter crp = crpixels.begin(), end = crpixels.end(); crp < end - 1 ; ++crp) {
        *image->at(crp->col - image->getX0(), crp->row - image->getY0()) = crp->val;
    }
}

/*!
 * @brief Find cosmic rays in an Image, and mask and remove them
 *
 * @return vector of CR's Footprints
 */
template <typename MaskedImageT>
std::vector<std::shared_ptr<detection::Footprint>>
findCosmicRays(MaskedImageT &mimage,      ///< Image to search
               detection::Psf const &psf, ///< the Image's PSF
               double const bkgd,         ///< unsubtracted background of frame, DN
               lsst::pex::policy::Policy const &policy, ///< Policy directing the behavior
               bool const keep                          ///< if true, don't remove the CRs
              ) {
    typedef typename MaskedImageT::Image ImageT;
    typedef typename ImageT::Pixel ImagePixel;
    typedef typename MaskedImageT::Mask::Pixel MaskPixel;

    // Parse the Policy
    double const minSigma = policy.getDouble("minSigma");    // min sigma over sky in pixel for CR candidate
    double const minDn = policy.getDouble("min_DN");         // min number of DN in an CRs
    double const cond3Fac = policy.getDouble("cond3_fac");   // fiddle factor for condition #3
    double const cond3Fac2 = policy.getDouble("cond3_fac2"); // 2nd fiddle factor for condition #3
    int const niteration = policy.getInt("niteration");      // Number of times to look for contaminated
                                                             // pixels near CRs
    int const nCrPixelMax = policy.getInt("nCrPixelMax");    // maximum number of contaminated pixels
/*
 * thresholds for 3rd condition
 *
 * Realise PSF at center of image
 */
    lsst::afw::math::Kernel::ConstPtr kernel = psf.getLocalKernel();
    if (!kernel) {
        throw LSST_EXCEPT(pexExcept::NotFoundError, "Psf is unable to return a kernel");
    }
    detection::Psf::Image psfImage = detection::Psf::Image(geom::ExtentI(kernel->getWidth(), kernel->getHeight()));
    kernel->computeImage(psfImage, true);

    int const xc = kernel->getCtrX();   // center of PSF
    int const yc = kernel->getCtrY();

    double const I0 = psfImage(xc, yc);
    double const thresH = cond3Fac2*(0.5*(psfImage(xc - 1, yc) + psfImage(xc + 1, yc)))/I0; // horizontal
    double const thresV = cond3Fac2*(0.5*(psfImage(xc, yc - 1) + psfImage(xc, yc + 1)))/I0; // vertical
    double const thresD = cond3Fac2*(0.25*(psfImage(xc - 1, yc - 1) + psfImage(xc + 1, yc + 1) +
                                           psfImage(xc - 1, yc + 1) + psfImage(xc + 1, yc - 1)))/I0; // diag
/*
 * Setup desired mask planes
 */
    MaskPixel const badBit = mimage.getMask()->getPlaneBitMask("BAD"); // Generic bad pixels
    MaskPixel const crBit = mimage.getMask()->getPlaneBitMask("CR"); // CR-contaminated pixels
    MaskPixel const interpBit = mimage.getMask()->getPlaneBitMask("INTRP"); // Interpolated pixels
    MaskPixel const saturBit = mimage.getMask()->getPlaneBitMask("SAT"); // Saturated pixels
    MaskPixel const nodataBit = mimage.getMask()->getPlaneBitMask("NO_DATA"); // Non data pixels

    MaskPixel const badMask = (badBit | interpBit | saturBit | nodataBit); // naughty pixels
/*
 * Go through the frame looking at each pixel (except the edge ones which we ignore)
 */
    int const ncol = mimage.getWidth();
    int const nrow = mimage.getHeight();

    std::vector<CRPixel<ImagePixel> > crpixels; // storage for detected CR-contaminated pixels
    typedef typename std::vector<CRPixel<ImagePixel> >::iterator crpixel_iter;
    typedef typename std::vector<CRPixel<ImagePixel> >::reverse_iterator crpixel_riter;

    for (int j = 1; j < nrow - 1; ++j) {
        typename MaskedImageT::xy_locator loc = mimage.xy_at(1, j); // locator for data

        for (int i = 1; i < ncol - 1; ++i, ++loc.x()) {
            ImagePixel corr = 0;
            if (!is_cr_pixel<MaskedImageT>(&corr, loc, minSigma,
                                           thresH, thresV, thresD, bkgd, cond3Fac)) {
                continue;
            }
/*
 * condition #4
 */
            if (loc.mask() & badMask) {
                continue;
            }
            if ((loc.mask(-1,  1) | loc.mask(0,  1) | loc.mask(1,  1) |
                 loc.mask(-1,  0) |                   loc.mask(1,  0) |
                 loc.mask(-1, -1) | loc.mask(0, -1) | loc.mask(1, -1)) & interpBit) {
                continue;
            }
/*
 * OK, it's a CR
 *
 * replace CR-contaminated pixels with reasonable values as we go through
 * image, which increases the detection rate
 */
            crpixels.push_back(CRPixel<ImagePixel>(i + mimage.getX0(), j + mimage.getY0(), loc.image()));
            loc.image() = corr;         /* just a preliminary estimate */

            if (static_cast<int>(crpixels.size()) > nCrPixelMax) {
                reinstateCrPixels(mimage.getImage().get(), crpixels);

                throw LSST_EXCEPT(lsst::pex::exceptions::LengthError,
                                  (boost::format("Too many CR pixels (max %d)") % nCrPixelMax).str());
            }
        }
    }
/*
 * We've found them on a pixel-by-pixel basis, now merge those pixels
 * into cosmic rays
 */
    std::vector<int> aliases;           // aliases for initially disjoint parts of CRs
    aliases.reserve(1 + crpixels.size()/2); // initial size of aliases

    std::vector<detection::IdSpan::Ptr> spans; // y:x0,x1 for objects
    spans.reserve(aliases.capacity());  // initial size of spans

    aliases.push_back(0);               // 0 --> 0

    /**
     In this loop, we look for strings of CRpixels on the same row and adjoining columns;
     each of these becomes a Span with a unique ID.
     */

    int ncr = 0;                        // number of detected cosmic rays
    if (!crpixels.empty()) {
        int id;                         // id number for a CR
        int x0 = -1, x1 = -1, y = -1;   // the beginning and end column, and row of this span in a CR

        // I am dummy
        CRPixel<ImagePixel> dummy(0, -1, 0, -1);
        crpixels.push_back(dummy);
        //printf("Created dummy CR: i %i, id %i, col %i, row %i, val %g\n", dummy.get_i(), dummy.id, dummy.col, dummy.row, (double)dummy.val);
        for (crpixel_iter crp = crpixels.begin(); crp < crpixels.end() - 1 ; ++crp) {
            //printf("Looking at CR: i %i, id %i, col %i, row %i, val %g\n", crp->get_i(), crp->id, crp->col, crp->row, (double)crp->val);

            if (crp->id < 0) {           // not already assigned
                crp->id = ++ncr;        // a new CR
                aliases.push_back(crp->id);
                y = crp->row;
                x0 = x1 = crp->col;
                //printf("  Assigned ID %i; looking at row %i, start col %i\n", crp->id, crp->row, crp->col);
            }
            id = crp->id;
            //printf("  Next CRpix has i=%i, id=%i, row %i, col %i\n", crp[1].get_i(), crp[1].id, crp[1].row, crp[1].col);

            if (crp[1].row == crp[0].row && crp[1].col == crp[0].col + 1) {
                //printf("  Adjoining!  Set next CRpix id = %i; x1=%i\n", crp[1].id, x1);
                crp[1].id = id;
                ++x1;
            } else {
                assert (y >= 0 && x0 >= 0 && x1 >= 0);
                spans.push_back(detection::IdSpan::Ptr(new detection::IdSpan(id, y, x0, x1)));
                //printf("  Not adjoining; adding span id=%i, y=%i, x = [%i, %i]\n", id, y, x0, x1);
            }
        }
    }

    // At the end of this loop, all crpixel entries have been assigned an ID,
    // except for the "dummy" entry at the end of the array.
    if (crpixels.size() > 0) {
        for (crpixel_iter cp = crpixels.begin(); cp != crpixels.end() - 1; cp++) {
            assert(cp->id >= 0);
            assert(cp->col >= 0);
            assert(cp->row >= 0);
        }
        // dummy:
        assert(crpixels[crpixels.size()-1].id == -1);
        assert(crpixels[crpixels.size()-1].col == 0);
        assert(crpixels[crpixels.size()-1].row == -1);
    }

    for (std::vector<detection::IdSpan::Ptr>::iterator sp = spans.begin(), end = spans.end(); sp != end; sp++) {
        assert((*sp)->id >= 0);
        assert((*sp)->y >= 0);
        assert((*sp)->x0 >= 0);
        assert((*sp)->x1 >= (*sp)->x0);
        for (std::vector<detection::IdSpan::Ptr>::iterator sp2 = sp + 1; sp2 != end; sp2++) {
            assert((*sp2)->y >= (*sp)->y);
            if ((*sp2)->y == (*sp)->y) {
                assert((*sp2)->x0 > (*sp)->x1);
            }
        }
    }

/*
 * See if spans touch each other
 */
    for (std::vector<detection::IdSpan::Ptr>::iterator sp = spans.begin(), end = spans.end();
         sp != end; ++sp) {
        int const y = (*sp)->y;
        int const x0 = (*sp)->x0;
        int const x1 = (*sp)->x1;

        // this loop will probably run for only a few steps
        for (std::vector<detection::IdSpan::Ptr>::iterator sp2 = sp + 1; sp2 != end; ++sp2) {
            if ((*sp2)->y == y) {
                // on this row (but not adjoining columns, since it would have been merged into this span);
                // keep looking.
                continue;
            } else if ((*sp2)->y != (y + 1)) {
                // sp2 is more than one row below; can't be connected.
                break;
            } else if ((*sp2)->x0 > (x1 + 1)) {
                // sp2 is more than one column away to the right; can't be connected
                break;
            } else if ((*sp2)->x1 >= (x0 - 1)) {
                // touches
                int r1 = detection::resolve_alias(aliases, (*sp)->id);
                int r2 = detection::resolve_alias(aliases, (*sp2)->id);
                aliases[r1] = r2;
            }
        }
    }





/*
 * Resolve aliases; first alias chains, then the IDs in the spans
 */
    for (unsigned int i = 0; i != spans.size(); ++i) {
        spans[i]->id = detection::resolve_alias(aliases, spans[i]->id);
    }

/*
 * Sort spans by ID, so we can sweep through them once
 */
    if (spans.size() > 0) {
        std::sort(spans.begin(), spans.end(), detection::IdSpanCompar());
    }

/*
 * Build Footprints from spans
 */
    std::vector<std::shared_ptr<detection::Footprint>> CRs; // our cosmic rays

    if (spans.size() > 0) {
        int id = spans[0]->id;
        unsigned int i0 = 0;            // initial value of i
        for (unsigned int i = i0; i <= spans.size(); ++i) { // <= size to catch the last object
            if (i == spans.size() || spans[i]->id != id) {
                std::shared_ptr<detection::Footprint> cr(new detection::Footprint());

                std::vector<geom::Span> spanList;
                spanList.reserve(i - i0);
                for (; i0 < i; ++i0) {
                    spanList.push_back(geom::Span(spans[i0]->y, spans[i0]->x0, spans[i0]->x1));
                }
                cr->setSpans(std::make_shared<geom::SpanSet>(std::move(spanList)));
                CRs.push_back(cr);
            }

            if (i < spans.size()) {
                id = spans[i]->id;
            }
        }
    }

    reinstateCrPixels(mimage.getImage().get(), crpixels);
/*
 * apply condition #1
 */
    CountsInCR<typename ImageT::Pixel> CountDN(bkgd);
    for (std::vector<std::shared_ptr<detection::Footprint>>::iterator cr = CRs.begin(), end = CRs.end();
         cr != end; ++cr) {
        // find the sum of pixel values within the CR
        (*cr)->getSpans()->applyFunctor(CountDN, *mimage.getImage());

        LOGL_DEBUG("TRACE4.algorithms.CR", "CR at (%d, %d) has %g DN",
                   (*cr)->getBBox().getMinX(), (*cr)->getBBox().getMinY(), CountDN.getCounts());
        if (CountDN.getCounts() < minDn) { /* not bright enough */
            LOGL_DEBUG("TRACE5.algorithms.CR", "Erasing CR");

            cr = CRs.erase(cr);
            --cr;                       // back up to previous CR (we're going to increment it)
            --end;
        }
        CountDN.reset();
    }
    ncr = CRs.size();           /* some may have been too faint */
/*
 * We've found them all, time to kill them all
 */
    bool const debias_values = true;
    bool grow = false;
    LOGL_DEBUG("TRACE2.algorithms.CR", "Removing initial list of CRs");
    removeCR(mimage, CRs, bkgd, crBit, saturBit, badMask, debias_values, grow);
#if 0                                   // Useful to see phase 2 in ds9; debugging only
    (void)setMaskFromFootprintList(mimage.getMask().get(), CRs,
                                   mimage.getMask()->getPlaneBitMask("DETECTED"));
#endif
/*
 * Now that we've removed them, go through image again, examining area around
 * each CR for extra bad pixels. Note that we set cond3Fac = 0 for this pass
 *
 * We iterate niteration times;  niter==1 was sufficient for SDSS data, but megacam
 * CCDs are different -- who knows for other devices?
 */
    bool too_many_crs = false;          // we've seen too many CR pixels
    int nextra = 0;                     // number of pixels added to list of CRs
    for (int i = 0; i != niteration && !too_many_crs; ++i) {
        LOGL_DEBUG("TRACE1.algorithms.CR", "Starting iteration %d", i);
        for (std::vector<std::shared_ptr<detection::Footprint>>::iterator fiter = CRs.begin();
             fiter != CRs.end(); fiter++) {
            std::shared_ptr<detection::Footprint> cr = *fiter;
/*
 * Are all those `CR' pixels interpolated?  If so, don't grow it
 */
            {
                // this work should be taken on in DM-9538
                //std::shared_ptr<detection::Footprint> om = footprintAndMask(cr, mimage.getMask(), interpBit);
                // Placeholder until DM-9538 then remove empty footprint
                auto om = std::make_shared<detection::Footprint>();
                int const npix = (om) ? om->getArea() : 0;

                if (static_cast<std::size_t>(npix) == cr->getArea()) {
                    continue;
                }
            }
/*
 * No; some of the suspect pixels aren't interpolated
 */
            detection::Footprint extra;                     // extra pixels added to cr
            for (auto siter = cr->getSpans()->begin();
                 siter != cr->getSpans()->end(); siter++) {
                auto const span = siter;

                /*
                 * Check the lines above and below the span.  We're going to check a 3x3 region around
                 * the pixels, so we need a buffer around the edge.  We check the pixels just to the
                 * left/right of the span, so the buffer needs to be 2 pixels (not just 1) in the
                 * column direction, but only 1 in the row direction.
                 */
                int const y = span->getY() - mimage.getY0();
                if (y < 2 || y >= nrow - 2) {
                    continue;
                }
                int x0 = span->getX0() - mimage.getX0();
                int x1 = span->getX1() - mimage.getX0();
                x0 = (x0 < 2) ? 2 : (x0 > ncol - 3) ? ncol - 3 : x0;
                x1 = (x1 < 2) ? 2 : (x1 > ncol - 3) ? ncol - 3 : x1;

                checkSpanForCRs(&extra, crpixels, y - 1, x0, x1, mimage,
                                minSigma/2, thresH, thresV, thresD, bkgd, 0, keep);
                checkSpanForCRs(&extra, crpixels, y,     x0, x1, mimage,
                                minSigma/2, thresH, thresV, thresD, bkgd, 0, keep);
                checkSpanForCRs(&extra, crpixels, y + 1, x0, x1, mimage,
                                minSigma/2, thresH, thresV, thresD, bkgd, 0, keep);
            }

            if (extra.getSpans()->size() > 0) {      // we added some pixels
                if (nextra + static_cast<int>(crpixels.size()) > nCrPixelMax) {
                    too_many_crs = true;
                    break;
                }

                nextra += extra.getArea();

                std::vector<geom::Span> tmpSpanList(cr->getSpans()->begin(),
                                                    cr->getSpans()->end());
                for (auto const & spn : (*extra.getSpans())) {
                    tmpSpanList.push_back(spn);
                }
                cr->setSpans(std::make_shared<geom::SpanSet>(std::move(tmpSpanList)));
            }
        }

        if (nextra == 0) {
            break;
        }
    }
/*
 * mark those pixels as CRs
 */
    if (!too_many_crs) {
        for (auto const & foot : CRs) {
            foot->getSpans()->setMask(*mimage.getMask().get(), crBit);
        }
    }
/*
 * Maybe reinstate initial values; n.b. the same pixel may appear twice, so we want the
 * first value stored (hence the uses of rbegin/rend)
 *
 * We have to do this if we decide _not_ to remove certain CRs,
 * for example those which lie next to saturated pixels
 */
    if (keep || too_many_crs) {
        if (crpixels.size() > 0) {
            int const imageX0 = mimage.getX0();
            int const imageY0 = mimage.getY0();

            std::sort(crpixels.begin(), crpixels.end()); // sort into birth order

            crpixel_riter rend = crpixels.rend();
            for (crpixel_riter crp = crpixels.rbegin(); crp != rend; ++crp) {
                if (crp->row == -1)
                    // dummy; skip it.
                    continue;
                mimage.at(crp->col - imageX0, crp->row - imageY0).image() = crp->val;
            }
        }
    } else {
        if (true || nextra > 0) {
            grow = true;
            LOGL_DEBUG("TRACE2.algorithms.CR", "Removing final list of CRs, grow = %d", grow);
            removeCR(mimage, CRs, bkgd, crBit, saturBit, badMask, debias_values, grow);
        }
/*
 * we interpolated over all CR pixels, so set the interp bits too
 */
        for (auto const & foot : CRs) {
            foot->getSpans()->setMask(*mimage.getMask().get(), static_cast<MaskPixel>(crBit | interpBit));
        }
    }

    if (too_many_crs) {                 // we've cleaned up, so we can throw the exception
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError,
                          (boost::format("Too many CR pixels (max %d)") % nCrPixelMax).str());
    }

    return CRs;
}

/*****************************************************************************/
namespace {
/*
 * Is condition 3 true?
 */
template<typename ImageT>
bool condition_3(ImageT *estimate,        // estimate of true value of pixel
                 double const peak,       // counts in central pixel (no sky)
                 double const mean_ns,    // mean in NS direction (no sky)
                 double const mean_we,    //  "   "  WE    "  "     "   "
                 double const mean_swne,  //  "   "  SW-NE "  "     "   "
                 double const mean_nwse,  //  "   "  NW-SE "  "     "   "
                 double const dpeak,      // standard deviation of peak value
                 double const dmean_ns,   //   s.d. of mean in NS direction
                 double const dmean_we,   //    "   "   "   "  WE    "  "
                 double const dmean_swne, //  "   "   "   "  SW-NE "  "
                 double const dmean_nwse, //  "   "   "   "  NW-SE "  "
                 double const thresH,     // horizontal threshold
                 double const thresV,     // vertical threshold
                 double const thresD,     // diagonal threshold
                 double const cond3Fac    // fiddle factor for noise
                )
{
   if (thresV*(peak - cond3Fac*dpeak) > mean_ns + cond3Fac*dmean_ns) {
       *estimate = (ImageT)mean_ns;
       return true;
   }

   if (thresH*(peak - cond3Fac*dpeak) > mean_we + cond3Fac*dmean_we) {
       *estimate = mean_we;
       return true;
   }

   if (thresD*(peak - cond3Fac*dpeak) > mean_swne + cond3Fac*dmean_swne) {
       *estimate = mean_swne;
       return true;
   }

   if (thresD*(peak - cond3Fac*dpeak) > mean_nwse + cond3Fac*dmean_nwse) {
       *estimate = mean_nwse;
       return true;
   }

   return false;
}

/************************************************************************************************************/
/*
 * Interpolate over a CR's pixels
 */
template <typename MaskedImageT>
class RemoveCR {
public:
    RemoveCR(MaskedImageT const& mimage,
             double const bkgd,
             typename MaskedImageT::Mask::Pixel badMask,
             bool const debias,
             lsst::afw::math::Random& rand
            ) :_image(mimage),
                _bkgd(bkgd),
                _ncol(mimage.getWidth()),
                _nrow(mimage.getHeight()),
                _badMask(badMask),
                _debias(debias),
                _rand(rand) {}

    void operator()(geom::Point2I const & point, int) {
        int const & x = point.getX();
        int const & y = point.getY();
        typename MaskedImageT::xy_locator loc = _image.xy_at(x - _image.getX0(),
                                                             y - _image.getY0());
        typedef typename MaskedImageT::Image::Pixel MImagePixel;
        MImagePixel min = std::numeric_limits<MImagePixel>::max();
        int ngood = 0;          // number of good values on min

        MImagePixel const minval = _bkgd - 2*sqrt(loc.variance()); // min. acceptable pixel value after interp
/*
 * W-E row
 */
        if (x - 2 >= 0 && x + 2 < _ncol) {
            if ((loc.mask(-2, 0) | _badMask) || (loc.mask(-1, 0) | _badMask) ||
                (loc.mask( 1, 0) | _badMask) || (loc.mask( 2, 0) | _badMask)) {
                ;                       // estimate is contaminated
            } else {
                MImagePixel const v_m2 = loc.image(-2, 0);
                MImagePixel const v_m1 = loc.image(-1, 0);
                MImagePixel const v_p1 = loc.image( 1, 0);
                MImagePixel const v_p2 = loc.image( 2, 0);

                MImagePixel const tmp =
                    interp::lpc_1_c1*(v_m1 + v_p1) + interp::lpc_1_c2*(v_m2 + v_p2);

                if (tmp > minval && tmp < min) {
                    min = tmp;
                    ngood++;
                }
            }
        }
/*
 * N-S column
 */
        if (y - 2 >= 0 && y + 2 < _nrow) {
            if ((loc.mask(0, -2) | _badMask) || (loc.mask(0, -1) | _badMask) ||
                (loc.mask(0,  1) | _badMask) || (loc.mask(0,  2) | _badMask)) {
                ;                       /* estimate is contaminated */
            } else {
                MImagePixel const v_m2 = loc.image(0, -2);
                MImagePixel const v_m1 = loc.image(0, -1);
                MImagePixel const v_p1 = loc.image(0,  1);
                MImagePixel const v_p2 = loc.image(0,  2);

                MImagePixel const tmp =
                    interp::lpc_1_c1*(v_m1 + v_p1) + interp::lpc_1_c2*(v_m2 + v_p2);

                if (tmp > minval && tmp < min) {
                    min = tmp;
                    ngood++;
                }
            }
        }
/*
 * SW--NE diagonal
 */
        if (x - 2 >= 0 && x + 2 < _ncol && y - 2 >= 0 && y + 2 < _nrow) {
            if ((loc.mask(-2, -2) | _badMask) || (loc.mask(-1, -1) | _badMask) ||
                (loc.mask( 1,  1) | _badMask) || (loc.mask( 2,  2) | _badMask)) {
                ;                       /* estimate is contaminated */
            } else {
                MImagePixel const v_m2 = loc.image(-2, -2);
                MImagePixel const v_m1 = loc.image(-1, -1);
                MImagePixel const v_p1 = loc.image( 1,  1);
                MImagePixel const v_p2 = loc.image( 2,  2);

                MImagePixel const tmp =
                    interp::lpc_1s2_c1*(v_m1 + v_p1) + interp::lpc_1s2_c2*(v_m2 + v_p2);

                if (tmp > minval && tmp < min) {
                    min = tmp;
                    ngood++;
                }
            }
        }
/*
 * SE--NW diagonal
 */
        if (x - 2 >= 0 && x + 2 < _ncol && y - 2 >= 0 && y + 2 < _nrow) {
            if ((loc.mask( 2, -2) | _badMask) || (loc.mask( 1, -1) | _badMask) ||
                (loc.mask(-1,  1) | _badMask) || (loc.mask(-2,  2) | _badMask)) {
                ;                       /* estimate is contaminated */
            } else {
                MImagePixel const v_m2 = loc.image( 2, -2);
                MImagePixel const v_m1 = loc.image( 1, -1);
                MImagePixel const v_p1 = loc.image(-1,  1);
                MImagePixel const v_p2 = loc.image(-2,  2);

                MImagePixel const tmp =
                    interp::lpc_1s2_c1*(v_m1 + v_p1) + interp::lpc_1s2_c2*(v_m2 + v_p2);

                if (tmp > minval && tmp < min) {
                    min = tmp;
                    ngood++;
                }
            }
        }
/*
 * Have we altogether failed to find an acceptable value? If so interpolate
 * using the full-up interpolation code both vertically and horizontally
 * and take the average. This can fail for large enough defects (e.g. CRs
 * lying in bad columns), in which case the interpolator returns -1. If
 * both directions fail, use the background value.
 */
        if (ngood == 0) {
            std::pair<bool, MImagePixel const> val_h =
                interp::singlePixel(x, y, _image, true,  minval);
            std::pair<bool, MImagePixel const> val_v =
                interp::singlePixel(x, y, _image, false, minval);

            if (!val_h.first) {
                if (!val_v.first) {    // Still no good value. Guess wildly
                    min = _bkgd + sqrt(loc.variance())*_rand.gaussian();
                } else {
                    min = val_v.second;
                }
            } else {
                if (val_v.first) {
                    min = val_h.second;
                } else {
                    min = (val_v.second + val_h.second)/2;
                }
            }
        }
/*
 * debias the minimum; If more than one uncontaminated estimate was
 * available, estimate the bias to be simply that due to choosing the
 * minimum of two Gaussians. In fact, even some of the "good" pixels
 * may have some extra charge, so even if ngood > 2, still use this
 * estimate
 */
        if (ngood > 0) {
            LOGL_DEBUG("TRACE3.algorithms.CR", "Adopted min==%g at (%d, %d) (ngood=%d)",
                                  static_cast<double>(min), x, y, ngood);
        }

        if (_debias && ngood > 1) {
            min -= interp::min2GaussianBias*sqrt(loc.variance())*_rand.gaussian();
        }

        loc.image() = min;
    }
private:
    MaskedImageT const & _image;
    double _bkgd;
    int _ncol, _nrow;
    typename MaskedImageT::Mask::Pixel _badMask;
    bool _debias;
    lsst::afw::math::Random& _rand;
};

/************************************************************************************************************/
/*
 * actually remove CRs from the frame
 */
template<typename ImageT, typename MaskT>
void removeCR(image::MaskedImage<ImageT, MaskT> & mi,  // image to search
              std::vector<std::shared_ptr<detection::Footprint>> & CRs, // list of cosmic rays
              double const bkgd, // non-subtracted background
              MaskT const , // Bit value used to label CRs
              MaskT const saturBit, // Bit value used to label saturated pixels
              MaskT const badMask, // Bit mask for bad pixels
              bool const debias, // statistically debias values?
              bool const grow   // Grow CRs?
             )
{
    lsst::afw::math::Random rand;    // a random number generator
    /*
     * replace the values of cosmic-ray contaminated pixels with 1-dim 2nd-order weighted means Cosmic-ray
     * contaminated pixels have already been given a mask value, crBit
     *
     * If there are no good options (i.e. all estimates are contaminated), try using just pixels that are not
     * CRs; failing that, interpolate in the row- or column direction over as large a distance as is required
     *
     * XXX SDSS (and we) go through this list backwards; why?
     */

    // a functor to remove a CR
    RemoveCR<image::MaskedImage<ImageT, MaskT> > removeCR(mi, bkgd, badMask, debias, rand);

    for (std::vector<std::shared_ptr<detection::Footprint>>::reverse_iterator fiter = CRs.rbegin();
         fiter != CRs.rend(); ++fiter) {
        std::shared_ptr<detection::Footprint> cr = *fiter;
/*
 * If I grow this CR does it touch saturated pixels?  If so, don't
 * interpolate and add CR pixels to saturated mask
 */
/* This is commented out pending work from DM-9538
        if (grow && cr->getArea() < 100) {
            try {
                bool const isotropic = false; // use a slow isotropic grow?
                auto gcr = std::make_shared<detection::Footprint>(cr->getSpans()->dilate(1), cr.getRegion());
                std::shared_ptr<detection::Footprint> const saturPixels = footprintAndMask(gcr, mi.getMask(), saturBit);
                auto const saturPixels = footprintAndMask(gcr, mi.getMask(), saturBit);

             if (saturPixels->getArea() > 0) { // pixel is adjacent to a saturation trail
                 setMaskFromFootprint(mi.getMask().get(), *saturPixels, saturBit);

                    continue;
                }
            } catch(lsst::pex::exceptions::LengthError &) {
                continue;
            }
        }
        */
/*
 * OK, fix it
 */
        // applyFunctor must have an argument other than the functor, however in this case
        // additional arguments are unneeded, so pass a constant 1 to be iterated over
        // and ignore
        cr->getSpans()->applyFunctor(removeCR, 1);
    }
}
}

/************************************************************************************************************/
//
// Explicit instantiations
// \cond
#define INSTANTIATE(TYPE) \
    template \
    std::vector<std::shared_ptr<detection::Footprint>> \
    findCosmicRays(lsst::afw::image::MaskedImage<TYPE> &image,  \
                   detection::Psf const &psf,                   \
                   double const bkgd,                           \
                   lsst::pex::policy::Policy const& policy,     \
                   bool const keep                              \
                  )

INSTANTIATE(float);
INSTANTIATE(double);                    // Why do we need double images?
// \endcond
}}} // namespace lsst::meas::algorithms
