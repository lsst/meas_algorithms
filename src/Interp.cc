// -*- LSST-C++ -*-
/**
 * \file
 *
 * \brief Interpolate over CCD defects
 *
 * \ingroup detection
 */
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <string>
#include <typeinfo>
#include <limits>
#include <boost/format.hpp>

#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/exceptions.h>
#include <lsst/afw/image/MaskedImage.h>
#include "lsst/detection/Interp.h"

namespace image = lsst::afw::image;
namespace detection = lsst::detection;

typedef std::vector<detection::Defect::Ptr>::iterator Defect_iterT;
typedef std::vector<detection::Defect::Ptr>::const_iterator Defect_citerT;

/************************************************************************************************************/
/*
 * Classify an vector of detection::Defect::Ptr for the given row, returning a vector of 1-D
 * Defects (i.e. y0 == y1).  In general we can merge in saturated pixels at
 * this step, although we don't currently do so.
 */
static std::vector<detection::Defect::Ptr>
classify_defects(std::vector<detection::Defect::Ptr> const & badList, // list of bad things
                 int const y,           // the row to process
                 int const ncol,        // number of columns in image
                 int = 0                // number of rows in image
                ) {

    std::vector<detection::Defect::Ptr> badList1D;
    
    for (Defect_citerT begin = badList.begin(), end = badList.end(), bri = begin; bri != end; ++bri) {
        detection::Defect::Ptr defect = *bri;

        if (y < defect->getY0() || y > defect->getY1() || ncol < defect->getX0()) {
            continue;
        }

        int const x0 = defect->getX0();
        int x1 = defect->getX1();
        //
        // Look for other defects that touch this one, and push them onto badList1d
        //
        for (++bri; bri != end; ++bri) {
            defect = *bri;

            if (y < defect->getY0() || y > defect->getY1() || x1 < defect->getX0() - 1) {
                --bri;
                break;
            }
            x1 = defect->getX1();
        }
        
        int const nbad = x1 - x0 + 1;
        assert(nbad >= 1);

        defect = detection::Defect::Ptr(new detection::Defect(image::BBox(image::PointI(x0, y), nbad, 1)));
        badList1D.push_back(defect);

        if (bri == end) {
            break;
        }
    }
    //
    // Now process our new list
    //
    for (Defect_citerT begin = badList1D.begin(), end = badList1D.end(), bri = begin; bri != end; ++bri) {
        detection::Defect::Ptr defect = *bri;

        int const nbad = defect->getX1() - defect->getX0() + 1;
        assert(nbad >= 1);        

        if(defect->getX0() == 0) {
            if(nbad >= detection::Defect::WIDE_DEFECT) {
                defect->classify(detection::Defect::WIDE_LEFT, 03);
            } else {
                defect->classify(detection::Defect::LEFT, 03 << nbad);
            }
        } else if(defect->getX0() == 1) {	/* only second column is usable */
            if(nbad >= detection::Defect::WIDE_DEFECT) {
                defect->classify(detection::Defect::WIDE_NEAR_LEFT, (01 << 2) | 03);
            } else {
                defect->classify(detection::Defect::NEAR_LEFT, (01 << (nbad + 2)) | 03);
            }
        } else if(defect->getX1() == ncol - 2) { /* use only penultimate column */
            if(nbad >= detection::Defect::WIDE_DEFECT) {
                defect->classify(detection::Defect::WIDE_NEAR_RIGHT, (03 << 2) | 02);
            } else {
                defect->classify(detection::Defect::NEAR_RIGHT, (03 << (nbad + 2)) | 02);
            }
        } else if(defect->getX1() == ncol - 1) {
            if(nbad >= detection::Defect::WIDE_DEFECT) {
                defect->classify(detection::Defect::WIDE_RIGHT, 03);
            } else {
                defect->classify(detection::Defect::RIGHT, 03 << nbad);
            }
        } else if(nbad >= detection::Defect::WIDE_DEFECT) {
            defect->classify(detection::Defect::WIDE, (03 << 2) | 03);
        } else {
            defect->classify(detection::Defect::MIDDLE, (03 << (nbad + 2)) | 03);
        }
/*
 * look for bad columns in regions that we'll get `good' values from.
 *
 * We know that no two Defects are adjacent.
 */
        if(bri != begin) {
            int nshift = 0;             // number of bits to shift to get to left edge of defect pattern
            switch (defect->getPos()) {
              case detection::Defect::WIDE:		// no bits
              case detection::Defect::WIDE_NEAR_LEFT:	//       are used to encode
              case detection::Defect::WIDE_NEAR_RIGHT:	//            the bad section of data
                nshift = 0; break;
              default:
                nshift = nbad; break;
            }

            if(defect->getPos() == detection::Defect::RIGHT || defect->getPos() == detection::Defect::WIDE_RIGHT) {
                ;				/* who cares about pixels to the left? */
            } else {
                detection::Defect::Ptr const defect_m = *(bri - 1);
                assert(defect_m->getX1() < defect->getX0());

                if(defect_m->getX1() == defect->getX0() - 2) {
                    defect->classify(defect->getPos(), (defect->getType() & ~(02 << (nshift + 2))));
                }
            }
        }

        if(bri + 1 != end) {
            if(defect->getPos() == detection::Defect::LEFT || defect->getPos() == detection::Defect::WIDE_LEFT) {
                ;				/* who cares about pixels to the right? */
            } else {
                detection::Defect::Ptr const defect_p = *(bri + 1);

                if(defect->getX1() == defect_p->getX0() - 2) {
                    defect->classify(defect->getPos(), (defect->getType() & ~01));
                }
            }
        }
    }

    return badList1D;
}

/*****************************************************************************/
/*
 * Interpolate over the defects in a given line of data. In the comments,
 * a bad pixel is written as . and a good one as #, so the simple case of
 * one bad pixel is ##.## (a type of 033)
 */
template<typename MaskedImageT>
static void do_defects(std::vector<detection::Defect::Ptr> const & badList, // list of bad things
                       int const y,     // Row that we should fix
                       MaskedImageT& mi, // data to fix
                       typename MaskedImageT::Mask::Pixel const interpBit, // bit to set when we interpolated
                       int min			// minimum acceptable value
                      ) {
    typedef typename MaskedImageT::Image::Pixel ImageT;
    
    int const ncol = mi.getWidth();
    
    ImageT out1_2, out1_1, out2_1, out2_2; // == out[bad_x1-2], ..., out[bad_x2+2]
    ImageT val;                         // unpack a pixel value
    //
    // Get pointer to this row of data
    //
    typename MaskedImageT::Image::x_iterator out = mi.getImage()->row_begin(y);
    typename MaskedImageT::Mask::x_iterator mask = mi.getMask()->row_begin(y);

    for (Defect_citerT begin = badList.begin(), end = badList.end(), bri = begin; bri != end; ++bri) {
        detection::Defect::Ptr const defect = *bri;
       
        if (y < defect->getY0() || y > defect->getY1()) {
            continue;
        }

        int const bad_x0 = defect->getX0();
        int const bad_x1 = defect->getX1();

        for (int c = bad_x0; c <= bad_x1; ++c) {
            mask[c] |= interpBit;
        }

        switch (defect->getPos()) {
          case detection::Defect::LEFT:
            assert(bad_x0 >= 0 && bad_x1 + 2 < ncol);

            out2_1 = out[bad_x1 + 1];
            out2_2 = out[bad_x1 + 2];

            switch (defect->getType()) {
              case 06:              /* .##, <noise^2> = 0 */
                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 014:             /* ..##, <noise^2> = 0 */
                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 030:             /* ...##, <noise^2> = 0 */
                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 060:             /* ....##, <noise^2> = 0 */
                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 0140:            /* .....##, <noise^2> = 0 */
                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 0300:            /* ......##, <noise^2> = 0 */
                val = 0.5003*out2_1 + 0.4997*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x0 + 2] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 0600:            /* .......##, <noise^2> = 0 */
                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.5003*out2_1 + 0.4997*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x0 + 2] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x1 - 3] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 01400:           /* ........##, <noise^2> = 0 */
                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 0.5003*out2_1 + 0.4997*out2_2;
                out[bad_x0 + 2] = (val < min) ? out2_1 : val;

                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x0 + 3] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x1 - 3] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 03000:           /* .........##, <noise^2> = 0 */
                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0 + 2] = (val < min) ? out2_1 : val;

                val = 0.5003*out2_1 + 0.4997*out2_2;
                out[bad_x0 + 3] = (val < min) ? out2_1 : val;

                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x1 - 4] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x1 - 3] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              case 06000:           /* ..........##, <noise^2> = 0 */
                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0] = (val < min) ? out2_1 : val;

                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0 + 1] = (val < min) ? out2_1 : val;

                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0 + 2] = (val < min) ? out2_1 : val;

                val = 0.5000*out2_1 + 0.5000*out2_2;
                out[bad_x0 + 3] = (val < min) ? out2_1 : val;

                val = 0.5003*out2_1 + 0.4997*out2_2;
                out[bad_x0 + 4] = (val < min) ? out2_1 : val;

                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x1 - 4] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x1 - 3] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              default:
                //shFatal("Unsupported defect type: LEFT 0%o", defect->getType());
                break;			/* NOTREACHED */
            }
            break;
          case detection::Defect::WIDE_LEFT:
            assert(bad_x0 >= 0);
            if(bad_x1 + 2 >= ncol) {	/* left defect extends near
                                           right edge of data! */
                if(bad_x1 == ncol - 2) {	/* one column remains */
                    val = out[ncol - 1];
                } else {
                    val = std::numeric_limits<ImageT>::max();		/* there is no information */
                }
                for(int j = bad_x0;j <= bad_x1;j++) {
                    out[j] = val;
                }
                break;
            }
            out2_1 = out[bad_x1 + 1];
            out2_2 = out[bad_x1 + 2];

            switch (defect->getType()) {
              case 02:            /* ?#., <noise^2> = 0 */
                val = 1.0000*out2_1;
                val = (val < min) ? out2_1 : val;

                for(int j = bad_x0;j <= bad_x1;j++) {
                    out[j] = val;
                }
                break;
              case 03:            /* ?##, <noise^2> = 0 */
                val = 0.5000*out2_1 + 0.5000*out2_2;
                if (val < min) {
                    val = out2_1;
                }

                for(int j = bad_x0;j < bad_x1 - 5;j++) {
                    out[j] = val;
                }

                val = 0.5003*out2_1 + 0.4997*out2_2;
                out[bad_x1 - 5] = (val < min) ? out2_1 : val;

                val = 0.5041*out2_1 + 0.4959*out2_2;
                out[bad_x1 - 4] = (val < min) ? out2_1 : val;

                val = 0.5370*out2_1 + 0.4630*out2_2;
                out[bad_x1 - 3] = (val < min) ? out2_1 : val;

                val = 0.6968*out2_1 + 0.3032*out2_2;
                out[bad_x1 - 2] = (val < min) ? out2_1 : val;

                val = 1.0933*out2_1 - 0.0933*out2_2;
                out[bad_x1 - 1] = (val < min) ? out2_1 : val;

                val = 1.4288*out2_1 - 0.4288*out2_2;
                out[bad_x1] = (val < min) ? out2_1 : val;

                break;
              default:
                //shFatal("Unsupported defect type: WIDE_LEFT 0%o",defect[i].type);
                break;			/* NOTREACHED */
            }
	 
            break;
          case detection::Defect::RIGHT:
            assert(bad_x0 >= 2 && bad_x1 < ncol);

            out1_2 = out[bad_x0 - 2];
            out1_1 = out[bad_x0 - 1];
	 
            switch (defect->getType()) {
              case 06:              /* ##., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 014:             /* ##.., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 030:             /* ##..., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 060:             /* ##...., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 0140:            /* ##....., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x1 - 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 0300:            /* ##......, <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x0 + 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x1 - 2] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.4997*out1_2 + 0.5003*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 0600:            /* ##......., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x0 + 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x1 - 3] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x1 - 2] = (val < min) ? out1_1 : val;

                val = 0.4997*out1_2 + 0.5003*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 01400:           /* ##........, <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x0 + 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x0 + 3] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x1 - 3] = (val < min) ? out1_1 : val;

                val = 0.4997*out1_2 + 0.5003*out1_1;
                out[bad_x1 - 2] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 03000:           /* ##........., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x0 + 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x0 + 3] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x1 - 4] = (val < min) ? out1_1 : val;

                val = 0.4997*out1_2 + 0.5003*out1_1;
                out[bad_x1 - 3] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1 - 2] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              case 06000:           /* ##.........., <noise^2> = 0 */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x0 + 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x0 + 3] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x0 + 4] = (val < min) ? out1_1 : val;

                val = 0.4997*out1_2 + 0.5003*out1_1;
                out[bad_x1 - 4] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1 - 3] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1 - 2] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1 - 1] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                out[bad_x1] = (val < min) ? out1_1 : val;

                break;
              default:
                //shFatal("Unsupported defect type: RIGHT 0%o",defect[i].type);
                break;			/* NOTREACHED */
            }
            break;
          case detection::Defect::WIDE_RIGHT:
            assert(bad_x1 < ncol);

            if(bad_x0 < 2) {		/* right defect extends near
                                           left edge of data! */
                if(bad_x0 == 1) {		/* one column remains */
                    val = out[0];
                } else {
                    val = std::numeric_limits<ImageT>::max();		/* there is no information */
                }
                for(int j = bad_x0;j <= bad_x1;j++) {
                    out[j] = val;
                }
                break;
            }
	 
            out1_2 = out[bad_x0 - 2];
            out1_1 = out[bad_x0 - 1];

            switch (defect->getType()) {
              case 03:			/* ##?, S/N = infty */
                val = -0.4288*out1_2 + 1.4288*out1_1;
                out[bad_x0] = (val < min) ? out1_1 : val;

                val = -0.0933*out1_2 + 1.0933*out1_1;
                out[bad_x0 + 1] = (val < min) ? out1_1 : val;

                val = 0.3032*out1_2 + 0.6968*out1_1;
                out[bad_x0 + 2] = (val < min) ? out1_1 : val;

                val = 0.4630*out1_2 + 0.5370*out1_1;
                out[bad_x0 + 3] = (val < min) ? out1_1 : val;

                val = 0.4959*out1_2 + 0.5041*out1_1;
                out[bad_x0 + 4] = (val < min) ? out1_1 : val;

                val = 0.4997*out1_2 + 0.5003*out1_1;
                out[bad_x0 + 5] = (val < min) ? out1_1 : val;

                val = 0.5000*out1_2 + 0.5000*out1_1;
                val = (val < min) ? out1_1 : val;

                for(int j = bad_x0 + 6;j <= bad_x1;j++) {
                    out[j] = val;
                }
                break;
              default:
                //shFatal("Unsupported defect type: WIDE_RIGHT 0%o",defect[i].type);
                break;			/* NOTREACHED */
            }
            break;
          case detection::Defect::MIDDLE:
          case detection::Defect::NEAR_LEFT:
          case detection::Defect::NEAR_RIGHT:
            if(defect->getPos() == detection::Defect::MIDDLE) {
                assert(bad_x0 >= 2 && bad_x1 + 2 < ncol);
                out1_2 = out[bad_x0 - 2];
                out2_2 = out[bad_x1 + 2];
            } else if(defect->getPos() == detection::Defect::NEAR_LEFT) {
                assert(bad_x0 >= 1 && bad_x1 + 2 < ncol);
                out1_2 = -1;		/* NOTUSED */
                out2_2 = out[bad_x1 + 2];
            } else if(defect->getPos() == detection::Defect::NEAR_RIGHT) {
                assert(bad_x0 >= 2 && bad_x1 + 1 < ncol);
                out1_2 = out[bad_x0 - 2];
                out2_2 = -1;		/* NOTUSED */
            } else {
                //shFatal("Unknown defect classification %d (%s:%d)",defect->getPos(), __FILE__,__LINE__);
                out1_2 = out2_2 = -1;	/* NOTUSED */
            }
            out1_1 = out[bad_x0 - 1];
            out2_1 = out[bad_x1 + 1];
	 
            switch (defect->getType()) {
              case 012:             /* #.#., <noise^2> = 0, sigma = 1 */
                val = 0.5000*out1_1 + 0.5000*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1): val;

                break;
              case 013:		/* #.##, <noise^2> = 0 */
                val = 0.4875*out1_1 + 0.8959*out2_1 - 0.3834*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 022:             /* #..#., <noise^2> = 0, sigma = 1 */
                val = 0.7297*out1_1 + 0.2703*out2_1;
                out[bad_x0] = (val < 0) ? 0 : val;

                val = 0.2703*out1_1 + 0.7297*out2_1;
                out[bad_x1] = (val < 0) ? 0 : val;

                break;
              case 023:		/* #..##, <noise^2> = 0 */
                val = 0.7538*out1_1 + 0.5680*out2_1 - 0.3218*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3095*out1_1 + 1.2132*out2_1 - 0.5227*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 032:		/* ##.#., <noise^2> = 0 */
                val = -0.3834*out1_2 + 0.8959*out1_1 + 0.4875*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 033:		/* ##.##, <noise^2> = 0 */
                /* These coefficients are also available as
                   interp::interp_1_c1 and interp::interp_1_c2 */
                val = -0.2737*out1_2 + 0.7737*out1_1 + 0.7737*out2_1 - 0.2737*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;
                
                break;
              case 042:			/* #...#., <noise^2> = 0, sigma = 1 */
                val = 0.8430*out1_1 + 0.1570*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                val = 0.5000*out1_1 + 0.5000*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                val = 0.1570*out1_1 + 0.8430*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                break;
              case 043:             /* #...##, <noise^2> = 0 */
                val = 0.8525*out1_1 + 0.2390*out2_1 - 0.0915*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5356*out1_1 + 0.8057*out2_1 - 0.3413*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2120*out1_1 + 1.3150*out2_1 - 0.5270*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;	    
              case 062:		/* ##..#., <noise^2> = 0 */
                val = -0.5227*out1_2 + 1.2132*out1_1 + 0.3095*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.3218*out1_2 + 0.5680*out1_1 + 0.7538*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 063:		/* ##..##, <noise^2> = 0 */
                val = -0.4793*out1_2 + 1.1904*out1_1 + 0.5212*out2_1 - 0.2323*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2323*out1_2 + 0.5212*out1_1 + 1.1904*out2_1 - 0.4793*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0102:            /* #....#., <noise^2> = 0, sigma = 1 */
                val = 0.8810*out1_1 + 0.1190*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.6315*out1_1 + 0.3685*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.3685*out1_1 + 0.6315*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.1190*out1_1 + 0.8810*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                break;
              case 0103:            /* #....##, <noise^2> = 0 */
                val = 0.8779*out1_1 + 0.0945*out2_1 + 0.0276*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6327*out1_1 + 0.3779*out2_1 - 0.0106*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4006*out1_1 + 0.8914*out2_1 - 0.2920*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1757*out1_1 + 1.3403*out2_1 - 0.5160*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0142:		/* ##...#., <noise^2> = 0 */
                val = -0.5270*out1_2 + 1.3150*out1_1 + 0.2120*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.3413*out1_2 + 0.8057*out1_1 + 0.5356*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.0915*out1_2 + 0.2390*out1_1 + 0.8525*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0143:		/* ##...##, <noise^2> = 0 */
                val = -0.5230*out1_2 + 1.3163*out1_1 + 0.2536*out2_1 - 0.0469*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.3144*out1_2 + 0.8144*out1_1 + 0.8144*out2_1 - 0.3144*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.0469*out1_2 + 0.2536*out1_1 + 1.3163*out2_1 - 0.5230*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0202:            /* #.....#., <noise^2> = 0, sigma = 1 */
                val = 0.8885*out1_1 + 0.1115*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                val = 0.6748*out1_1 + 0.3252*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                val = 0.5000*out1_1 + 0.5000*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                val = 0.3252*out1_1 + 0.6748*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                val = 0.1115*out1_1 + 0.8885*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;
	    
                break;
              case 0203:            /* #.....##, <noise^2> = 0 */
                val = 0.8824*out1_1 + 0.0626*out2_1 + 0.0549*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6601*out1_1 + 0.2068*out2_1 + 0.1331*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4938*out1_1 + 0.4498*out2_1 + 0.0564*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3551*out1_1 + 0.9157*out2_1 - 0.2708*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1682*out1_1 + 1.3447*out2_1 - 0.5129*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0302:		/* ##....#., <noise^2> = 0 */
                val = -0.5160*out1_2 + 1.3403*out1_1 + 0.1757*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2920*out1_2 + 0.8914*out1_1 + 0.4006*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.0106*out1_2 + 0.3779*out1_1 + 0.6327*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0276*out1_2 + 0.0945*out1_1 + 0.8779*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0303:		/* ##....##, <noise^2> = 0 */
                val = -0.5197*out1_2 + 1.3370*out1_1 + 0.1231*out2_1 + 0.0596*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2924*out1_2 + 0.8910*out1_1 + 0.3940*out2_1 + 0.0074*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0074*out1_2 + 0.3940*out1_1 + 0.8910*out2_1 - 0.2924*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0596*out1_2 + 0.1231*out1_1 + 1.3370*out2_1 - 0.5197*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0402:            /* #......#., <noise^2> = 0, sigma = 1 */
                val = 0.8893*out1_1 + 0.1107*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.6830*out1_1 + 0.3170*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5435*out1_1 + 0.4565*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4565*out1_1 + 0.5435*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.3170*out1_1 + 0.6830*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.1107*out1_1 + 0.8893*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                break;
              case 0403:            /* #......##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0588*out2_1 + 0.0583*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6649*out1_1 + 0.1716*out2_1 + 0.1635*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5212*out1_1 + 0.2765*out2_1 + 0.2024*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4477*out1_1 + 0.4730*out2_1 + 0.0793*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3465*out1_1 + 0.9201*out2_1 - 0.2666*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0602:		/* ##.....#., <noise^2> = 0 */
                val = -0.5129*out1_2 + 1.3447*out1_1 + 0.1682*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2708*out1_2 + 0.9157*out1_1 + 0.3551*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0564*out1_2 + 0.4498*out1_1 + 0.4938*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1331*out1_2 + 0.2068*out1_1 + 0.6601*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0549*out1_2 + 0.0626*out1_1 + 0.8824*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 0603:		/* ##.....##, <noise^2> = 0 */
                val = -0.5179*out1_2 + 1.3397*out1_1 + 0.0928*out2_1 + 0.0854*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2796*out1_2 + 0.9069*out1_1 + 0.2231*out2_1 + 0.1495*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1+ out2_1) : val;

                val = 0.0533*out1_2 + 0.4467*out1_1 + 0.4467*out2_1 + 0.0533*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1+ out2_1) : val;

                val = 0.1495*out1_2 + 0.2231*out1_1 + 0.9069*out2_1 - 0.2796*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1+ out2_1) : val;

                val = 0.0854*out1_2 + 0.0928*out1_1 + 1.3397*out2_1 - 0.5179*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 01002:		/* #.......#., <noise^2> = 0, sigma = 1 */
                val = 0.8894*out1_1 + 0.1106*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.6839*out1_1 + 0.3161*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5517*out1_1 + 0.4483*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5000*out1_1 + 0.5000*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4483*out1_1 + 0.5517*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.3161*out1_1 + 0.6839*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.1106*out1_1 + 0.8894*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                break;
              case 01003:           /* #.......##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1676*out2_1 + 0.1670*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5260*out1_1 + 0.2411*out2_1 + 0.2329*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4751*out1_1 + 0.2995*out2_1 + 0.2254*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4390*out1_1 + 0.4773*out2_1 + 0.0836*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3456*out1_1 + 0.9205*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 01402:           /* ##......#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2666*out1_2 + 0.9201*out1_1 + 0.3465*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0793*out1_2 + 0.4730*out1_1 + 0.4477*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2024*out1_2 + 0.2765*out1_1 + 0.5212*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1635*out1_2 + 0.1716*out1_1 + 0.6649*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0583*out1_2 + 0.0588*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 01403:		/* ##......##, <noise^2> = 0 */
                val = -0.5177*out1_2 + 1.3400*out1_1 + 0.0891*out2_1 + 0.0886*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2771*out1_2 + 0.9095*out1_1 + 0.1878*out2_1 + 0.1797*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0677*out1_2 + 0.4614*out1_1 + 0.2725*out2_1 + 0.1984*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1+ out2_1) : val;

                val = 0.1984*out1_2 + 0.2725*out1_1 + 0.4614*out2_1 + 0.0677*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1+ out2_1) : val;

                val = 0.1797*out1_2 + 0.1878*out1_1 + 0.9095*out2_1 - 0.2771*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1+ out2_1) : val;

                val = 0.0886*out1_2 + 0.0891*out1_1 + 1.3400*out2_1 - 0.5177*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 02002:           /* #........#., <noise^2> = 0, sigma = 1 */
                val = 0.8894*out1_1 + 0.1106*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.6839*out1_1 + 0.3161*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5526*out1_1 + 0.4474*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5082*out1_1 + 0.4918*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4918*out1_1 + 0.5082*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4474*out1_1 + 0.5526*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.3161*out1_1 + 0.6839*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.1106*out1_1 + 0.8894*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                break;
              case 02003:           /* #........##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1673*out2_1 + 0.1673*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5265*out1_1 + 0.2370*out2_1 + 0.2365*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4799*out1_1 + 0.2641*out2_1 + 0.2560*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4664*out1_1 + 0.3038*out2_1 + 0.2298*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4381*out1_1 + 0.4778*out2_1 + 0.0841*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3455*out1_1 + 0.9206*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 03002:           /* ##.......#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2661*out1_2 + 0.9205*out1_1 + 0.3456*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0836*out1_2 + 0.4773*out1_1 + 0.4390*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2254*out1_2 + 0.2995*out1_1 + 0.4751*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2329*out1_2 + 0.2411*out1_1 + 0.5260*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1670*out1_2 + 0.1676*out1_1 + 0.6654*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0585*out1_2 + 0.0585*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 03003:		/* ##.......##, <noise^2> = 0 */
                val = -0.5177*out1_2 + 1.3400*out1_1 + 0.0889*out2_1 + 0.0888*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2768*out1_2 + 0.9098*out1_1 + 0.1838*out2_1 + 0.1832*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0703*out1_2 + 0.4639*out1_1 + 0.2370*out2_1 + 0.2288*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2130*out1_2 + 0.2870*out1_1 + 0.2870*out2_1 + 0.2130*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2288*out1_2 + 0.2370*out1_1 + 0.4639*out2_1 + 0.0703*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1832*out1_2 + 0.1838*out1_1 + 0.9098*out2_1 - 0.2768*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0888*out1_2 + 0.0889*out1_1 + 1.3400*out2_1 - 0.5177*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 04002:           /* #.........#., <noise^2> = 0 */
                val = 0.8894*out1_1 + 0.1106*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6839*out1_1 + 0.3161*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5527*out1_1 + 0.4473*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5091*out1_1 + 0.4909*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5000*out1_1 + 0.5000*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4909*out1_1 + 0.5091*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4473*out1_1 + 0.5527*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3161*out1_1 + 0.6839*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1106*out1_1 + 0.8894*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 04003:           /* #.........##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1673*out2_1 + 0.1673*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5265*out1_1 + 0.2368*out2_1 + 0.2367*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4804*out1_1 + 0.2601*out2_1 + 0.2595*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4712*out1_1 + 0.2685*out2_1 + 0.2603*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4654*out1_1 + 0.3043*out2_1 + 0.2302*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4380*out1_1 + 0.4778*out2_1 + 0.0842*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3455*out1_1 + 0.9206*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 06002:           /* ##........#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2661*out1_2 + 0.9206*out1_1 + 0.3455*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0841*out1_2 + 0.4778*out1_1 + 0.4381*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2298*out1_2 + 0.3038*out1_1 + 0.4664*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2560*out1_2 + 0.2641*out1_1 + 0.4799*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2365*out1_2 + 0.2370*out1_1 + 0.5265*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_2 + 0.1673*out1_1 + 0.6654*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0585*out1_2 + 0.0585*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 06003:		/* ##........##, <noise^2> = 0 */
                val = -0.5177*out1_2 + 1.3400*out1_1 + 0.0888*out2_1 + 0.0888*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2768*out1_2 + 0.9098*out1_1 + 0.1835*out2_1 + 0.1835*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0705*out1_2 + 0.4642*out1_1 + 0.2329*out2_1 + 0.2324*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2155*out1_2 + 0.2896*out1_1 + 0.2515*out2_1 + 0.2434*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2434*out1_2 + 0.2515*out1_1 + 0.2896*out2_1 + 0.2155*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2324*out1_2 + 0.2329*out1_1 + 0.4642*out2_1 + 0.0705*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1835*out1_2 + 0.1835*out1_1 + 0.9098*out2_1 - 0.2768*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0888*out1_2 + 0.0888*out1_1 + 1.3400*out2_1 - 0.5177*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 010002:          /* #..........#., <noise^2> = 0, sigma = 1 */
                val = 0.8894*out1_1 + 0.1106*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.6839*out1_1 + 0.3161*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5527*out1_1 + 0.4473*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5092*out1_1 + 0.4908*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.5009*out1_1 + 0.4991*out2_1;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4991*out1_1 + 0.5009*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4908*out1_1 + 0.5092*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.4473*out1_1 + 0.5527*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.3161*out1_1 + 0.6839*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                val = 0.1106*out1_1 + 0.8894*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1):  val;

                break;
              case 010003:          /* #..........##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1673*out2_1 + 0.1673*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5265*out1_1 + 0.2367*out2_1 + 0.2367*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4804*out1_1 + 0.2598*out2_1 + 0.2598*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4717*out1_1 + 0.2644*out2_1 + 0.2639*out2_2;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4703*out1_1 + 0.2690*out2_1 + 0.2608*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4654*out1_1 + 0.3043*out2_1 + 0.2303*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4380*out1_1 + 0.4778*out2_1 + 0.0842*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3455*out1_1 + 0.9206*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 014002:          /* ##.........#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2661*out1_2 + 0.9206*out1_1 + 0.3455*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0842*out1_2 + 0.4778*out1_1 + 0.4380*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2302*out1_2 + 0.3043*out1_1 + 0.4654*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2603*out1_2 + 0.2685*out1_1 + 0.4712*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2595*out1_2 + 0.2601*out1_1 + 0.4804*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2367*out1_2 + 0.2368*out1_1 + 0.5265*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_2 + 0.1673*out1_1 + 0.6654*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0585*out1_2 + 0.0585*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 014003:          /* ##.........##, <noise^2> = 0 */
                val = -0.5177*out1_2 + 1.3400*out1_1 + 0.0888*out2_1 + 0.0888*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2768*out1_2 + 0.9098*out1_1 + 0.1835*out2_1 + 0.1835*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0705*out1_2 + 0.4642*out1_1 + 0.2326*out2_1 + 0.2326*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2158*out1_2 + 0.2899*out1_1 + 0.2474*out2_1 + 0.2469*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2459*out1_2 + 0.2541*out1_1 + 0.2541*out2_1 + 0.2459*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2469*out1_2 + 0.2474*out1_1 + 0.2899*out2_1 + 0.2158*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2326*out1_2 + 0.2326*out1_1 + 0.4642*out2_1 + 0.0705*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1835*out1_2 + 0.1835*out1_1 + 0.9098*out2_1 - 0.2768*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0888*out1_2 + 0.0888*out1_1 + 1.3400*out2_1 - 0.5177*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 020003:          /* #...........##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1673*out2_1 + 0.1673*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5265*out1_1 + 0.2367*out2_1 + 0.2367*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4804*out1_1 + 0.2598*out2_1 + 0.2598*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4718*out1_1 + 0.2641*out2_1 + 0.2641*out2_2;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4708*out1_1 + 0.2649*out2_1 + 0.2644*out2_2;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4702*out1_1 + 0.2690*out2_1 + 0.2608*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4654*out1_1 + 0.3044*out2_1 + 0.2303*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4380*out1_1 + 0.4778*out2_1 + 0.0842*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3455*out1_1 + 0.9206*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 030002:          /* ##..........#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2661*out1_2 + 0.9206*out1_1 + 0.3455*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0842*out1_2 + 0.4778*out1_1 + 0.4380*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2303*out1_2 + 0.3043*out1_1 + 0.4654*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2608*out1_2 + 0.2690*out1_1 + 0.4703*out2_1;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2639*out1_2 + 0.2644*out1_1 + 0.4717*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2598*out1_2 + 0.2598*out1_1 + 0.4804*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2367*out1_2 + 0.2367*out1_1 + 0.5265*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_2 + 0.1673*out1_1 + 0.6654*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0585*out1_2 + 0.0585*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 030003:          /* ##..........##, <noise^2> = 0 */
                val = -0.5177*out1_2 + 1.3400*out1_1 + 0.0888*out2_1 + 0.0888*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2768*out1_2 + 0.9098*out1_1 + 0.1835*out2_1 + 0.1835*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0705*out1_2 + 0.4642*out1_1 + 0.2326*out2_1 + 0.2326*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2158*out1_2 + 0.2899*out1_1 + 0.2472*out2_1 + 0.2471*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2462*out1_2 + 0.2544*out1_1 + 0.2500*out2_1 + 0.2495*out2_2;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2495*out1_2 + 0.2500*out1_1 + 0.2544*out2_1 + 0.2462*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2471*out1_2 + 0.2472*out1_1 + 0.2899*out2_1 + 0.2158*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2326*out1_2 + 0.2326*out1_1 + 0.4642*out2_1 + 0.0705*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1835*out1_2 + 0.1835*out1_1 + 0.9098*out2_1 - 0.2768*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0888*out1_2 + 0.0888*out1_1 + 1.3400*out2_1 - 0.5177*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 040003:          /* #............##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1673*out2_1 + 0.1673*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5265*out1_1 + 0.2367*out2_1 + 0.2367*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4804*out1_1 + 0.2598*out2_1 + 0.2598*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4718*out1_1 + 0.2641*out2_1 + 0.2641*out2_2;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4708*out1_1 + 0.2646*out2_1 + 0.2646*out2_2;
                out[bad_x0 + 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4707*out1_1 + 0.2649*out2_1 + 0.2644*out2_2;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4702*out1_1 + 0.2690*out2_1 + 0.2608*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4654*out1_1 + 0.3044*out2_1 + 0.2303*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4380*out1_1 + 0.4778*out2_1 + 0.0842*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3455*out1_1 + 0.9206*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 060002:          /* ##...........#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2661*out1_2 + 0.9206*out1_1 + 0.3455*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0842*out1_2 + 0.4778*out1_1 + 0.4380*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2303*out1_2 + 0.3044*out1_1 + 0.4654*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2608*out1_2 + 0.2690*out1_1 + 0.4702*out2_1;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2644*out1_2 + 0.2649*out1_1 + 0.4708*out2_1;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2641*out1_2 + 0.2641*out1_1 + 0.4718*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2598*out1_2 + 0.2598*out1_1 + 0.4804*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2367*out1_2 + 0.2367*out1_1 + 0.5265*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_2 + 0.1673*out1_1 + 0.6654*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0585*out1_2 + 0.0585*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              default:
                //shFatal("Unsupported defect type: MIDDLE 0%o",defect[i].type);
                break;			/* NOTREACHED */
            }
            break;
          case detection::Defect::WIDE:
          case detection::Defect::WIDE_NEAR_LEFT:
          case detection::Defect::WIDE_NEAR_RIGHT:
            if(defect->getPos() == detection::Defect::WIDE_NEAR_LEFT) {
                assert(bad_x0 >= 1);

                if(bad_x1 + 2 >= ncol) {	/* left defect extends near
                                                   right edge of data! */
                    if(bad_x1 == ncol - 2) {	/* one column remains */
                        val = out[ncol - 1];
                    } else {
                        val = std::numeric_limits<ImageT>::max();		/* there is no information */
                    }
                    for(int j = bad_x0;j <= bad_x1;j++) {
                        out[j] = val;
                    }
                    break;
                }
                out1_2 = -1;		/* NOTUSED */
                out2_2 = out[bad_x1 + 2];
            } else if(defect->getPos() == detection::Defect::WIDE) {
                assert(bad_x0 >= 2 && bad_x1 + 2 < ncol);
                out1_2 = out[bad_x0 - 2];
                out2_2 = out[bad_x1 + 2];
            } else if(defect->getPos() == detection::Defect::WIDE_NEAR_RIGHT) {
                assert(bad_x1 + 1 < ncol);

                if(bad_x0 < 2) {		/* right defect extends near
                                                   left edge of data! */
                    if(bad_x0 == 1) {	/* one column remains */
                        val = out[0];
                    } else {
                        val = std::numeric_limits<ImageT>::max();	/* there is no information */
                    }
                    for(int j = bad_x0;j <= bad_x1;j++) {
                        out[j] = val;
                    }
                    break;
                }
                out1_2 = out[bad_x0 - 2];
                out2_2 = -1;		/* NOTUSED */
            } else {
                //shFatal("Unknown defect classification %d (%s:%d)",defect->getPos(), __FILE__,__LINE__);
                out1_2 = out2_2 = -1;	/* NOTUSED */
            }

            out1_1 = out[bad_x0 - 1];
            out2_1 = out[bad_x1 + 1];

            switch (defect->getType()) {

              case 06:			/* #?#., <noise^2> = 0 */
                val = 0.8894*out1_1 + 0.1106*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6839*out1_1 + 0.3161*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5527*out1_1 + 0.4473*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5092*out1_1 + 0.4908*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5010*out1_1 + 0.4990*out2_1;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5001*out1_1 + 0.4999*out2_1;
                out[bad_x0 + 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;


                val = 0.5000*out1_1 + 0.5000*out2_1;

                for(int j = bad_x0 + 6;j < bad_x1 - 5;j++) {
                    out[j] = val;
                }

                val = 0.4999*out1_1 + 0.5001*out2_1;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4990*out1_1 + 0.5010*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4908*out1_1 + 0.5092*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4473*out1_1 + 0.5527*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3161*out1_1 + 0.6839*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1106*out1_1 + 0.8894*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 07:           /* #?##, <noise^2> = 0 */
                val = 0.8829*out1_1 + 0.0585*out2_1 + 0.0585*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.6654*out1_1 + 0.1673*out2_1 + 0.1673*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.5265*out1_1 + 0.2367*out2_1 + 0.2367*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4804*out1_1 + 0.2598*out2_1 + 0.2598*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4718*out1_1 + 0.2641*out2_1 + 0.2641*out2_2;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4708*out1_1 + 0.2646*out2_1 + 0.2646*out2_2;
                out[bad_x0 + 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4707*out[bad_x0 - 1] + 0.2646*out[bad_x1 + 1] + 0.2646*out[bad_x1 + 2];

                for(int j = bad_x0 + 6;j < bad_x1 - 5;j++) {
                    out[j] = val;
                }

                val = 0.4707*out1_1 + 0.2649*out2_1 + 0.2644*out2_2;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4702*out1_1 + 0.2690*out2_1 + 0.2608*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4654*out1_1 + 0.3044*out2_1 + 0.2303*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.4380*out1_1 + 0.4778*out2_1 + 0.0842*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.3455*out1_1 + 0.9206*out2_1 - 0.2661*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_1 + 1.3452*out2_1 - 0.5125*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;
	    
                break;
              case 016:         /* ##?#., <noise^2> = 0 */
                val = -0.5125*out1_2 + 1.3452*out1_1 + 0.1673*out2_1;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2661*out1_2 + 0.9206*out1_1 + 0.3455*out2_1;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0842*out1_2 + 0.4778*out1_1 + 0.4380*out2_1;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2303*out1_2 + 0.3044*out1_1 + 0.4654*out2_1;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2608*out1_2 + 0.2690*out1_1 + 0.4702*out2_1;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2644*out1_2 + 0.2649*out1_1 + 0.4707*out2_1;
                out[bad_x0 + 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2646*out1_2 + 0.2646*out1_1 + 0.4707*out2_1;
                val = (val < min) ? 0.5*(out1_1 + out2_1) : val;
	    
                for(int j = bad_x0 + 6;j < bad_x1 - 5;j++) {
                    out[j] = val;
                }

                val = 0.2646*out1_2 + 0.2646*out1_1 + 0.4708*out2_1;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2641*out1_2 + 0.2641*out1_1 + 0.4718*out2_1;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2598*out1_2 + 0.2598*out1_1 + 0.4804*out2_1;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2367*out1_2 + 0.2367*out1_1 + 0.5265*out2_1;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1673*out1_2 + 0.1673*out1_1 + 0.6654*out2_1;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0585*out1_2 + 0.0585*out1_1 + 0.8829*out2_1;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              case 017:			/* ##?##, S/N = infty */
                val = -0.5177*out1_2 + 1.3400*out1_1 + 0.0888*out2_1 + 0.0888*out2_2;
                out[bad_x0] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = -0.2768*out1_2 + 0.9098*out1_1 + 0.1835*out2_1 + 0.1835*out2_2;
                out[bad_x0 + 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0705*out1_2 + 0.4642*out1_1 + 0.2326*out2_1 + 0.2326*out2_2;
                out[bad_x0 + 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2158*out1_2 + 0.2899*out1_1 + 0.2472*out2_1 + 0.2472*out2_2;
                out[bad_x0 + 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2462*out1_2 + 0.2544*out1_1 + 0.2497*out2_1 + 0.2497*out2_2;
                out[bad_x0 + 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2497*out1_2 + 0.2503*out1_1 + 0.2500*out2_1 + 0.2500*out2_2;
                out[bad_x0 + 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2500*out1_2 + 0.2500*out1_1 + 0.2500*out2_1 + 0.2500*out2_2;
                val = (val < min) ? 0.5*(out1_1 + out2_1) : val;
	    
                for(int j = bad_x0 + 6;j < bad_x1 - 5;j++) {
                    out[j] = val;
                }

                val = 0.2500*out1_2 + 0.2500*out1_1 + 0.2503*out2_1 + 0.2497*out2_2;
                out[bad_x1 - 5] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2497*out1_2 + 0.2497*out1_1 + 0.2544*out2_1 + 0.2462*out2_2;
                out[bad_x1 - 4] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2472*out1_2 + 0.2472*out1_1 + 0.2899*out2_1 + 0.2158*out2_2;
                out[bad_x1 - 3] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.2326*out1_2 + 0.2326*out1_1 + 0.4642*out2_1 + 0.0705*out2_2;
                out[bad_x1 - 2] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.1835*out1_2 + 0.1835*out1_1 + 0.9098*out2_1 - 0.2768*out2_2;
                out[bad_x1 - 1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                val = 0.0888*out1_2 + 0.0888*out1_1 + 1.3400*out2_1 - 0.5177*out2_2;
                out[bad_x1] = (val < min) ? 0.5*(out1_1 + out2_1) : val;

                break;
              default:
                //shFatal("Unsupported defect type: WIDE 0%o",defect[i].type);
                break;			/* NOTREACHED */
            }
            break;
        }
    }
}

/************************************************************************************************************/
/*!
 * \brief Process a set of known bad pixels in an image
 *
 * \return void
 */
namespace {
    struct Sort_DefectPtr {
        bool operator() (detection::Defect::Ptr const a, detection::Defect::Ptr const b) const {
            return a->getX0() < b->getX0();
        }
    };
}

template<typename MaskedImageT>
void lsst::detection::interpolateOverDefects(MaskedImageT& mimage, ///< Image to patch
                                             detection::PSF const &, ///< the Image's PSF
                                             std::vector<detection::Defect::Ptr> &_badList ///< List of Defects to patch
                                            ) {
/*
 * Setup desired mask planes
 */
    typename MaskedImageT::Mask::Pixel const interpBit = mimage.getMask()->getPlaneBitMask("INTRP"); // interp'd pixels
/*
 * Allow for image's origin
 */
    std::vector<detection::Defect::Ptr> badList;
    badList.reserve(_badList.size());
    for (std::vector<detection::Defect::Ptr>::iterator ptr = _badList.begin(), end = _badList.end(); ptr != end; ++ptr) {
        Defect::Ptr ndefect = Defect::Ptr(new Defect(**ptr));

        ndefect->shift(-mimage.getX0(), -mimage.getY0()); // allow for image's origin

        badList.push_back(ndefect);
    }
    
    {
        Sort_DefectPtr sort_defects;
        sort(badList.begin(), badList.end(), sort_defects);
    }
/*
 * Go through the frame looking at each pixel (except the edge ones which we ignore)
 */
    int const width = mimage.getWidth();
    int const height = mimage.getHeight();

    for (int y = 0; y != height; y++) {
        std::vector<detection::Defect::Ptr> badList1D = classify_defects(badList, y, width);
        do_defects(badList1D, y, mimage, interpBit, std::numeric_limits<typename MaskedImageT::Image::Pixel>::min());
    }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Return the interpolated value for a pixel, ignoring pixels given
 * by badmask. Interolation can either be vertical or horizontal
 *
 * Note that this is a pretty expensive routine, so use only after
 * suitable thought.
 */
template <typename MaskedImageT>
typename MaskedImageT::Image::Pixel detection::interp::singlePixel(
	int x,                          ///< column coordinate of the pixel in question
        int y,                          ///< row coordinate of the pixel in question
        MaskedImageT const& image,      ///< in this image
        bool horizontal,                ///< interpolate horizontally?
        double minval                   ///< minimum acceptable value
                                                                  ) {
#if SDSS
   BADCOLUMN defect;			/* describe a bad column */
   PIX *data;				/* temp array to interpolate in */
   int i;
   int i0, i1;				/* data corresponds to range of
					   {row,col} == [i0,i1] */
   int ndata;				/* dimension of data */
   static int ndatamax = 40;		/* largest allowable defect. XXX */
   int nrow, ncol;			/* == reg->n{row,col} */
   PIX *val;				/* pointer to pixel (rowc, colc) */
   int z1, z2;				/* range of bad {row,columns} */

   shAssert(badmask != NULL && badmask->type == shTypeGetFromName("OBJMASK"));
   shAssert(reg != NULL && reg->type == TYPE_PIX);
   nrow = reg->nrow;
   ncol = reg->ncol;

   if(horizontal) {
      for(z1 = colc - 1; z1 >= 0; z1--) {
	 if(!phPixIntersectMask(badmask, z1, rowc)) {
	    break;
	 }
      }
      z1++;
      
      for(z2 = colc + 1; z2 < ncol; z2++) {
	 if(!phPixIntersectMask(badmask, z2, rowc)) {
	    break;
	 }
      }
      z2--;
      
      i0 = (z1 > 2) ? z1 - 2 : 0;	/* origin of available required data */
      i1 = (z2 < ncol - 2) ? z2 + 2 : ncol - 1;	/* end of "      "   "    "  */

      if(i0 < 2 || i1 >= ncol - 2) {	/* interpolation will fail */
	 return(-1);			/* failure */
      }

      ndata = (i1 - i0 + 1);
      if(ndata > ndatamax) {
	 return(-1);			/* failure */
      }
      
      data = alloca(ndata*sizeof(PIX));
      for(i = i0; i <= i1; i++) {
	 data[i - i0] = reg->ROWS[rowc][i];
      }
      val = &data[colc - i0];
   } else {
      for(z1 = rowc - 1; z1 >= 0; z1--) {
	 if(!phPixIntersectMask(badmask, colc, z1)) {
	    break;
	 }
      }
      z1++;
      
      for(z2 = rowc + 1; z2 < nrow; z2++) {
	 if(!phPixIntersectMask(badmask, colc, z2)) {
	    break;
	 }
      }
      z2--;
      
      i0 = (z1 > 2) ? z1 - 2 : 0;	/* origin of available required data */
      i1 = (z2 < nrow - 2) ? z2 + 2 : nrow - 1;	/* end of "      "   "    "  */

      if(i0 < 2 || i1 >= ncol - 2) {	/* interpolation will fail */
	 return(-1);			/* failure */
      }

      ndata = (i1 - i0 + 1);
      if(ndata > ndatamax) {
	 return(-1);			/* failure */
      }
      
      data = alloca(ndata*sizeof(PIX));
      for(i = i0; i <= i1; i++) {
	 data[i - i0] = reg->ROWS[i][colc];
      }
      val = &data[rowc - i0];
   }
      
   defect.x1 = z1 - i0;
   defect.x2 = z2 - i0;
   classify_defects(&defect, 1, ndata);
   do_defect(&defect, 1, data, ndata, minval);

   return(*val);
#endif

   return std::numeric_limits<typename MaskedImageT::Image::Pixel>::min();   
}

/************************************************************************************************************/
//
// Explicit instantiations
//
typedef float imagePixelType;

template
void detection::interpolateOverDefects(image::MaskedImage<imagePixelType, image::MaskPixel> &image,
                                       detection::PSF const &psf,
                                       std::vector<detection::Defect::Ptr> &badList
                                      );
template
imagePixelType detection::interp::singlePixel(int x, int y,
                                              image::MaskedImage<imagePixelType, image::MaskPixel> const& image,
                                              bool horizontal, double minval);
//
// Why do we need double images?
//
#if 1
template
void lsst::detection::interpolateOverDefects(image::MaskedImage<double, image::MaskPixel> &image,
                                             detection::PSF const &psf,
                                             std::vector<detection::Defect::Ptr> &badList
                                            );

template
double detection::interp::singlePixel(int x, int y,
                                      image::MaskedImage<double, image::MaskPixel> const& image,
                                      bool horizontal, double minval);
#endif
