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

#include "lsst/meas/algorithms/shapelet/Pixel.h"
#include "lsst/meas/algorithms/shapelet/Params.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

#if 0
    void getPixList(
        const Image<double>& im, PixelList& pix,
        const Position cen, double sky, double noise, double gain,
        const Image<double>* weightImage, const Transformation& trans,
        double aperture, long& flag)
    {
        xdbg<<"Start GetPixList\n";
        if (weightImage) {
            xdbg<<"Using weight image for pixel noise.\n";
        } else {
            xdbg<<"noise = "<<noise<<std::endl;
            xdbg<<"gain = "<<gain<<std::endl;
        }

        DSmallMatrix22 localD;
        trans.getDistortion(cen,localD);

        double det = std::abs(localD.TMV_det());
        double pixScale = sqrt(det); // arcsec/pixel
        xdbg<<"pixscale = "<<pixScale<<std::endl;

        // xAp,yAp are the maximum deviation from the center in x,y
        // such that u^2+v^2 = aperture^2
        double xAp = aperture / det * 
            sqrt(localD(0,0)*localD(0,0) + localD(0,1)*localD(0,1));
        double yAp = aperture / det * 
            sqrt(localD(1,0)*localD(1,0) + localD(1,1)*localD(1,1));
        xdbg<<"aperture = "<<aperture<<std::endl;
        xdbg<<"xap = "<<xAp<<", yap = "<<yAp<<std::endl;

        int xMin = im.getXMin();
        int yMin = im.getYMin();

        double xCen = cen.getX();
        double yCen = cen.getY();
        xdbg<<"cen = "<<xCen<<"  "<<yCen<<std::endl;
        xdbg<<"xmin, ymin = "<<xMin<<"  "<<yMin<<std::endl;
        // xCen,yCen are given on a 1-based grid.
        // ie. where the lower left corner pixel is (1,1), rather than (0,0).
        // The easiest way to do this is to just decrease xCen,yCen by 1 each:
        //--xCen; --yCen;
        // === This is now handled by x_offset, y_offset in ReadCatalog

        int i1 = int(floor(xCen-xAp-xMin));
        int i2 = int(ceil(xCen+xAp-xMin));
        int j1 = int(floor(yCen-yAp-yMin));
        int j2 = int(ceil(yCen+yAp-yMin));
        xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

        if (i1 < 0) { i1 = 0; flag |= EDGE; }
        if (i2 > int(im.getMaxI())) { i2 = im.getMaxI(); flag |= EDGE; }
        if (j1 < 0) { j1 = 0; flag |= EDGE; }
        if (j2 > int(im.getMaxJ())) { j2 = im.getMaxJ(); flag |= EDGE; }
        xdbg<<"i1,i2,j1,j2 => "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

        double apsq = aperture*aperture;

        // Do this next loop in two passes.  First figure out which 
        // pixels we want to use.  Then we can resize pix to the full size
        // we will need, and go back through and enter the pixels.
        // This saves us a lot of resizing calls in vector, which are
        // both slow and can fragment the memory.
        xdbg<<"nx = "<<i2-i1+1<<std::endl;
        xdbg<<"ny = "<<j2-j1+1<<std::endl;
        Assert(i2-i1+1 >= 0);
        Assert(j2-j1+1 >= 0);
        std::vector<std::vector<bool> > shouldUsePix(
            i2-i1+1,std::vector<bool>(j2-j1+1,false));
        int nPix = 0;

        double chipX = xMin+i1-xCen;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = yMin+j1-yCen;
            double u = localD(0,0)*chipX+localD(0,1)*chipY;
            double v = localD(1,0)*chipX+localD(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=localD(0,1),v+=localD(1,1)) {
                // u,v are in arcsec
                double rsq = u*u + v*v;
                if (rsq <= apsq) {
                    shouldUsePix[i-i1][j-j1] = true;
                    ++nPix;
                }
            }
        }

        xdbg<<"npix = "<<nPix<<std::endl;
        pix.resize(nPix);

        xdbg<<"pixlist size = "<<nPix<<" = "<<nPix*sizeof(Pixel)<<" bytes = "<<nPix*sizeof(Pixel)/1024.<<" KB\n";

        int k=0;
        chipX = xMin+i1-xCen;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = yMin+j1-yCen;
            double u = localD(0,0)*chipX+localD(0,1)*chipY;
            double v = localD(1,0)*chipX+localD(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=localD(0,1),v+=localD(1,1)) {
                if (shouldUsePix[i-i1][j-j1]) {
                    double flux = im(i,j)-sky;
                    double inverseVariance;
                    if (weightImage) {
                        inverseVariance = (*weightImage)(i,j);
                    } else {
                        double var = noise;
                        if (gain != 0.) var += im(i,j)/gain;
                        inverseVariance = 1./var;
                    }
                    if (inverseVariance > 0.0) {
                        double inverseSigma = sqrt(inverseVariance);
                        Assert(k < int(pix.size()));
                        pix[k++] = Pixel(u,v,flux,inverseSigma);
                    }
                }
            }
        }
        Assert(k <= int(pix.size()));
        // Not necessarily == because we skip pixels with 0.0 variance
        pix.resize(k);
        Assert(k == int(pix.size()));
        nPix = pix.size(); // may have changed.
        xdbg<<"npix => "<<nPix<<std::endl;
        if (nPix < 10) flag |= LT10PIX;
    }

    double getLocalSky(
        const Image<float>& bkg, 
        const Position cen, const Transformation& trans,
        double aperture, long& flag)
    {
        // This function is very similar in structure to the above getPixList
        // function.  It does the same thing with the distortion and the 
        // aperture and such.  
        // The return value is the mean sky value within the aperture.

        xdbg<<"Start GetLocalSky\n";

        DSmallMatrix22 localD;
        trans.getDistortion(cen,localD);

        double det = std::abs(localD.TMV_det());
        double pixScale = sqrt(det); // arcsec/pixel
        xdbg<<"pixscale = "<<pixScale<<std::endl;

        // xAp,yAp are the maximum deviation from the center in x,y
        // such that u^2+v^2 = aperture^2
        double xAp = aperture / det * 
            sqrt(localD(0,0)*localD(0,0) + localD(0,1)*localD(0,1));
        double yAp = aperture / det * 
            sqrt(localD(1,0)*localD(1,0) + localD(1,1)*localD(1,1));
        xdbg<<"aperture = "<<aperture<<std::endl;
        xdbg<<"xap = "<<xAp<<", yap = "<<yAp<<std::endl;

        int xMin = bkg.getXMin();
        int yMin = bkg.getYMin();

        double xCen = cen.getX();
        double yCen = cen.getY();

        int i1 = int(floor(xCen-xAp-xMin));
        int i2 = int(ceil(xCen+xAp-xMin));
        int j1 = int(floor(yCen-yAp-yMin));
        int j2 = int(ceil(yCen+yAp-yMin));
        xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;
        if (i1 < 0) { i1 = 0; }
        if (i2 > int(bkg.getMaxI())) { i2 = bkg.getMaxI(); }
        if (j1 < 0) { j1 = 0; }
        if (j2 > int(bkg.getMaxJ())) { j2 = bkg.getMaxJ(); }
        xdbg<<"i1,i2,j1,j2 => "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

        double apsq = aperture*aperture;

        xdbg<<"nx = "<<i2-i1+1<<std::endl;
        xdbg<<"ny = "<<j2-j1+1<<std::endl;
        Assert(i2-i1+1 >= 0);
        Assert(j2-j1+1 >= 0);

        double meanSky = 0.;
        int nPix = 0;

        double chipX = xMin+i1-xCen;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = yMin+j1-yCen;
            double u = localD(0,0)*chipX+localD(0,1)*chipY;
            double v = localD(1,0)*chipX+localD(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=localD(0,1),v+=localD(1,1)) {
                // u,v are in arcsec
                double rsq = u*u + v*v;
                if (rsq <= apsq) {
                    meanSky += bkg(i,j);
                    ++nPix;
                }
            }
        }

        xdbg<<"nPix = "<<nPix<<std::endl;
        if (nPix == 0) { flag |= BKG_NOPIX; return 0.; }

        meanSky /= nPix;
        xdbg<<"meansky = "<<meanSky<<std::endl;
        return meanSky;
    }
#endif

    void getSubPixList(PixelList& pix, const PixelList& allpix,
                       double aperture, long& flag)
    {
        // Select a subset of allpix that are within the given aperture
        xdbg<<"Start GetSubPixList\n";
        xdbg<<"allpix has "<<allpix.size()<<" objects\n";
        xdbg<<"new apertur size is "<<aperture<<std::endl;

        double apsq = aperture*aperture;

        pix.clear();
        const int nPix = allpix.size();
        for(int i=0;i<nPix;++i) {
            if (std::norm(allpix[i].getPos()) < apsq) {
                pix.push_back(allpix[i]);
            }
        }
        xdbg<<"done: npix = "<<pix.size()<<std::endl;
        if (pix.size() < 10) flag |= LT10PIX;
    }

}}}}
