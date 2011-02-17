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
        double aperture, double xOffset, double yOffset, long& flag)
    {
        xdbg<<"Start GetPixList\n";
        if (weightImage) {
            xdbg<<"Using weight image for pixel noise.\n";
        } else {
            xdbg<<"noise = "<<noise<<std::endl;
            xdbg<<"gain = "<<gain<<std::endl;
        }

        DSmallMatrix22 D;
        trans.getDistortion(cen,D);
        // This gets us from chip coordinates to sky coordinates:
        // ( u ) = D ( x )
        // ( v )     ( y )
        xdbg<<"D = "<<D<<std::endl;
        double det = std::abs(D.TMV_det());
        double pixScale = sqrt(det); // arsec/pixel
        xdbg<<"pixscale = "<<pixScale<<std::endl;

        // xAp,yAp are the maximum deviation from the center in x,y
        // such that u'^2+v'^2 = aperture^2
        double xAp = aperture * sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1))/det;
        double yAp = aperture * sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1))/det;
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
        // This is now handled by parameters xOffset and yOffset, which are
        // set from the parameters: cat_x_offset and cat_y_offset
        xCen -= xOffset;
        yCen -= yOffset;

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
        double peak = 0.;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = yMin+j1-yCen;
            double u = D(0,0)*chipX+D(0,1)*chipY;
            double v = D(1,0)*chipX+D(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
                double rsq = u*u+v*v;
                if (rsq <= apsq) {
                    shouldUsePix[i-i1][j-j1] = true;
                    ++nPix;
                    if (im(i,j) > peak) peak = im(i,j);
                }
            }
        }

        xdbg<<"npix = "<<nPix<<std::endl;
        pix.resize(nPix);

        xdbg<<"pixlist size = "<<nPix<<" = "<<nPix*sizeof(Pixel)<<
            " bytes = "<<nPix*sizeof(Pixel)/1024.<<" KB\n";

        int k=0;
        chipX = xMin+i1-xCen;
        xdbg<<"Bright pixels are:\n";
        peak -= sky;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = yMin+j1-yCen;
            double u = D(0,0)*chipX+D(0,1)*chipY;
            double v = D(1,0)*chipX+D(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
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
                        Pixel p(u,v,flux,inverseSigma);
                        pix[k++] = p;
                        if (flux > peak / 10.) {
                            xdbg<<p.getPos()<<"  "<<p.getFlux()<<std::endl;
                        }
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
        const Image<double>& bkg, 
        const Position cen, const Transformation& trans, double aperture,
        double xOffset, double yOffset, long& flag)
    {
        // This function is very similar in structure to the above getPixList
        // function.  It does the same thing with the distortion and the 
        // aperture and such.  
        // The return value is the mean sky value within the aperture.

        xdbg<<"Start GetLocalSky\n";

        DSmallMatrix22 D;
        trans.getDistortion(cen,D);

        double det = std::abs(D.TMV_det());
        double pixScale = sqrt(det); // arcsec/pixel
        xdbg<<"pixscale = "<<pixScale<<std::endl;

        // xAp,yAp are the maximum deviation from the center in x,y
        // such that u^2+v^2 = aperture^2
        double xAp = aperture / det * 
            sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1));
        double yAp = aperture / det * 
            sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1));
        xdbg<<"aperture = "<<aperture<<std::endl;
        xdbg<<"xap = "<<xAp<<", yap = "<<yAp<<std::endl;

        int xMin = bkg.getXMin();
        int yMin = bkg.getYMin();

        double xCen = cen.getX();
        double yCen = cen.getY();
        xdbg<<"cen = "<<xCen<<"  "<<yCen<<std::endl;
        xdbg<<"xmin, ymin = "<<xMin<<"  "<<yMin<<std::endl;
        xCen -= xOffset;
        yCen -= yOffset;

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
            double u = D(0,0)*chipX+D(0,1)*chipY;
            double v = D(1,0)*chipX+D(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
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

    void getSubPixList(
        PixelList& pix, const PixelList& allPix,
        std::complex<double> cen_offset, std::complex<double> shear,
        double aperture, long& flag)
    {
        // Select a subset of allPix that are within the given aperture
        const int nTot = allPix.size();
        xdbg<<"Start GetSubPixList\n";
        xdbg<<"allPix has "<<nTot<<" objects\n";
        xdbg<<"new aperture = "<<aperture<<std::endl;
        xdbg<<"cen_offset = "<<cen_offset<<std::endl;
        xdbg<<"shear = "<<shear<<std::endl;

        double normg = norm(shear);
        double g1 = real(shear);
        double g2 = imag(shear);
        double apsq = aperture*aperture;

        // Do this next loop in two passes.  First figure out which 
        // pixels we want to use.  Then we can resize pix to the full size
        // we will need, and go back through and enter the pixels.
        // This saves us a lot of resizing calls in vector, which are
        // both slow and can fragment the memory.
        std::vector<bool> shouldUsePix(nTot,false);
        int nPix = 0;

        double peak = 0.;
        for(int i=0;i<nTot;++i) {
            std::complex<double> z = allPix[i].getPos() - cen_offset;
            double u = real(z);
            double v = imag(z);
            // (1 + |g|^2) (u^2+v^2) - 2g1 (u^2-v^2) - 2g2 (2uv)
            // u,v are in arcsec
            double usq = u*u;
            double vsq = v*v;
            double rsq = (1.+normg)*(usq+vsq) - 2.*g1*(usq-vsq) - 4.*g2*u*v;
            rsq /= (1.-normg);
            if (rsq <= apsq) {
                //xdbg<<"u,v = "<<u<<','<<v<<"  rsq = "<<rsq<<std::endl;
                shouldUsePix[i] = true;
                ++nPix;
                if (allPix[i].getFlux() > peak) peak = allPix[i].getFlux();
            }
        }

        xdbg<<"npix = "<<nPix<<std::endl;
        pix.resize(nPix);

        xdbg<<"pixlist size = "<<nPix<<" = "<<nPix*sizeof(Pixel)<<
            " bytes = "<<nPix*sizeof(Pixel)/1024.<<" KB\n";

        int k=0;
        xdbg<<"Bright pixels are:\n";
        for(int i=0;i<nTot;++i) if(shouldUsePix[i]) {
            Pixel p = allPix[i];
            p.setPos(p.getPos() - cen_offset);
            pix[k++] = p;
            if (p.getFlux() > peak / 10.) {
                xdbg<<p.getPos()<<"  "<<p.getFlux()<<std::endl;
            }
        }
        Assert(k == int(pix.size()));

        if (nPix < 10) flag |= LT10PIX;
    }

}}}}
