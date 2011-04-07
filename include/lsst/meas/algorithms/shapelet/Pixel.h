#ifndef MeasAlgoShapeletPixel_H
#define MeasAlgoShapeletPixel_H

//#define PIXELLIST_USE_POOL

#include <complex>
#include <string>
#if 0
#include "lsst/meas/algorithms/shapelet/Image.h"
#include "lsst/meas/algorithms/shapelet/Transformation.h"
#endif
#include "lsst/meas/algorithms/shapelet/ConfigFile.h"
#include "lsst/meas/algorithms/shapelet/Bounds.h"

#ifdef __INTEL_COMPILER
#pragma warning (disable : 1418)
#endif
#include "boost/shared_ptr.hpp"
#ifdef __INTEL_COMPILER
#pragma warning (default : 1418)
#endif

#ifdef PIXELLIST_USE_POOL
#define PIXELLIST_BLOCK 1024*1024*100  // 100 MB per block
#include "PoolAllocator.h"
#endif

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class Pixel 
    { 
    public :

        Pixel() : _pos(0.), _flux(0.), _inverseSigma(0.) {}

        Pixel(double u, double v, double flux, double inverseSigma) :
            _pos(u,v), _flux(flux), _inverseSigma(inverseSigma) {}

        Pixel(std::complex<double> z, double flux, double inverseSigma) :
            _pos(z), _flux(flux), _inverseSigma(inverseSigma) {}

        ~Pixel() {}

        std::complex<double> getPos() const { return _pos; }

        double getFlux() const { return _flux; }

        double getInverseSigma() const { return _inverseSigma; }

        void setPos(const std::complex<double>& pos) { _pos = pos; }

        void setFlux(const double flux) { _flux = flux; }

        void setInverseSigma(const double inverseSigma) 
        { _inverseSigma = inverseSigma; }

    private :

        std::complex<double> _pos;
        double _flux;
        double _inverseSigma;
    };

    // Most of these methods are not (intrinsically) thread-safe, 
    // since they might be using my pool allocator, so they need to 
    // be wrapped in a critical block.
    // Therefore, all the methods for PixelList are defined in Pixel_omp.cpp.

    class PixelList
    {
    public :

        PixelList();
        PixelList(const int n);
        PixelList(const PixelList& rhs);
        PixelList& operator=(const PixelList& rhs);
        ~PixelList();

        void usePool();

        // These mimic the same functionality of a std::vector<Pixel>
        size_t size() const;
        void reserve(const int n);
        size_t capacity() const;
        void resize(const int n);
        void clear();
        void push_back(const Pixel& p);
        Pixel& operator[](const int i);
        const Pixel& operator[](const int i) const;
        void sort(const Position& cen);

    private :

        bool _shouldUsePool;
        boost::shared_ptr<std::vector<Pixel> > _v1;
#ifdef PIXELLIST_USE_POOL
        typedef PoolAllocator<Pixel,PIXELLIST_BLOCK> PoolAllocPixel;
        boost::shared_ptr<std::vector<Pixel,PoolAllocPixel> > _v2;
#else
        boost::shared_ptr<std::vector<Pixel> > _v2;
#endif

    };

#if 0
    void getPixList(
        const Image<double>& im, PixelList& pix,
        const Position cen, double sky, double noise, double gain,
        const Image<double>* weightImage, const Transformation& trans,
        double aperture, double xOffset, double yOffset, long& flag);

    double getLocalSky(
        const Image<double>& bkg, 
        const Position cen, const Transformation& trans,
        double aperture, double xOffset, double yOffset, long& flag);
#endif

    void getSubPixList(
        PixelList& pix, const PixelList& allPix,
        std::complex<double> cen_offset, std::complex<double> shear,
        double aperture, long& flag);

    inline void getSubPixList(
        PixelList& pix, const PixelList& allPix,
        std::complex<double> cen_offset, double aperture, long& flag)
    { getSubPixList(pix,allPix,cen_offset,0.,aperture,flag); }

    inline void getSubPixList(
        PixelList& pix, const PixelList& allPix,
        double aperture, long& flag)
    { getSubPixList(pix,allPix,0.,0.,aperture,flag); }


}}}}

#endif
