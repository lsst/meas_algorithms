
#include "lsst/meas/algorithms/shapelet/Pixel.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

#ifdef PIXELLIST_USE_POOL
    PixelList::PixelList() :
        _shouldUsePool(false), _v1(new std::vector<Pixel>()) 
    {}

    PixelList::PixelList(const int n) :
        _shouldUsePool(false), _v1(new std::vector<Pixel>(n)) 
    {}

    PixelList::PixelList(const PixelList& rhs) :
        _shouldUsePool(false), _v1(new std::vector<Pixel>(rhs.size())) 
    {
        if (rhs._shouldUsePool) 
            std::copy(rhs._v2->begin(),rhs._v2->end(),_v1->begin());
        else *_v1 = *rhs._v1;
    }

    PixelList& PixelList::operator=(const PixelList& rhs)
    {
        if (size() != rhs.size()) resize(rhs.size());

        if (_shouldUsePool) {
            if (rhs._shouldUsePool) *_v2 = *rhs._v2;
            else std::copy(rhs._v1->begin(),rhs._v1->end(),_v2->begin());
        } else {
            if (rhs._shouldUsePool) 
                std::copy(rhs._v2->begin(),rhs._v2->end(),_v1->begin());
            else *_v1 = *rhs._v1;
        }
        return *this;
    }

    PixelList::~PixelList()
    {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            _v2.reset();
        }
    }

    void PixelList::usePool() 
    {
        // This should be done before any elements are added.
        if (_v1.get()) Assert(_v1->size() == 0);
        if (_v2.get()) Assert(_v2->size() == 0);
        _v1.reset();
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            _v2.reset(new std::vector<Pixel,PoolAllocPixel>());
        }
        _shouldUsePool = true; 
    }

    size_t PixelList::size() const
    {
        if (_shouldUsePool) return _v2->size();
        else return _v1->size();
    }

    void PixelList::reserve(const int n)
    {
        if (_shouldUsePool) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                _v2->reserve(n);
            }
        } else {
            _v1->reserve(n);
        }
    }

    size_t PixelList::capacity() const
    { return _shouldUsePool ? _v2->capacity() : _v1->capacity(); }

    void PixelList::resize(const int n)
    {
        if (_shouldUsePool) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                _v2->resize(n);
            }
        } else {
            _v1->resize(n);
        }
    }

    void PixelList::clear()
    {
        if (_shouldUsePool) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                _v2->clear();
            }
        } else {
            _v1->clear();
        }
    }

    void PixelList::push_back(const Pixel& p)
    {
        if (_shouldUsePool) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                _v2->push_back(p);
            }
        } else {
            _v1->push_back(p);
        }
    }

    Pixel& PixelList::operator[](const int i)
    {
        if (_shouldUsePool) return (*_v2)[i];
        else return (*_v1)[i];
    }

    const Pixel& PixelList::operator[](const int i) const
    {
        if (_shouldUsePool) return (*_v2)[i];
        else return (*_v1)[i];
    }

#else

    PixelList::PixelList() : _v1(new std::vector<Pixel>()) {}

    PixelList::PixelList(const int n) : _v1(new std::vector<Pixel>(n)) {}

    PixelList::PixelList(const PixelList& rhs) :
        _v1(new std::vector<Pixel>(*rhs._v1)) 
    {}

    PixelList& PixelList::operator=(const PixelList& rhs)
    {
        if (size() != rhs.size()) resize(rhs.size());
        *_v1 = *rhs._v1;
        return *this;
    }

    PixelList::~PixelList()
    { _v1.reset(); }

    size_t PixelList::size() const
    { return _v1->size(); }

    void PixelList::reserve(const int n)
    { _v1->reserve(n); }

    size_t PixelList::capacity() const
    { return _v1->capacity(); }

    void PixelList::resize(const int n)
    { _v1->resize(n); }

    void PixelList::clear()
    { _v1->clear(); }

    void PixelList::push_back(const Pixel& p)
    { _v1->push_back(p); }

    Pixel& PixelList::operator[](const int i)
    { return (*_v1)[i]; }

    const Pixel& PixelList::operator[](const int i) const
    { return (*_v1)[i]; }

#endif
}}}}
