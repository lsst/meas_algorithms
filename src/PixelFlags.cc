// -*- lsst-c++ -*-

#include "boost/make_shared.hpp"

#include "lsst/meas/algorithms/PixelFlags.h"
#include "lsst/afw/detection/FootprintFunctor.h"

namespace lsst { namespace meas { namespace algorithms {

namespace {

template <typename MaskedImageT>
class FootprintBits : public afw::detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintBits(MaskedImageT const& mimage) :
        afw::detection::FootprintFunctor<MaskedImageT>(mimage), _bits(0)
    {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _bits = 0x0;
    }
    virtual void reset(afw::detection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        _bits |= loc.mask(0, 0);
    }

    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    typename MaskedImageT::Mask::Pixel _bits;
};

class PixelFlagAlgorithm : public Algorithm {
public:

    enum {
        EDGE=0,
        BAD,
        INTERPOLATED,
        INTERPOLATED_CENTER,
        SATURATED,
        SATURATED_CENTER,
        CR,
        CR_CENTER,
        N_FLAGS
    };

    typedef boost::array< afw::table::Key<afw::table::Flag>, N_FLAGS > KeyArray;

    PixelFlagAlgorithm(PixelFlagControl const & ctrl, afw::table::Schema & schema);

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(PixelFlagAlgorithm);

    KeyArray _keys;
};

PixelFlagAlgorithm::PixelFlagAlgorithm(PixelFlagControl const & ctrl, afw::table::Schema & schema) :
    Algorithm(ctrl) 
{
    _keys[EDGE] = schema.addField<afw::table::Flag>(
        ctrl.name + ".edge", "source is in region labeled EDGE"
    );
    _keys[INTERPOLATED] = schema.addField<afw::table::Flag>(
        ctrl.name + ".interpolated.any", "source's footprint includes interpolated pixels"
    );
    _keys[INTERPOLATED_CENTER] = schema.addField<afw::table::Flag>(
        ctrl.name + ".interpolated.center", "source's center is close to interpolated pixels"
    );
    _keys[SATURATED] = schema.addField<afw::table::Flag>(
        ctrl.name + ".saturated.any", "source's footprint includes saturated pixels"
    );
    _keys[SATURATED_CENTER] = schema.addField<afw::table::Flag>(
        ctrl.name + ".saturated.center", "source's center is close to saturated pixels"
    );
    _keys[CR] = schema.addField<afw::table::Flag>(
        ctrl.name + ".cr.any", "source's footprint includes suspected CR pixels"
    );
    _keys[CR_CENTER] = schema.addField<afw::table::Flag>(
        ctrl.name + ".cr.center", "source's center is close to suspected CR pixels"
    );
    _keys[BAD] = schema.addField<afw::table::Flag>(
        ctrl.name + ".bad", "source is in region labeled BAD"
    );
}

template <typename PixelT>
void PixelFlagAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;
    FootprintBits<MaskedImageT> func(exposure.getMaskedImage());

    // Catch centroids off the image or NAN
    if (!exposure.getMaskedImage().getBBox().contains(afw::geom::Point2I(center) -
                                                      afw::geom::Extent2I(exposure.getXY0()))) {
        source.set(_keys[EDGE], true);
        return;                         // Can't continue safely
    }

    // Check for bits set in the source's Footprint
    PTR(afw::detection::Footprint) foot = source.getFootprint();
    if (foot) {
        func.apply(*foot);
        if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
            source.set(_keys[EDGE], true);
        }
        if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("BAD")) {
            source.set(_keys[BAD], true);
        }
        if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
            source.set(_keys[INTERPOLATED], true);
        }
        if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
            source.set(_keys[SATURATED], true);
        }
        if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("CR")) {
            source.set(_keys[CR], true);
        }
    }

    // Check for bits set in the 3x3 box around the center
    afw::geom::Point2I llc(afw::image::positionToIndex(center.getX()) - 1,
                           afw::image::positionToIndex(center.getY()) - 1);
    afw::detection::Footprint const middle(afw::geom::BoxI(llc, afw::geom::ExtentI(3))); // central 3x3
    func.apply(middle);
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        source.set(_keys[INTERPOLATED_CENTER], true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        source.set(_keys[SATURATED_CENTER], true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("CR")) {
        source.set(_keys[CR_CENTER], true);
    }

}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(PixelFlagAlgorithm);

} // anonymous

PTR(AlgorithmControl) PixelFlagControl::_clone() const {
    return boost::make_shared<PixelFlagControl>(*this);
}

PTR(Algorithm) PixelFlagControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<PixelFlagAlgorithm>(*this, boost::ref(schema));
}

}}} // namespace lsst::meas::algorithms
