// -*- lsst-c++ -*-

#include "boost/make_shared.hpp"

#include "lsst/meas/algorithms/SkyCoord.h"

namespace lsst { namespace meas { namespace algorithms {

namespace {

class SkyCoordAlgorithm : public Algorithm {
public:

    explicit SkyCoordAlgorithm(SkyCoordControl const & ctrl);

public:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

private:
    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SkyCoordAlgorithm);
};

SkyCoordAlgorithm::SkyCoordAlgorithm(SkyCoordControl const & ctrl) : Algorithm(ctrl) {}

template <typename PixelT>
void SkyCoordAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    if (!exposure.hasWcs()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicError,
            "Cannot run SkyCoord algorithm with no WCS."
        );
    }
    source.updateCoord(*exposure.getWcs());
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(SkyCoordAlgorithm);

} // anonymous

PTR(AlgorithmControl) SkyCoordControl::_clone() const {
    return boost::make_shared<SkyCoordControl>(*this);
}

PTR(Algorithm) SkyCoordControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<SkyCoordAlgorithm>(*this);
}

}}} // namespace lsst::meas::algorithms
