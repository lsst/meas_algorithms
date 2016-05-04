// -*- LSST-C++ -*-

#include "lsst/meas/algorithms/KernelPsf.h"
#include "lsst/meas/algorithms/KernelPsfFactory.h"

namespace lsst { namespace meas { namespace algorithms {

PTR(afw::detection::Psf::Image) KernelPsf::doComputeKernelImage(
    afw::geom::Point2D const & position, afw::image::Color const& color
) const {
    PTR(Psf::Image) im = std::make_shared<Psf::Image>(_kernel->getDimensions());
    _kernel->computeImage(*im, true, position.getX(), position.getY());
    return im;
}

KernelPsf::KernelPsf(afw::math::Kernel const & kernel, afw::geom::Point2D const & averagePosition) :
    ImagePsf(!kernel.isSpatiallyVarying()), _kernel(kernel.clone()), _averagePosition(averagePosition) {}

KernelPsf::KernelPsf(PTR(afw::math::Kernel) kernel, afw::geom::Point2D const & averagePosition) :
    ImagePsf(!kernel->isSpatiallyVarying()), _kernel(kernel), _averagePosition(averagePosition) {}

PTR(afw::detection::Psf) KernelPsf::clone() const { return std::make_shared<KernelPsf>(*this); }

afw::geom::Point2D KernelPsf::getAveragePosition() const { return _averagePosition; }

namespace {

KernelPsfFactory<> registration("KernelPsf");

} // anonymous

KernelPsfPersistenceHelper const & KernelPsfPersistenceHelper::get() {
    static KernelPsfPersistenceHelper instance;
    return instance;
}

KernelPsfPersistenceHelper::KernelPsfPersistenceHelper() :
    schema(),
    kernel(schema.addField<int>("kernel", "archive ID of nested kernel object")),
    averagePosition(afw::table::PointKey<double>::addFields(
                        schema, "averagePosition", "average position of stars used to make the PSF", "pixels"
                    ))
{
    schema.getCitizen().markPersistent();
}

bool KernelPsf::isPersistable() const { return _kernel->isPersistable(); }

std::string KernelPsf::getPersistenceName() const { return "KernelPsf"; }

std::string KernelPsf::getPythonModule() const { return "lsst.meas.algorithms"; }

void KernelPsf::write(OutputArchiveHandle & handle) const {
    static KernelPsfPersistenceHelper const & keys = KernelPsfPersistenceHelper::get();
    afw::table::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    PTR(afw::table::BaseRecord) record = catalog.addNew();
    record->set(keys.kernel, handle.put(_kernel));
    record->set(keys.averagePosition, _averagePosition);
    handle.saveCatalog(catalog);
}

}}} // namespace lsst::meas::algorithms
