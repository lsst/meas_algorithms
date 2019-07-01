// -*- LSST-C++ -*-

#include "lsst/geom/Box.h"
#include "lsst/afw/table/io/Persistable.cc"
#include "lsst/meas/algorithms/KernelPsf.h"
#include "lsst/meas/algorithms/KernelPsfFactory.h"

namespace lsst {
namespace afw {
namespace table {
namespace io {

template std::shared_ptr<meas::algorithms::KernelPsf>
PersistableFacade<meas::algorithms::KernelPsf>::dynamicCast(std::shared_ptr<Persistable> const&);

}  // namespace io
}  // namespace table
}  // namespace afw
namespace meas {
namespace algorithms {

PTR(afw::detection::Psf::Image)
KernelPsf::doComputeKernelImage(geom::Point2D const& position, afw::image::Color const& color) const {
    PTR(Psf::Image) im = std::make_shared<Psf::Image>(_kernel->getDimensions());
    _kernel->computeImage(*im, true, position.getX(), position.getY());
    return im;
}

geom::Box2I KernelPsf::doComputeBBox(geom::Point2D const& position, afw::image::Color const& color) const {
    return _kernel->getBBox();
}

KernelPsf::KernelPsf(afw::math::Kernel const& kernel, geom::Point2D const& averagePosition)
        : ImagePsf(!kernel.isSpatiallyVarying()),
          _kernel(kernel.clone()),
          _averagePosition(averagePosition) {}

KernelPsf::KernelPsf(PTR(afw::math::Kernel) kernel, geom::Point2D const& averagePosition)
        : ImagePsf(!kernel->isSpatiallyVarying()), _kernel(kernel), _averagePosition(averagePosition) {}

PTR(afw::detection::Psf) KernelPsf::clone() const { return std::make_shared<KernelPsf>(*this); }

PTR(afw::detection::Psf) KernelPsf::resized(int width, int height) const {
    return std::make_shared<KernelPsf>(*_kernel->resized(width, height), _averagePosition);
}

geom::Point2D KernelPsf::getAveragePosition() const { return _averagePosition; }

namespace {

KernelPsfFactory<> registration("KernelPsf");

}  // namespace

KernelPsfPersistenceHelper const& KernelPsfPersistenceHelper::get() {
    static KernelPsfPersistenceHelper instance;
    return instance;
}

KernelPsfPersistenceHelper::KernelPsfPersistenceHelper()
        : schema(),
          kernel(schema.addField<int>("kernel", "archive ID of nested kernel object")),
          averagePosition(afw::table::PointKey<double>::addFields(
                  schema, "averagePosition", "average position of stars used to make the PSF", "pixel")) {}

bool KernelPsf::isPersistable() const noexcept { return _kernel->isPersistable(); }

std::string KernelPsf::getPersistenceName() const { return "KernelPsf"; }

std::string KernelPsf::getPythonModule() const { return "lsst.meas.algorithms"; }

void KernelPsf::write(OutputArchiveHandle& handle) const {
    static KernelPsfPersistenceHelper const& keys = KernelPsfPersistenceHelper::get();
    afw::table::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    PTR(afw::table::BaseRecord) record = catalog.addNew();
    record->set(keys.kernel, handle.put(_kernel));
    record->set(keys.averagePosition, _averagePosition);
    handle.saveCatalog(catalog);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
