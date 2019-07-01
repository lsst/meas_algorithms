// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include <cmath>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/Persistable.cc"
#include "lsst/afw/table/aggregates.h"
#include "lsst/meas/algorithms/DoubleGaussianPsf.h"

namespace lsst {
namespace afw {
namespace table {
namespace io {

template std::shared_ptr<meas::algorithms::DoubleGaussianPsf>
PersistableFacade<meas::algorithms::DoubleGaussianPsf>::dynamicCast(std::shared_ptr<Persistable> const&);

}  // namespace io
}  // namespace table
}  // namespace afw
namespace meas {
namespace algorithms {

namespace {

// Read-only singleton struct containing the schema and keys that a double-Gaussian Psf is mapped
// to in record persistence.
struct DoubleGaussianPsfPersistenceHelper {
    afw::table::Schema schema;
    afw::table::PointKey<int> dimensions;
    afw::table::Key<double> sigma1;
    afw::table::Key<double> sigma2;
    afw::table::Key<double> b;

    static DoubleGaussianPsfPersistenceHelper const& get() {
        static DoubleGaussianPsfPersistenceHelper instance;
        return instance;
    }

    // No copying
    DoubleGaussianPsfPersistenceHelper(const DoubleGaussianPsfPersistenceHelper&) = delete;
    DoubleGaussianPsfPersistenceHelper& operator=(const DoubleGaussianPsfPersistenceHelper&) = delete;

    // No moving
    DoubleGaussianPsfPersistenceHelper(DoubleGaussianPsfPersistenceHelper&&) = delete;
    DoubleGaussianPsfPersistenceHelper& operator=(DoubleGaussianPsfPersistenceHelper&&) = delete;

private:
    DoubleGaussianPsfPersistenceHelper()
            : schema(),
              dimensions(afw::table::PointKey<int>::addFields(schema, "dimensions", "width/height of kernel",
                                                              "pixel")),
              sigma1(schema.addField<double>("sigma1", "radius of inner Gaussian", "pixel")),
              sigma2(schema.addField<double>("sigma2", "radius of outer Gaussian", "pixel")),
              b(schema.addField<double>("b", "central amplitude of outer Gaussian (inner amplitude == 1)")) {}
};

class DoubleGaussianPsfFactory : public afw::table::io::PersistableFactory {
public:
    virtual PTR(afw::table::io::Persistable)
            read(InputArchive const& archive, CatalogVector const& catalogs) const {
        static DoubleGaussianPsfPersistenceHelper const& keys = DoubleGaussianPsfPersistenceHelper::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        afw::table::BaseRecord const& record = catalogs.front().front();
        LSST_ARCHIVE_ASSERT(record.getSchema() == keys.schema);
        return std::make_shared<DoubleGaussianPsf>(
                record.get(keys.dimensions.getX()), record.get(keys.dimensions.getY()),
                record.get(keys.sigma1), record.get(keys.sigma2), record.get(keys.b));
    }

    DoubleGaussianPsfFactory(std::string const& name) : afw::table::io::PersistableFactory(name) {}
};

// Helper function for ctor: need to construct the kernel to pass to KernelPsf, because we
// can't change it after construction.
PTR(afw::math::Kernel)
makeDoubleGaussianKernel(int width, int height, double sigma1, double& sigma2, double b) {
    if (b == 0.0 && sigma2 == 0.0) {
        sigma2 = 1.0;  // avoid 0/0 at centre of Psf
    }
    if (sigma1 <= 0 || sigma2 <= 0) {
        throw LSST_EXCEPT(pex::exceptions::DomainError,
                          (boost::format("sigma may not be 0: %g, %g") % sigma1 % sigma2).str());
    }
    afw::math::DoubleGaussianFunction2<double> dg(sigma1, sigma2, b);
    PTR(afw::math::Kernel) kernel(new afw::math::AnalyticKernel(width, height, dg));
    return kernel;
}

std::string getDoubleGaussianPsfPersistenceName() { return "DoubleGaussianPsf"; }

DoubleGaussianPsfFactory registration(getDoubleGaussianPsfPersistenceName());

}  // namespace

DoubleGaussianPsf::DoubleGaussianPsf(int width, int height, double sigma1, double sigma2, double b)
        : KernelPsf(makeDoubleGaussianKernel(width, height, sigma1, sigma2, b)),
          _sigma1(sigma1),
          _sigma2(sigma2),
          _b(b) {}

PTR(afw::detection::Psf) DoubleGaussianPsf::clone() const {
    return std::make_shared<DoubleGaussianPsf>(getKernel()->getWidth(), getKernel()->getHeight(), _sigma1,
                                               _sigma2, _b);
}

PTR(afw::detection::Psf) DoubleGaussianPsf::resized(int width, int height) const {
    return std::make_shared<DoubleGaussianPsf>(width, height, _sigma1, _sigma2, _b);
}

std::string DoubleGaussianPsf::getPersistenceName() const { return getDoubleGaussianPsfPersistenceName(); }

void DoubleGaussianPsf::write(OutputArchiveHandle& handle) const {
    static DoubleGaussianPsfPersistenceHelper const& keys = DoubleGaussianPsfPersistenceHelper::get();
    afw::table::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    PTR(afw::table::BaseRecord) record = catalog.addNew();
    (*record).set(keys.dimensions.getX(), getKernel()->getWidth());
    (*record).set(keys.dimensions.getY(), getKernel()->getHeight());
    (*record).set(keys.sigma1, getSigma1());
    (*record).set(keys.sigma2, getSigma2());
    (*record).set(keys.b, getB());
    handle.saveCatalog(catalog);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

namespace lsst {
namespace afw {
namespace detection {

}  // namespace detection
}  // namespace afw
}  // namespace lsst
