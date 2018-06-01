// -*- lsst-c++ -*-

/**
 *  @file lsst/meas/algorithms/KernelPsfFactory.h
 *
 *  Utilities for persisting KernelPsf and subclasses thereof.  Should only be included
 *  directly in source files and never swigged.
 *
 *  Implementations are in KernelPsf.cc.
 */

#include "lsst/meas/algorithms/KernelPsf.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/aggregates.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  @brief A read-only singleton struct containing the schema and key used in persistence for KernelPsf.
 */
struct KernelPsfPersistenceHelper {
    afw::table::Schema schema;
    afw::table::Key<int> kernel;
    afw::table::PointKey<double> averagePosition;

    static KernelPsfPersistenceHelper const& get();

    // No copying
    KernelPsfPersistenceHelper(const KernelPsfPersistenceHelper&) = delete;
    KernelPsfPersistenceHelper& operator=(const KernelPsfPersistenceHelper&) = delete;

    // No moving
    KernelPsfPersistenceHelper(KernelPsfPersistenceHelper&&) = delete;
    KernelPsfPersistenceHelper& operator=(KernelPsfPersistenceHelper&&) = delete;

private:
    KernelPsfPersistenceHelper();
};

/**
 *  @brief A PersistableFactory for KernelPsf and its subclasses.
 *
 *  If a KernelPsf subclass has no data members other than its kernel, table persistence for
 *  it can be implemented simply by reimplementing getPersistenceName() and registering
 *  a specialization of KernelPsfFactory.
 *
 *  @tparam T    KernelPsf subclass the factory will construct.
 *  @tparam K    Kernel subclass the Psf constructor requires.
 */
template <typename T = KernelPsf, typename K = afw::math::Kernel>
class KernelPsfFactory : public afw::table::io::PersistableFactory {
public:
    virtual PTR(afw::table::io::Persistable) read(afw::table::io::InputArchive const& archive,
                                                  afw::table::io::CatalogVector const& catalogs) const {
        static KernelPsfPersistenceHelper const& keys = KernelPsfPersistenceHelper::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        afw::table::BaseRecord const& record = catalogs.front().front();
        LSST_ARCHIVE_ASSERT(record.getSchema() == keys.schema);
        return PTR(T)(new T(archive.get<K>(record.get(keys.kernel)), record.get(keys.averagePosition)));
    }

    KernelPsfFactory(std::string const& name) : afw::table::io::PersistableFactory(name) {}
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
