/*
 * LSST Data Management System
 * Copyright 2018 LSST/AURA.
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

#include "lsst/meas/algorithms/CoaddTransmissionCurve.h"
#include "lsst/geom/Box.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

struct CoaddInputData {
    std::shared_ptr<afw::image::TransmissionCurve const> transmission;
    std::shared_ptr<afw::geom::SkyWcs const> sensorWcs;
    std::shared_ptr<afw::geom::polygon::Polygon const> validPolygon;
    geom::Box2D bbox;
    double weight;
};

class CoaddTransmissionCurve : public afw::image::TransmissionCurve {
public:
    CoaddTransmissionCurve(std::shared_ptr<afw::geom::SkyWcs const> coaddWcs,
                           afw::table::ExposureCatalog const& inputSensors)
            : _coaddWcs(coaddWcs),
              _wavelengthBounds(std::numeric_limits<double>::infinity(),  // min = +inf, max = -inf *for now*
                                -std::numeric_limits<double>::infinity()),
              _throughputAtBounds(0.0, 0.0) {
        _inputs.reserve(inputSensors.size());
        afw::table::Key<double> weightKey = inputSensors.getSchema()["weight"];
        double weightSum = 0.0;
        for (auto const& record : inputSensors) {
            double const weight = record.get(weightKey);
            weightSum += weight;
            auto transmission = record.getTransmissionCurve();
            if (transmission == nullptr) {
                throw LSST_EXCEPT(
                        pex::exceptions::InvalidParameterError,
                        (boost::format("No TransmissionCurve for input with ID %d") % record.getId()).str());
            }
            if (record.getWcs() == nullptr) {
                throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                                  (boost::format("No Wcs for input with ID %d") % record.getId()).str());
            }
            // Set wavelength bounds from the overall minimum and maximum wavelengths from all inputs.
            auto const inputWavelengthBounds = transmission->getWavelengthBounds();
            _wavelengthBounds.first = std::min(_wavelengthBounds.first, inputWavelengthBounds.first);
            _wavelengthBounds.second = std::max(_wavelengthBounds.second, inputWavelengthBounds.second);
            // Set throughput at and outside bounds from weighted average of throughput at and outside bounds.
            auto const inputThroughputAtBounds = transmission->getThroughputAtBounds();
            _throughputAtBounds.first += inputThroughputAtBounds.first * weight;
            _throughputAtBounds.second += inputThroughputAtBounds.second * weight;
            // Add an element to the vector with all the stuff we need to store for each epoch.
            CoaddInputData input = {transmission, record.getWcs(), record.getValidPolygon(),
                                    geom::Box2D(record.getBBox()), weight};
            _inputs.push_back(std::move(input));
        }
        _throughputAtBounds.first /= weightSum;
        _throughputAtBounds.second /= weightSum;
    }

    // Private constructor used only for persistence
    CoaddTransmissionCurve(std::shared_ptr<afw::geom::SkyWcs> coaddWcs,
                           std::pair<double, double> wavelengthBounds,
                           std::pair<double, double> throughputAtBounds, std::vector<CoaddInputData> inputs)
            : _coaddWcs(std::move(coaddWcs)),
              _wavelengthBounds(std::move(wavelengthBounds)),
              _throughputAtBounds(std::move(throughputAtBounds)),
              _inputs(std::move(inputs)) {}

    // All TransmissionCurves are immutable and noncopyable.
    CoaddTransmissionCurve(CoaddTransmissionCurve const&) = delete;
    CoaddTransmissionCurve(CoaddTransmissionCurve&&) = delete;
    CoaddTransmissionCurve& operator=(CoaddTransmissionCurve const&) = delete;
    CoaddTransmissionCurve& operator=(CoaddTransmissionCurve&&) = delete;

    std::pair<double, double> getWavelengthBounds() const override { return _wavelengthBounds; }

    std::pair<double, double> getThroughputAtBounds() const override { return _throughputAtBounds; }

    void sampleAt(geom::Point2D const& position, ndarray::Array<double const, 1, 1> const& wavelengths,
                  ndarray::Array<double, 1, 1> const& out) const override {
        auto const coord = _coaddWcs->pixelToSky(position);
        ndarray::Array<double, 1, 1> tmp = ndarray::allocate(out.getShape());
        tmp.deep() = 0.0;
        out.deep() = 0.0;
        double weightSum = 0.0;
        for (auto const& input : _inputs) {
            geom::Point2D const inputPosition = input.sensorWcs->skyToPixel(coord);
            if (!input.bbox.contains(inputPosition)) {
                continue;
            }
            if (input.validPolygon && !input.validPolygon->contains(inputPosition)) {
                continue;
            }
            // note that `tmp` is an output argument here
            input.transmission->sampleAt(inputPosition, wavelengths, tmp);
            tmp.deep() *= input.weight;
            out.deep() += tmp;
            weightSum += input.weight;
        }
        if (weightSum == 0.0) {
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                              (boost::format("No input TransmissionCurves at point (%s, %s)") %
                               position.getX() % position.getY())
                                      .str());
        }
        out.deep() /= weightSum;
    }

    bool isPersistable() const noexcept override {
        for (auto const& input : _inputs) {
            if (!input.transmission->isPersistable()) {
                return false;
            }
        }
        return true;
    }

private:
    std::string getPersistenceName() const override { return "CoaddTransmissionCurve"; }

    std::string getPythonModule() const override { return "lsst.meas.algorithms"; }

    struct PersistenceHelper;

    class Factory;

    void write(OutputArchiveHandle& handle) const override;

    static Factory registration;

    std::shared_ptr<afw::geom::SkyWcs const> _coaddWcs;
    std::pair<double, double> _wavelengthBounds;
    std::pair<double, double> _throughputAtBounds;
    std::vector<CoaddInputData> _inputs;
};

struct CoaddTransmissionCurve::PersistenceHelper {
    afw::table::Schema mainSchema;
    afw::table::Key<int> coaddWcs;
    afw::table::Key<double> wavelengthMin;
    afw::table::Key<double> wavelengthMax;
    afw::table::Key<double> throughputAtMin;
    afw::table::Key<double> throughputAtMax;
    afw::table::Schema inputDataSchema;
    afw::table::Key<int> transmission;
    afw::table::Key<int> sensorWcs;
    afw::table::Key<int> validPolygon;
    afw::table::Box2DKey bbox;
    afw::table::Key<double> weight;

    static PersistenceHelper const& get() {
        static PersistenceHelper const instance;
        return instance;
    }

private:
    PersistenceHelper()
            : mainSchema(),
              coaddWcs(mainSchema.addField<int>("coaddWcs", "archive ID for the coadd's WCS")),
              wavelengthMin(mainSchema.addField<double>("wavelengthMin", "getWavelengthBounds().min")),
              wavelengthMax(mainSchema.addField<double>("wavelengthMax", "getWavelengthBounds().max")),
              throughputAtMin(mainSchema.addField<double>("throughputAtMin", "throughputAtBounds().first")),
              throughputAtMax(mainSchema.addField<double>("throughputAtMax", "throughputAtBounds().second")),
              inputDataSchema(),
              transmission(inputDataSchema.addField<int>("transmission",
                                                         "archive ID for this sensor's TransmissionCurve")),
              sensorWcs(inputDataSchema.addField<int>("sensorWcs", "archive ID for this sensor's WCS")),
              validPolygon(inputDataSchema.addField<int>("validPolygon",
                                                         "archive ID for this sensor's ValidPolygon")),
              bbox(afw::table::Box2DKey::addFields(inputDataSchema, "bbox", "bounding box of the sensor",
                                                   "pixel")),
              weight(inputDataSchema.addField<double>("weight",
                                                      "relative weight for this sensor in the average")) {}
};

void CoaddTransmissionCurve::write(OutputArchiveHandle& handle) const {
    auto const& keys = PersistenceHelper::get();
    auto mainCat = handle.makeCatalog(keys.mainSchema);
    auto mainRecord = mainCat.addNew();
    mainRecord->set(keys.coaddWcs, handle.put(_coaddWcs));
    mainRecord->set(keys.wavelengthMin, _wavelengthBounds.first);
    mainRecord->set(keys.wavelengthMax, _wavelengthBounds.second);
    mainRecord->set(keys.throughputAtMin, _throughputAtBounds.first);
    mainRecord->set(keys.throughputAtMax, _throughputAtBounds.second);
    handle.saveCatalog(mainCat);
    auto inputDataCat = handle.makeCatalog(keys.inputDataSchema);
    for (auto const& input : _inputs) {
        auto inputDataRecord = inputDataCat.addNew();
        inputDataRecord->set(keys.transmission, handle.put(input.transmission));
        inputDataRecord->set(keys.sensorWcs, handle.put(input.sensorWcs));
        inputDataRecord->set(keys.validPolygon, handle.put(input.validPolygon));
        inputDataRecord->set(keys.bbox, input.bbox);
        inputDataRecord->set(keys.weight, input.weight);
    }
    handle.saveCatalog(inputDataCat);
}

class CoaddTransmissionCurve::Factory : public afw::table::io::PersistableFactory {
public:
    std::shared_ptr<afw::table::io::Persistable> read(InputArchive const& archive,
                                                      CatalogVector const& catalogs) const override {
        auto const& keys = PersistenceHelper::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 2u);
        LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys.mainSchema);
        LSST_ARCHIVE_ASSERT(catalogs.back().getSchema() == keys.inputDataSchema);
        auto const& mainRecord = catalogs.front().front();
        std::vector<CoaddInputData> inputs;
        inputs.reserve(catalogs.back().size());
        for (auto const& inputDataRecord : catalogs.back()) {
            CoaddInputData input = {
                    archive.get<afw::image::TransmissionCurve>(inputDataRecord.get(keys.transmission)),
                    archive.get<afw::geom::SkyWcs>(inputDataRecord.get(keys.sensorWcs)),
                    archive.get<afw::geom::polygon::Polygon>(inputDataRecord.get(keys.validPolygon)),
                    inputDataRecord.get(keys.bbox), inputDataRecord.get(keys.weight)};
            inputs.push_back(std::move(input));
        }
        return std::make_shared<CoaddTransmissionCurve>(
                archive.get<afw::geom::SkyWcs>(mainRecord.get(keys.coaddWcs)),
                std::make_pair(mainRecord.get(keys.wavelengthMin), mainRecord.get(keys.wavelengthMax)),
                std::make_pair(mainRecord.get(keys.throughputAtMin), mainRecord.get(keys.throughputAtMax)),
                std::move(inputs));
    }

    Factory(std::string const& name) : afw::table::io::PersistableFactory(name) {}
};

CoaddTransmissionCurve::Factory CoaddTransmissionCurve::registration("CoaddTransmissionCurve");

}  // namespace

std::shared_ptr<afw::image::TransmissionCurve const> makeCoaddTransmissionCurve(
        std::shared_ptr<afw::geom::SkyWcs const> coaddWcs, afw::table::ExposureCatalog const& inputSensors) {
    return std::make_shared<CoaddTransmissionCurve>(coaddWcs, inputSensors);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
