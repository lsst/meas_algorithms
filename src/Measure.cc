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
 
/// \file

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/Measure.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace pexExceptions = lsst::pex::exceptions;    
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;

namespace {
    /*
     * Return the numeric value of name as double
     */
    double getNumeric(lsst::pex::policy::Policy const& policy, std::string const& name)
    {
        return policy.isDouble(name) ? policy.getDouble(name) : policy.getInt(name);
    }
}

/************************************************************************************************************/
/**
 * @brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintCentroid : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintCentroid(MaskedImageT const& mimage ///< The image the source lives in
                              ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                                  _n(0), _sum(0), _sumx(0), _sumy(0),
                                  _min( std::numeric_limits<double>::max()), _xmin(0), _ymin(0),
                                  _max(-std::numeric_limits<double>::max()), _xmax(0), _ymax(0),
                                  _bits(0) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _n = 0;
        _sum = _sumx = _sumy = 0.0;
        _min =  std::numeric_limits<double>::max();
        _xmin = _ymin = 0;
        _max = -std::numeric_limits<double>::max();
        _xmax = _ymax = 0;
        _bits = 0x0;
    }
    virtual void reset(afwDetection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

        _n++;
        _sum += val;
        _sumx += afwImage::indexToPosition(x)*val;
        _sumy += afwImage::indexToPosition(y)*val;
        _bits |= loc.mask(0, 0);

        if (val < _min) {
            _min = val;
            _xmin = x;
            _ymin = y;
        }
        if (val > _max) {
            _max = val;
            _xmax = x;
            _ymax = y;
        }
    }

    /// Return the number of pixels
    int getN() const { return _n; }
    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the Footprint's column centroid
    double getX() const { return _sumx/_sum; }
    /// Return the Footprint's row centroid
    double getY() const { return _sumy/_sum; }
    /// Return the Footprint's peak pixel
    PTR(afwDetection::Peak) makePeak(bool isNegative) const {
        return boost::make_shared<afwDetection::Peak>(isNegative ? afwDetection::Peak(_xmin, _ymin) :
                                                      afwDetection::Peak(_xmax, _ymax));
    }
    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    int _n;
    double _sum, _sumx, _sumy;
    double _min;
    int _xmin, _ymin;
    double _max;
    int _xmax, _ymax;
    typename MaskedImageT::Mask::Pixel _bits;
};

/************************************************************************************************************/
/**
 * @brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintBits : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintBits(MaskedImageT const& mimage ///< The image the source lives in
                              ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                                  _bits(0) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _bits = 0x0;
    }
    virtual void reset(afwDetection::Footprint const&) {}

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


template<typename ExposureT>
void detail::checkFootprint(ExposurePatch<ExposureT>& patch, 
                            typename ExposureT::MaskedImageT::Mask::Pixel const bits // Bits in footprint
    ) {
    patch.setFlags(ExposurePatch<ExposureT>::NONE);
    
    // Check for bits set in the Footprint
    if (bits & ExposureT::MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
        patch.orFlag(ExposurePatch<ExposureT>::EDGE);
    }
    if (bits & ExposureT::MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        patch.orFlag(ExposurePatch<ExposureT>::INTERP);
    }
    if (bits & ExposureT::MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        patch.orFlag(ExposurePatch<ExposureT>::SAT);
    }

    // Check for bits set near the centroid
    afwGeom::Point2D const& center = patch.getCenter(); // Center in appropriate coordinates
    afwGeom::Point2I llc(afwImage::positionToIndex(center.getX()) - 1,
                         afwImage::positionToIndex(center.getY()) - 1);
    afwDet::Footprint const middle(afwGeom::BoxI(llc, afwGeom::ExtentI(3))); // central 3x3
    
    FootprintBits<typename ExposureT::MaskedImageT> bitsFunctor(patch.getExposure()->getMaskedImage());
    bitsFunctor.apply(middle);
    if (bitsFunctor.getBits() & ExposureT::MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        patch.orFlag(ExposurePatch<ExposureT>::INTERP_CENTER);
    }
    if (bitsFunctor.getBits() & ExposureT::MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        patch.orFlag(ExposurePatch<ExposureT>::SAT_CENTER);
    }
}



#if 0
#include "boost/graph/topological_sort.hpp"
// http://www.boost.org/doc/libs/1_47_0/libs/graph/doc/kevin_bacon.html
// http://www.boost.org/doc/libs/1_47_0/libs/graph/example/file_dependencies.cpp
// http://www.boost.org/doc/libs/1_47_0/libs/graph/doc/adjacency_list.html
PTR(std::vector<PTR(AlgorithmT)>) MeasureSources::topologicalSort() {
    if (sorted) {
        return sorted;
    }
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        boost::property<boost::vertex_name_t, std::string> > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef std::map<std::string, Vertex> NameVertexMap;

    Graph graph;
    NameVertexMap nvm;
    
    // Create a vertex for each algorithm
    for (typename std::vector<PTR(AlgorithmT)>::const_iterator algIter = _active.begin(); 
         algIter != _active.end(); ++algIter) {
        PTR(AlgorithmT) alg = *algIter;
        std::string const& name = alg->getName();

        NameVertexMap::iterator pos; 
        bool inserted;
        nvm.insert(std::make_pair(name, Vertex()));
    }

    // Link vertices with edges, defining the dependencies
    for (typename std::vector<PTR(AlgorithmT)>::const_iterator algIter = _active.begin(); 
         algIter != _active.end(); ++algIter) {
        PTR(AlgorithmT) alg = *algIter;
        std::string const& name = alg->getName();

        std::vector<std::string> const& deps = alg->getDependencies();
        for (std::vector<std::string>::const_iterator depIter = deps.begin();
             depIter != deps.end(); ++depIter) {
            boost::add_edge(nvm[name], nvm[*depIter], graph);
        }
    }

    // Topological sort to get correct ordering
    std::vector<Vertex> vertices(numAlgorithms);
    boost::topological_sort(graph, vertices.begin());

    // Translate to list of algorithms
    boost::property_map<Graph, vertex_name_t>::type algorithms = boost::get(boost::vertex_name, graph);

    sorted = boost::make_shared<std::vector<AlgorithmT> >();
    for (std::vector<Vertex>::const_iterator iter = vertices.begin(); iter != vertices.end(); ++iter) {
        sorted->push_back(algorithms[*iter]);
    }
   
    _sorted = sorted;

    return sorted;
}
#endif

/************************************************************************************************************/

static void nullAstrom(afwDet::Source& target, afwDet::Source const& source) {
    target.setXAstrom(source.getXAstrom());
    target.setYAstrom(source.getYAstrom());
    target.setFlagForDetection(target.getFlagForDetection() | Flags::PEAKCENTER);
}

/// Extractors to call the right extraction method 
struct ApPhotExtractor {
    typedef afwDet::Photometry MeasurementT;
    static std::string name() { return "source.apFlux"; }
    static PTR(afwDet::Measurement<MeasurementT>) measurements(afwDet::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDet::Source& source, afwDet::Photometry const& phot) {
        source.extractApPhotometry(phot);
    }
};
struct PsfPhotExtractor {
    typedef afwDet::Photometry MeasurementT;
    static std::string name() { return "source.psfFlux"; }
    static PTR(afwDet::Measurement<MeasurementT>) measurements(afwDet::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDet::Source& source, afwDet::Photometry const& phot) {
        source.extractPsfPhotometry(phot);
    }
};
struct ModelPhotExtractor {
    typedef afwDet::Photometry MeasurementT;
    static std::string name() { return "source.modelFlux"; }
    static PTR(afwDet::Measurement<MeasurementT>) measurements(afwDet::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDet::Source& source, afwDet::Photometry const& phot) {
        source.extractModelPhotometry(phot);
    }
};
struct InstPhotExtractor {
    typedef afwDet::Photometry MeasurementT;
    static std::string name() { return "source.instFlux"; }
    static PTR(afwDet::Measurement<MeasurementT>) measurements(afwDet::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDet::Source& source, afwDet::Photometry const& phot) {
        source.extractInstPhotometry(phot);
    }
};
struct AstrometryExtractor {
    typedef afwDet::Astrometry MeasurementT;
    static std::string name() { return "source.astrom"; }
    static PTR(afwDet::Measurement<MeasurementT>) measurements(afwDet::Source& source) {
        return source.getAstrometry();
    }
    static void extract(afwDet::Source& source, afwDet::Astrometry const& astrom) {
        source.extractAstrometry(astrom);
    }
};
struct ShapeExtractor {
    typedef afwDet::Shape MeasurementT;
    static std::string name() { return "source.shape"; }
    static PTR(afwDet::Measurement<MeasurementT>) measurements(afwDet::Source& source) {
        return source.getShape();
    }
    static void extract(afwDet::Source& source, afwDet::Shape const& shape) {
        source.extractShape(shape);
    }
};

/// Templated function to extract the correct measurement
template<class ExtractorT>
static void extractMeasurements(afwDet::Source& source,
                                pexPolicy::Policy const& policy)
{
    std::string const name = ExtractorT::name();
    if (policy.isString(name)) {
        std::string const& alg = policy.getString(name);
        if (alg != "NONE") {
            try {
                typename afwDet::Measurement<typename ExtractorT::MeasurementT>::TPtr meas = 
                    ExtractorT::measurements(source)->find(alg);
                if (!meas) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, 
                                      (boost::format("Can't find measurement from algorithm %s "
                                                     "to extract %s") % alg % name).str());
                }
                ExtractorT::extract(source, *meas);
            } catch (pexExceptions::Exception& e) {
                LSST_EXCEPT_ADD(e, (boost::format("Extracting algorithm %s for %s") % 
                                    alg % name).str());
                throw e;
            }
        }
    }
}



/**
 * Use *this to measure the Footprint foot, setting fields in src
 */
template<typename ExposureT>
template<class Measurer> 
void MeasureSources<ExposureT>::_measure(
    afwDet::Source& target,
    afwDet::Source& source, 
    afwImage::Wcs const& wcs,
    typename Measurer::ExposureContainerT& patches
    ) {

    // For each algorithm in topological order
    //  try
    //   measure
    //  catch: set null
    //   add result
    // Extract astrometry, photometry, shape
    // Set star/galaxy

    Measurer::check(source, patches);

    // Convert to const 
    typename Measurer::ConstExposureContainerT& constPatches = Measurer::constify(patches);

    // Astrometry
    if (!getMeasureAstrom()) {
        nullAstrom(target, source);
    } else {
        PTR(afwDet::Astrometry) astrom = 
            Measurer::template measure<afwDet::Astrometry>(getMeasureAstrom(), target, source, constPatches);
        target.setAstrometry(astrom);
        extractMeasurements<AstrometryExtractor>(target, _policy);
        if (lsst::utils::isnan(target.getXAstrom()) || lsst::utils::isnan(target.getYAstrom())) {
            nullAstrom(target, source);
        } else {
            Measurer::updateAstrom(target, patches);
        }
    }

    // Shapes
    if (getMeasureShape()) {
        PTR(afwDet::Shape) shapes =
            Measurer::template measure<afwDet::Shape>(getMeasureShape(), target, source, constPatches);
        target.setShape(shapes);
        extractMeasurements<ShapeExtractor>(target, _policy);
    }

    // Photometry
    if (getMeasurePhotom()) {
        PTR(afwDet::Photometry) phot = 
            Measurer::template measure<afwDet::Photometry>(getMeasurePhotom(), target, source, constPatches);
        target.setPhotometry(phot);
        extractMeasurements<ApPhotExtractor>(target, _policy);
        extractMeasurements<PsfPhotExtractor>(target, _policy);
        extractMeasurements<ModelPhotExtractor>(target, _policy);
        extractMeasurements<InstPhotExtractor>(target, _policy);

        // Set photometry flags
        boost::int64_t flag = target.getFlagForDetection();
        for (afwDet::Measurement<afwDet::Photometry>::const_iterator i = phot->begin();
             i != phot->end(); ++i) {
            flag |= (*i)->getFlag();
        }
        target.setFlagForDetection(flag);

        // Add some star/galaxy information.  The "extendedness" parameter is supposed to be the
        // probability of being extended
        std::vector<float> fac(3);// Fiddle factors for star/galaxy separation
        fac[0] = _policy.getDouble("classification.sg_fac1");
        fac[1] = _policy.getDouble("classification.sg_fac2");
        fac[2] = _policy.getDouble("classification.sg_fac3");

        bool const isStar = ((fac[0]*target.getInstFlux() + fac[1]*target.getInstFluxErr()) <
                             (target.getPsfFlux() + fac[2]*target.getPsfFluxErr()) ? 0.0 : 1.0);
#if 0
        target.setExtendedness(isStar ? 0.0 : 1.0);
#else
        target.setApDia(isStar ? 0.0 : 1.0);
#endif
    }

    // Translate ExposurePatch flags to Source flags.
    target.setFlagForDetection(target.getFlagForDetection() |
                               ExposurePatch<ExposureT>::sourceFlags(Measurer::flags(constPatches)));

}

//
// Explicit instantiations
//
// \cond
template class MeasureSources<afwImage::Exposure<float> >;
template class MeasureSources<afwImage::Exposure<int> >;
// \endcond
}}}
