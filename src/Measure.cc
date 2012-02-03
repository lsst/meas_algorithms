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

namespace pexPolicy = lsst::pex::policy;
namespace pexExceptions = lsst::pex::exceptions;    
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;


namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/*
 * Return the numeric value of name as double
 */
double getNumeric(lsst::pex::policy::Policy const& policy, std::string const& name)
{
    return policy.isDouble(name) ? policy.getDouble(name) : policy.getInt(name);
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


/// Check footprint for bad pixels
template<typename ExposureT>
void checkFootprint(ExposurePatch<ExposureT>& patch,                         // Patch to check
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
    afwDetection::Footprint const middle(afwGeom::BoxI(llc, afwGeom::ExtentI(3))); // central 3x3
    
    FootprintBits<typename ExposureT::MaskedImageT> bitsFunctor(patch.getExposure()->getMaskedImage());
    bitsFunctor.apply(middle);
    if (bitsFunctor.getBits() & ExposureT::MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        patch.orFlag(ExposurePatch<ExposureT>::INTERP_CENTER);
    }
    if (bitsFunctor.getBits() & ExposureT::MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        patch.orFlag(ExposurePatch<ExposureT>::SAT_CENTER);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions to cover different measurement scenarios
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Measuring sources on the same image on which they were detected.
template<typename ExposureT>
struct SingleMeasurer {
    typedef ExposurePatch<ExposureT> ExposureContainerT;
    typedef ExposurePatch<ExposureT> const ConstExposureContainerT;

    /// Check pixels in the footprint, setting appropriate flags, and get a rough starting x,y position
    static void check(lsst::afw::detection::Source& source, ExposurePatch<ExposureT>& patch) {
        FootprintCentroid<typename ExposureT::MaskedImageT> centroider(patch.getExposure()->getMaskedImage());
        centroider.apply(*patch.getFootprint());
        double const x = centroider.getX();
        double const y = centroider.getY();
        source.setXAstrom(x);
        source.setYAstrom(y);
        patch.setCenter(lsst::afw::geom::Point2D(x, y));
        checkFootprint(patch, centroider.getBits());
    }

    /// Make the exposure container carry const members
    static ExposurePatch<ExposureT> const& constify(ExposurePatch<ExposureT> const& patch) {
        /// No change needed in this case
        return patch;
    }

    /// Make the measurement
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(CONST_PTR(MeasureQuantity<MeasurementT, ExposureT>) measurer,
                                     lsst::afw::detection::Source& target, 
                                     lsst::afw::detection::Source const& source,
                                     ExposurePatch<ExposureT> const& patch) {
        return measurer->measureSingle(target, source, patch);
    }

    /// Execute the algorithm
    template<typename MeasurementT>
    static PTR(MeasurementT) algorithm(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                       lsst::afw::detection::Source const& target,
                                       lsst::afw::detection::Source const& source,
                                       ExposurePatch<ExposureT> const& patch) {
        return alg->measureSingle(target, source, patch);
    }                                       

    /// Update the astrometry
    static void updateAstrom(lsst::afw::detection::Source const& target, ExposurePatch<ExposureT>& patch) {
        lsst::afw::geom::Point2D const center(target.getXAstrom(), target.getYAstrom());
        patch.setCenter(patch.fromStandard()(center));
    }

    /// Get combined flags
    static boost::int64_t flags(ExposurePatch<ExposureT> const& patch) {
        return patch.getFlags();
    }
};

/// Measuring a single source on multiple images
template<typename ExposureT>
struct MultipleMeasurer {
    typedef std::vector<PTR(ExposurePatch<ExposureT>)> ExposureContainerT;
    typedef std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const ConstExposureContainerT;

    /// Check pixels in the footprint, setting appropriate flags, and get a rough starting x,y position
    static void check(lsst::afw::detection::Source& source, 
                      std::vector<PTR(ExposurePatch<ExposureT>)>& patches)
    {
        for (size_t i = 0; i < patches.size(); ++i) {
            PTR(ExposurePatch<ExposureT>) p = patches[i];
            FootprintBits<typename ExposureT::MaskedImageT> bitsFunctor(p->getExposure()->getMaskedImage());
            bitsFunctor.apply(*p->getFootprint());
            checkFootprint(*p, bitsFunctor.getBits());
        }
    }

    /// Make the exposure container carry const members
    static std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const&
    constify(std::vector<PTR(ExposurePatch<ExposureT>)> const& patches) {
        /// The compiler doesn't know how to automatically convert
        /// std::vector<PTR(T)> to std::vector<CONST_PTR(T)> because the way the
        /// template system works means that in theory the two may be
        /// specialised differently.  This is an explicit conversion.
        ///
        /// see e.g., http://stackoverflow.com/questions/2102244/vector-and-const
        return reinterpret_cast<std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const&>(patches);
    }

    /// Make the measurement
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(CONST_PTR(MeasureQuantity<MeasurementT, ExposureT>) measurer,
                                     lsst::afw::detection::Source& target, 
                                     lsst::afw::detection::Source const& source,
                                     std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        return measurer->measureMultiple(target, source, patches);
    }
    
    /// Execute the algorithm
    template<typename MeasurementT>
    static PTR(MeasurementT) algorithm(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                       lsst::afw::detection::Source const& target, 
                                       lsst::afw::detection::Source const& source,
                                       std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        return alg->measureMultiple(target, source, patches);
    }                                       

    /// Update the astrometry
    static void updateAstrom(lsst::afw::detection::Source const& target,
                             std::vector<PTR(ExposurePatch<ExposureT>)>& patches) {
        lsst::afw::geom::Point2D const center(target.getXAstrom(), target.getYAstrom());
        for (size_t i = 0; i < patches.size(); ++i) {
            patches[i]->setCenter(patches[i]->fromStandard()(center));
        }
    }

    /// Get combined flags
    static boost::int64_t flags(std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        boost::int64_t eFlags = ExposurePatch<ExposureT>::ALL;
        for (typename ConstExposureContainerT::const_iterator iter = patches.begin(); 
             iter != patches.end(); ++iter) {
            eFlags &= (*iter)->getFlags();
        }
        return eFlags;
    }
};



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

/// Set 'null' astrometry values
void nullAstrom(afwDetection::Source& target, afwDetection::Source const& source) {
    target.setXAstrom(source.getXAstrom());
    target.setYAstrom(source.getYAstrom());
    target.setFlagForDetection(target.getFlagForDetection() | Flags::PEAKCENTER);
}

/// Extractors to call the right extraction method 
struct ApPhotExtractor {
    typedef afwDetection::Photometry MeasurementT;
    static std::string name() { return "source.apFlux"; }
    static PTR(afwDetection::Measurement<MeasurementT>) measurements(afwDetection::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDetection::Source& source, afwDetection::Photometry const& phot) {
        source.extractApPhotometry(phot);
    }
};
struct PsfPhotExtractor {
    typedef afwDetection::Photometry MeasurementT;
    static std::string name() { return "source.psfFlux"; }
    static PTR(afwDetection::Measurement<MeasurementT>) measurements(afwDetection::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDetection::Source& source, afwDetection::Photometry const& phot) {
        source.extractPsfPhotometry(phot);
    }
};
struct ModelPhotExtractor {
    typedef afwDetection::Photometry MeasurementT;
    static std::string name() { return "source.modelFlux"; }
    static PTR(afwDetection::Measurement<MeasurementT>) measurements(afwDetection::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDetection::Source& source, afwDetection::Photometry const& phot) {
        source.extractModelPhotometry(phot);
    }
};
struct InstPhotExtractor {
    typedef afwDetection::Photometry MeasurementT;
    static std::string name() { return "source.instFlux"; }
    static PTR(afwDetection::Measurement<MeasurementT>) measurements(afwDetection::Source& source) {
        return source.getPhotometry();
    }
    static void extract(afwDetection::Source& source, afwDetection::Photometry const& phot) {
        source.extractInstPhotometry(phot);
    }
};
struct AstrometryExtractor {
    typedef afwDetection::Astrometry MeasurementT;
    static std::string name() { return "source.astrom"; }
    static PTR(afwDetection::Measurement<MeasurementT>) measurements(afwDetection::Source& source) {
        return source.getAstrometry();
    }
    static void extract(afwDetection::Source& source, afwDetection::Astrometry const& astrom) {
        source.extractAstrometry(astrom);
    }
};
struct ShapeExtractor {
    typedef afwDetection::Shape MeasurementT;
    static std::string name() { return "source.shape"; }
    static PTR(afwDetection::Measurement<MeasurementT>) measurements(afwDetection::Source& source) {
        return source.getShape();
    }
    static void extract(afwDetection::Source& source, afwDetection::Shape const& shape) {
        source.extractShape(shape);
    }
};

/// Templated function to extract the correct measurement
template<class ExtractorT>
void extractMeasurements(afwDetection::Source& source,
                         pexPolicy::Policy const& policy
    )
{
    std::string const name = ExtractorT::name();
    if (policy.isString(name)) {
        std::string const& alg = policy.getString(name);
        if (alg != "NONE") {
            try {
                typename afwDetection::Measurement<typename ExtractorT::MeasurementT>::TPtr meas = 
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

/// Engine for measuring a source
template<class Measurer, typename ExposureT>
void doMeasure(
    MeasureSources<ExposureT> const& ms,
    afwDetection::Source& target,
    afwDetection::Source& source, 
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
    if (!ms.getMeasureAstrom()) {
        nullAstrom(target, source);
    } else {
        PTR(afwDetection::Astrometry) astrom = 
            Measurer::template measure<afwDetection::Astrometry>(ms.getMeasureAstrom(), target,
                                                           source, constPatches);
        target.setAstrometry(astrom);
        extractMeasurements<AstrometryExtractor>(target, ms.getPolicy());
        if (lsst::utils::isnan(target.getXAstrom()) || lsst::utils::isnan(target.getYAstrom())) {
            nullAstrom(target, source);
        } else {
            Measurer::updateAstrom(target, patches);
        }
    }

    // Shapes
    if (ms.getMeasureShape()) {
        PTR(afwDetection::Shape) shapes =
            Measurer::template measure<afwDetection::Shape>(ms.getMeasureShape(), target, source, 
                                                            constPatches);
        target.setShape(shapes);
        extractMeasurements<ShapeExtractor>(target, ms.getPolicy());
    }

    // Photometry
    if (ms.getMeasurePhotom()) {
        PTR(afwDetection::Photometry) phot = 
            Measurer::template measure<afwDetection::Photometry>(ms.getMeasurePhotom(), target,
                                                           source, constPatches);
        target.setPhotometry(phot);
        extractMeasurements<ApPhotExtractor>(target, ms.getPolicy());
        extractMeasurements<PsfPhotExtractor>(target, ms.getPolicy());
        extractMeasurements<ModelPhotExtractor>(target, ms.getPolicy());
        extractMeasurements<InstPhotExtractor>(target, ms.getPolicy());

        // Set photometry flags
        boost::int64_t flag = target.getFlagForDetection();
        for (afwDetection::Measurement<afwDetection::Photometry>::const_iterator i = phot->begin();
             i != phot->end(); ++i) {
            flag |= (*i)->getFlag();
        }
        target.setFlagForDetection(flag);

        // Add some star/galaxy information.  The "extendedness" parameter is supposed to be the
        // probability of being extended
        std::vector<float> fac(3);// Fiddle factors for star/galaxy separation
        fac[0] = ms.getPolicy().getDouble("classification.sg_fac1");
        fac[1] = ms.getPolicy().getDouble("classification.sg_fac2");
        fac[2] = ms.getPolicy().getDouble("classification.sg_fac3");

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

} // anonymous namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MeasureSources implementation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename ExposureT>
MeasureSources<ExposureT>::MeasureSources(pexPolicy::Policy const& policy) :
    _policy( policy),
    _moLog(pexLogging::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                           pexLogging::Log::INFO)),
    _measureAstrom(boost::make_shared<MeasureAstrometryT>()),
    _measurePhotom(boost::make_shared<MeasurePhotometryT>()),
    _measureShape(boost::make_shared<MeasureShapeT>())
{}

template<typename ExposureT>
void MeasureSources<ExposureT>::measure(afwDetection::Source& target, CONST_PTR(ExposureT) exp) const
{
    CONST_PTR(afwImage::Wcs) wcs = exp->getWcs();
    CONST_PTR(afwDetection::Footprint) foot = target.getFootprint();
    bool negative = target.getFlagForDetection() & Flags::DETECT_NEGATIVE;
    // Get highest peak
    afwDetection::Footprint::PeakList const& peakList = foot->getPeaks();
    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, 
                          (boost::format("No peak for source %d") % target.getId()).str());
    }
    PTR(afwDetection::Peak) peak = peakList[0];
    for (size_t i = 1; i < peakList.size(); ++i) {
        float value = peakList[i]->getPeakValue();
        if (negative) {
            value *= -1;
        }
        if (value > peak->getPeakValue()) {
            peak = peakList[i];
        }
    }
    afwGeom::Point2D center(peak->getFx(), peak->getFy());
    ExposurePatch<ExposureT> patch(exp, foot, center);
    doMeasure<SingleMeasurer<ExposureT> >(*this, target, target, patch);
}

template<typename ExposureT>
void MeasureSources<ExposureT>::measure(afwDetection::Source& target, CONST_PTR(ExposureT) exp, 
                                        afwGeom::Point2D const& center) const
{
    CONST_PTR(afwImage::Wcs) wcs = exp->getWcs();
    ExposurePatch<ExposureT> patch(exp, target.getFootprint(), center);
    doMeasure<SingleMeasurer<ExposureT> >(*this, target, target, patch);
}

template<typename ExposureT>
void MeasureSources<ExposureT>::measure(afwDetection::Source& target, afwDetection::Source const& source,
                                        afwImage::Wcs const& wcs, CONST_PTR(ExposureT) exp) const
{
    std::vector<CONST_PTR(ExposureT)> exposures(1);
    exposures[0] = exp;
    measure(target, source, wcs, exposures);
}

template<typename ExposureT>
void MeasureSources<ExposureT>::measure(
    afwDetection::Source& target,
    afwDetection::Source const& source,
    afwImage::Wcs const& wcs,
    std::vector<CONST_PTR(ExposureT)> const& exposures
    ) const
{
    size_t size = exposures.size();
    std::vector<PTR(ExposurePatch<ExposureT>)> patches(size);
    afwGeom::Point2D center(source.getXAstrom(), source.getYAstrom());
    for (size_t i = 0; i != size; ++i) {
        CONST_PTR(afwImage::Wcs) expWcs = exposures[i]->getWcs();
        if (!expWcs) {
            throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, 
                              (boost::format("No WCS for exposure %d") % i).str());
        }
        patches[i] = makeExposurePatch(exposures[i], *source.getFootprint(), center, wcs);
    }
    doMeasure<MultipleMeasurer<ExposureT> >(*this, target, const_cast<afwDetection::Source&>(source), 
                                            patches);
}



//
// Explicit instantiations
//
// \cond

#define INSTANTIATE(PIXEL) \
    template class MeasureSources<afwImage::Exposure<PIXEL> >; \
    template PTR(MeasureSources<afwImage::Exposure<PIXEL> >) \
        makeMeasureSources(afwImage::Exposure<PIXEL> const&, pexPolicy::Policy const&);

INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);

// \endcond
}}}
