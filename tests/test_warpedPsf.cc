/*
 * This file is part of meas_algorithms.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DISTORTION
#include <memory>
#include <random>

#include <boost/test/unit_test.hpp>
#include "astshim.h"

#include "lsst/afw/table/io/Persistable.cc"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/meas/algorithms/WarpedPsf.h"

using namespace std;
using namespace Eigen;
using namespace lsst::geom;
using namespace lsst::afw::math;
using namespace lsst::afw::image;
using namespace lsst::afw::detection;
using namespace lsst::afw::geom::ellipses;
using namespace lsst::meas::algorithms;

static std::mt19937 rng(0);  // RNG deliberately initialized with same seed every time
static std::uniform_int_distribution<> uni_int(0, 100);
static std::uniform_real_distribution<> uni_double(0.0, 1.0);

// -------------------------------------------------------------------------------------------------
//
// Helper functions

static inline Point2D randpt() {
    // returns randomly located point in [-100,100] x [-100,100]
    return Point2D(200 * uni_double(rng) - 100, 200 * uni_double(rng) - 100);
}

static inline double dist(const Point2D &p1, const Point2D &p2) {
    double dx = p1.getX() - p2.getX();
    double dy = p1.getY() - p2.getY();
    return sqrt(dx * dx + dy * dy);
}

static inline double dist(const AffineTransform &a1, const AffineTransform &a2) {
    double ret = 0.0;
    for (int i = 0; i < 6; i++) ret += (a1[i] - a2[i]) * (a1[i] - a2[i]);
    return sqrt(ret);
}

static inline double compare(const Image<double> &im1, const Image<double> &im2) {
    assert(im1.getWidth() == im2.getWidth());
    assert(im1.getHeight() == im2.getHeight());
    assert(im1.getX0() == im2.getX0());
    assert(im1.getY0() == im2.getY0());

    double t11 = 0.0;
    double t12 = 0.0;
    double t22 = 0.0;

    int nx = im1.getWidth();
    int ny = im1.getHeight();

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = im1(i, j);
            double y = im2(i, j);
            t11 += x * x;
            t12 += (x - y) * (x - y);
            t22 += y * y;
        }
    }

    assert(t11 > 0.0);
    assert(t22 > 0.0);
    return sqrt(fabs(t12) / sqrt(t11 * t22));
}

//
// Make a transform with the following form, using random "reasonable" coefficients:
//   x' = x + Ax + By + Cx^2 + Dxy + Ey^2
//   y' = y + Fx + Gy + Hx^2 + Ixy + Jy^2
//
std::shared_ptr<lsst::afw::geom::TransformPoint2ToPoint2> makeRandomToyTransform() {
    double A = 0.1 * (uni_double(rng) - 0.5);
    double B = 0.1 * (uni_double(rng) - 0.5);
    double C = 0.0001 * (uni_double(rng) - 0.5);
    double D = 0.0001 * (uni_double(rng) - 0.5);
    double E = 0.0001 * (uni_double(rng) - 0.5);
    double F = 0.1 * (uni_double(rng) - 0.5);
    double G = 0.1 * (uni_double(rng) - 0.5);
    double H = 0.0001 * (uni_double(rng) - 0.5);
    double I = 0.0001 * (uni_double(rng) - 0.5);
    double J = 0.0001 * (uni_double(rng) - 0.5);

    // clang-format off
    // each 4 entries are: coefficient, output index, power of x, power of y
    std::vector<double> const coeffVec = {
        1 + A, 1, 1, 0,
        B, 1, 0, 1,
        C, 1, 2, 0,
        D, 1, 1, 1,
        E, 1, 0, 2,

        F, 2, 1, 0,
        1 + G, 2, 0, 1,
        H, 2, 2, 0,
        I, 2, 1, 1,
        J, 2, 0, 2
    };
    // clang-format on
    int const nOut = 2;
    int const nCoeffs = coeffVec.size() / (nOut + 2);
    auto const coeffArr = ast::arrayFromVector(coeffVec, nCoeffs);
    auto const mapping = ast::PolyMap(coeffArr, nOut, "IterInverse=1, TolInverse=1e-8, NIterInverse=20");
    return std::make_shared<lsst::afw::geom::TransformPoint2ToPoint2>(mapping);
}

// -------------------------------------------------------------------------------------------------
//
// ToyPsf: general PDF of the form
//   exp(-ax^2/2 - bxy - cy^2/2)
//
// where
//   a = 0.1 (1 + Ax + By)
//   b = 0.1 (Cx + Dy)
//   c = 0.1 (1 + Ex + Fy)
//

//
// Helper function which fills an image with a normalized 2D Gaussian of the form
//   exp(-a(x-px)^2/2 - b(x-px)(y-py) - c(y-py)^2/2)
//
static std::shared_ptr<Image<double>> fill_gaussian(double a, double b, double c, double px, double py,
                                                    int nx, int ny, int x0, int y0) {
    // smallest eigenvalue
    double lambda = 0.5 * (a + c + sqrt((a - c) * (a - c) + b * b));

    // approximate size of box needed to hold kernel
    double width = sqrt(5. / lambda);

    assert(lambda > 1.0e-10);
    assert(x0 - px <= -width && x0 - px + nx - 1 >= width);
    assert(y0 - py <= -width && y0 - py + ny - 1 >= width);

    auto im = std::make_shared<Image<double>>(nx, ny);
    im->setXY0(x0, y0);

    double imSum = 0.0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = i + x0 - px;
            double y = j + y0 - py;
            double t = exp(-0.5 * a * x * x - b * x * y - 0.5 * c * y * y);
            (*im)(i, j) = t;
            imSum += t;
        }
    }

    (*im) /= imSum;
    return im;
}

struct ToyPsf : public ImagePsf {
    double _A, _B, _C, _D, _E, _F, _ksize;

    ToyPsf(double A, double B, double C, double D, double E, double F, int ksize)
            : _A(A), _B(B), _C(C), _D(D), _E(E), _F(F), _ksize(ksize) {}

    virtual ~ToyPsf() {}

    std::shared_ptr<Psf> clone() const override {
        return std::make_shared<ToyPsf>(_A, _B, _C, _D, _E, _F, _ksize);
    }

    std::shared_ptr<Psf> resized(int width, int height) const override {
        return std::make_shared<ToyPsf>(_A, _B, _C, _D, _E, _F, width);
    }

    void evalABC(double &a, double &b, double &c, Point2D const &p) const {
        double x = p.getX();
        double y = p.getY();

        a = 0.1 * (1.0 + _A * x + _B * y);
        b = 0.1 * (_C * x + _D * y);
        c = 0.1 * (1.0 + _E * x + _F * y);
    }

    Box2I doComputeBBox(Point2D const &, Color const &) const override {
        return Box2I(Point2I(-_ksize, -_ksize), Extent2I(2 * _ksize + 1, 2 * _ksize + 1));
    }

    std::shared_ptr<Image> doComputeKernelImage(Point2D const &ccdXY, Color const &) const override {
        double a, b, c;
        this->evalABC(a, b, c, ccdXY);

        Box2I bbox = computeBBox(getAveragePosition());
        return fill_gaussian(a, b, c, 0, 0, bbox.getWidth(), bbox.getHeight(), bbox.getMinX(),
                             bbox.getMinY());
    }

    // factory function
    static std::shared_ptr<ToyPsf> makeRandom(int ksize) {
        double A = 0.005 * (uni_double(rng) - 0.5);
        double B = 0.005 * (uni_double(rng) - 0.5);
        double C = 0.005 * (uni_double(rng) - 0.5);
        double D = 0.005 * (uni_double(rng) - 0.5);
        double E = 0.005 * (uni_double(rng) - 0.5);
        double F = 0.005 * (uni_double(rng) - 0.5);

        return std::make_shared<ToyPsf>(A, B, C, D, E, F, ksize);
    }

    bool isPersistable() const noexcept override { return true; }

    std::string getPersistenceName() const override { return "ToyPsf"; }

    std::string getPythonModule() const override { return "lsst.meas.algorithms"; };

    void write(OutputArchiveHandle &handle) const override {
        lsst::afw::table::Schema schema;
        lsst::afw::table::Key<double> keyA = schema.addField<double>("A", "");
        lsst::afw::table::Key<double> keyB = schema.addField<double>("B", "");
        lsst::afw::table::Key<double> keyC = schema.addField<double>("C", "");
        lsst::afw::table::Key<double> keyD = schema.addField<double>("D", "");
        lsst::afw::table::Key<double> keyE = schema.addField<double>("E", "");
        lsst::afw::table::Key<double> keyF = schema.addField<double>("F", "");
        lsst::afw::table::Key<double> keySize = schema.addField<double>("Size", "");

        lsst::afw::table::BaseCatalog catalog = handle.makeCatalog(schema);
        std::shared_ptr<lsst::afw::table::BaseRecord> record = catalog.addNew();
        record->set(keyA, _A);
        record->set(keyB, _B);
        record->set(keyC, _C);
        record->set(keyD, _D);
        record->set(keyE, _E);
        record->set(keyF, _F);
        record->set(keySize, _ksize);
        handle.saveCatalog(catalog);
    }
};

namespace {

// A PersistableFactory for ToyPsf
class : public lsst::afw::table::io::PersistableFactory {
public:
    std::shared_ptr<lsst::afw::table::io::Persistable> read(
            lsst::afw::table::io::InputArchive const &archive,
            lsst::afw::table::io::CatalogVector const &catalogs) const override {
        lsst::afw::table::Schema schema;
        lsst::afw::table::Key<double> keyA = schema.addField<double>("A", "");
        lsst::afw::table::Key<double> keyB = schema.addField<double>("B", "");
        lsst::afw::table::Key<double> keyC = schema.addField<double>("C", "");
        lsst::afw::table::Key<double> keyD = schema.addField<double>("D", "");
        lsst::afw::table::Key<double> keyE = schema.addField<double>("E", "");
        lsst::afw::table::Key<double> keyF = schema.addField<double>("F", "");
        lsst::afw::table::Key<double> keySize = schema.addField<double>("Size", "");

        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        lsst::afw::table::BaseRecord const &record = catalogs.front().front();
        LSST_ARCHIVE_ASSERT(record.getSchema() == schema);
        return std::make_shared<ToyPsf>(record.get(keyA), record.get(keyB), record.get(keyC),
                                        record.get(keyD), record.get(keyE), record.get(keyF),
                                        record.get(keySize));
    }

    using lsst::afw::table::io::PersistableFactory::PersistableFactory;
} registration("ToyPsf");

}  // namespace

namespace lsst::afw::table::io {

template std::shared_ptr<ToyPsf> PersistableFacade<ToyPsf>::dynamicCast(std::shared_ptr<Persistable> const &);

}  // namespace lsst::afw::table::io

BOOST_AUTO_TEST_CASE(warpedPsf) {
    auto distortion = makeRandomToyTransform();

    std::shared_ptr<ToyPsf> unwarped_psf = ToyPsf::makeRandom(100);
    std::shared_ptr<WarpedPsf> warped_psf = std::make_shared<WarpedPsf>(unwarped_psf, distortion);

    Point2D p = randpt();
    Point2D q = distortion->applyInverse(p);
    // Check inverse at this point
    BOOST_CHECK(dist(distortion->applyForward(q), p) < 1e-7);

    // warped image
    std::shared_ptr<Image<double>> im = warped_psf->computeImage(p);
    int nx = im->getWidth();
    int ny = im->getHeight();
    int x0 = im->getX0();
    int y0 = im->getY0();

    double a, b, c;
    unwarped_psf->evalABC(a, b, c, q);

    Eigen::Matrix2d m0;
    // clang-format off
    m0 << a, b,
          b, c;
    // clang-format on

    AffineTransform atr = lsst::afw::geom::linearizeTransform(*distortion->inverted(), p);

    Eigen::Matrix2d md;
    // clang-format off
    md << atr.getLinear()[0], atr.getLinear()[2],
          atr.getLinear()[1], atr.getLinear()[3];  // LinearTransform transposed index convention
    // clang-format on

    Eigen::Matrix2d m1 = md.transpose() * m0 * md;

    // this should be the same as the warped image, up to artifacts from warping/pixelization
    std::shared_ptr<Image<double>> im2 =
            fill_gaussian(m1(0, 0), m1(0, 1), m1(1, 1), p.getX(), p.getY(), nx, ny, x0, y0);

    // TODO: improve this test; the ideal thing would be to repeat with
    // finer resolutions and more stringent threshold
    BOOST_CHECK(compare(*im, *im2) < 0.006);

    // Check that computeBBox returns same dimensions as image
    BOOST_CHECK(warped_psf->computeBBox(p).getWidth() == nx);
    BOOST_CHECK(warped_psf->computeBBox(p).getHeight() == ny);
}

// Test that WarpedPsf properly pads original Psf images before warping, so that
// the warped image extends all the way to the edges.
// Because afw.math.warpImage will set unfilled pixels to exactly 0 by default,
// this test case checks that each of the 4 edges of a WarpedPsf is non-zero.
BOOST_AUTO_TEST_CASE(warpedPsfPadding) {
    auto distortion = makeRandomToyTransform();

    // Make psf with small kernel size  so that lack of padding is more apparent
    std::shared_ptr<ToyPsf> unwarped_psf = ToyPsf::makeRandom(7);
    std::shared_ptr<WarpedPsf> warped_psf = std::make_shared<WarpedPsf>(unwarped_psf, distortion);
    std::shared_ptr<Image<double>> warpedImage = warped_psf->computeKernelImage(Point2D(-10., 150.));

    // The outer columns and rows must test non-zero
    // Tolerance should be very low, because edges of small PSFs with large BBoxes
    // can legitimately have pixel values on order of minimum subnormal numbers (1e-323).
    // Tolerance may be as low as zero (which has an exact representation).
    const double zero = 0.;
    double sumRow = 0.;

    // Check first and last row
    for (const int &y : {0, warpedImage->getHeight() - 1}) {
        sumRow = 0.;
        for (Image<double>::x_iterator ptr = warpedImage->row_begin(y), end = warpedImage->row_end(y);
             ptr != end; ++ptr) {
            sumRow += *ptr;
        }
        BOOST_CHECK(std::abs(sumRow) > zero);
    }

    // Check first and last column
    for (const int &x : {0, warpedImage->getWidth() - 1}) {
        sumRow = 0.;
        for (Image<double>::y_iterator ptr = warpedImage->col_begin(x), end = warpedImage->col_end(x);
             ptr != end; ++ptr) {
            sumRow += *ptr;
        }
        BOOST_CHECK(std::abs(sumRow) > zero);
    }
}

template <class T>
std::shared_ptr<T> roundTrip(std::shared_ptr<T> original) {
    // writeFits(MemFileManager) and writeFits(string) are non-overrideable methods
    // that delegate to the same WarpedPsf code, so we only need to test one.
    lsst::afw::fits::MemFileManager tempFile;
    original->writeFits(tempFile);
    return lsst::afw::table::io::PersistableFacade<T>::readFits(tempFile);
}

BOOST_AUTO_TEST_CASE(persistence) {
    auto distortion = makeRandomToyTransform();

    std::shared_ptr<ToyPsf> unwarped_psf = ToyPsf::makeRandom(100);
    std::shared_ptr<WarpedPsf> warped_psf = std::make_shared<WarpedPsf>(unwarped_psf, distortion);

    // If ToyPsf's persistence is buggy, can't test WarpedPsf
    BOOST_TEST_REQUIRE(unwarped_psf->isPersistable());
    std::shared_ptr<ToyPsf> copy_unwarped_psf = roundTrip(unwarped_psf);
    for (int i = 0; i < 10; ++i) {
        Point2D p = randpt();
        // Threshold copied from warpedPsf test, no idea where it comes from
        //BOOST_TEST_REQUIRE(compare(*(unwarped_psf->computeImage(p)), *(copy_unwarped_psf->computeImage(p))) <
        //                   0.006);
    }

    //BOOST_TEST_REQUIRE(warped_psf->isPersistable());
    std::shared_ptr<WarpedPsf> copy_warped_psf = roundTrip(warped_psf);

    for (int i = 0; i < 10; ++i) {
        Point2D p = randpt();
        // Threshold copied from warpedPsf test, no idea where it comes from
        BOOST_TEST(compare(*(warped_psf->computeImage(p)), *(copy_warped_psf->computeImage(p))) < 0.006);
    }
}
