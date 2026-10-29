# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import numpy as np
import tempfile

from lsst.meas.algorithms import stamps
from lsst.afw import image as afwImage
from lsst.afw.fits import Fits, readMetadata
from lsst.afw.geom.testUtils import TransformTestBaseClass
from lsst.daf.base import PropertyList, PropertySet
import lsst.geom as geom
import lsst.afw.geom.transformFactory as tF
import lsst.utils.tests

_RNG = np.random.Generator(np.random.MT19937(5))


def image_hdu_metadata(filename):
    """Read the header of every image extension, ordered by stamp.

    Parameters
    ----------
    filename : `str`
        Name of the file to read.

    Returns
    -------
    headers : `list` [`lsst.daf.base.PropertyList`]
        The header of each image extension, in stamp order.
    """
    with Fits(filename, 'r') as f:
        nExtensions = f.countHdus()
    headers = {}
    for idx in range(1, nExtensions):
        md = readMetadata(filename, hdu=idx)
        if md['EXTNAME'] == 'IMAGE':
            headers[md['EXTVER']] = md
    return [headers[extver] for extver in sorted(headers)]


def make_stamps(n_stamps=3, use_archive=False):
    stampSize = 25
    # create dummy stamp stamps
    stampImages = [afwImage.MaskedImageF(stampSize, stampSize)
                   for _ in range(n_stamps)]
    for stampIm in stampImages:
        stampImArray = stampIm.image.array
        stampImArray += _RNG.random((stampSize, stampSize))
        stampMaskArray = stampIm.mask.array
        stampMaskArray += 10
        stampVarArray = stampIm.variance.array
        stampVarArray += 1000.
    ras = _RNG.random(n_stamps)*360.
    decs = _RNG.random(n_stamps)*180 - 90
    archive_elements = [tF.makeTransform(geom.AffineTransform(_RNG.random(2))) if use_archive else None
                        for _ in range(n_stamps)]

    metadata = PropertyList()
    metadata['RA_DEG'] = ras
    metadata['DEC_DEG'] = decs

    stamp_list = [stamps.Stamp(stamp_im=stampIm,
                               position=geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                                         geom.Angle(dec, geom.degrees)),
                               archive_element=ae,
                               metadata=metadata)
                  for stampIm, ra, dec, ae in zip(stampImages, ras, decs, archive_elements)]

    return stamps.Stamps(stamp_list, metadata=metadata, use_archive=True)


class StampsBaseTestCase(lsst.utils.tests.TestCase):
    """Test StampsBase.
    """
    def testReadFitsWithOptionsNotImplementedErrorRaised(self):
        """
        Test that subclasses have their own version
        of this implemented or an error is raised.
        """
        class FakeStampsBase(stamps.StampsBase):
            def __init__(self):
                return

        with self.assertRaises(NotImplementedError):
            FakeStampsBase.readFitsWithOptions('noFile', {})

    def testReadFitsWithOptionsMetadataError(self):
        """Test that error is raised when STAMPCLS returns None
        """
        with tempfile.NamedTemporaryFile() as f:
            ss = make_stamps()
            emptyMetadata = PropertyList()
            stamps.writeFits(
                f.name, [ss[0]], emptyMetadata, None, True, True
            )
            with self.assertRaises(RuntimeError):
                stamps.StampsBase.readFits(f.name)

    def testReadFitsReturnsNewClass(self):
        """Test that readFits will return subclass
        """
        class FakeStampsBase(stamps.StampsBase):
            def __init__(self):
                self._metadata = {}
                return

            @classmethod
            def readFitsWithOptions(cls, filename, options):
                return cls()

            def _refresh_metadata(self):
                self._metadata = {}

        fakeStamps = FakeStampsBase.readFitsWithOptions('noFile', {})
        self.assertEqual(type(fakeStamps), FakeStampsBase)


class StampsTestCase(lsst.utils.tests.TestCase):
    """Test Stamps.
    """
    def testAppend(self):
        """Test ability to append to a Stamps object
        """
        ss = make_stamps()
        s = ss[-1]
        ss.append(s)
        self.roundtrip(ss)
        # check if appending something other than a Stamp raises
        with self.assertRaises(ValueError):
            ss.append('hello world')

    def testExtend(self):
        ss = make_stamps()
        ss2 = make_stamps()
        ss.extend([s for s in ss2])
        # check if extending with something other than a Stamps
        # object raises
        with self.assertRaises(ValueError):
            ss.extend(['hello', 'world'])

    def testIO(self):
        """Test the class' write and readFits methods.
        """
        self.roundtrip(make_stamps())

    def testIOone(self):
        """Test the class' write and readFits methods for the special case of
           one stamp.
        """
        self.roundtrip(make_stamps(1))

    def testIOsub(self):
        """Test the class' write and readFits when passing on a bounding box.
        """
        bbox = geom.Box2I(geom.Point2I(3, 9), geom.Extent2I(11, 7))
        ss = make_stamps()
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            options = {'bbox': bbox}
            subStamps = stamps.Stamps.readFitsWithOptions(f.name, options)
            for s1, s2 in zip(ss, subStamps):
                self.assertEqual(bbox.getDimensions(), s2.stamp_im.getDimensions())
                self.assertMaskedImagesAlmostEqual(s1.stamp_im[bbox], s2.stamp_im)

    def testIOarchive(self):
        """Test the class' write and readFits when Stamps contain Persistables.
        """
        self.roundtripWithArchive(make_stamps(use_archive=True))

    def testMetadata(self):
        """Test that metadata is correctly written and read.
        """
        stamps = make_stamps()
        for stamp in stamps:
            self.assertTrue(stamp.metadata is not None)
            self.assertIn('RA_DEG', stamp.metadata)
            self.assertIn('DEC_DEG', stamp.metadata)

    def testStampMetadataTypes(self):
        """Test the types accepted for the metadata of a single stamp.

        A mapping is convenient to write at the call site, but the FITS writer
        needs a `~lsst.daf.base.PropertySet`, so it is normalised on the way in.
        """
        stampIm = afwImage.MaskedImageF(10, 10)
        propertySet = PropertySet()
        propertySet['SRCID'] = 1
        for metadata in ({'SRCID': 1}, PropertyList(), propertySet):
            stamp = stamps.Stamp(stamp_im=stampIm, metadata=metadata)
            self.assertIsInstance(stamp.metadata, PropertyList)
        # A PropertyList is kept as it is, rather than needlessly copied.
        propertyList = PropertyList()
        stamp = stamps.Stamp(stamp_im=stampIm, metadata=propertyList)
        self.assertIs(stamp.metadata, propertyList)
        self.assertIsNone(stamps.Stamp(stamp_im=stampIm).metadata)
        with self.assertRaises(TypeError):
            stamps.Stamp(stamp_im=stampIm, metadata='not a mapping')

    def testStampMetadata(self):
        """Test that metadata held by a single stamp round trips.

        Note that the keys become FITS header keywords, and so are upper-cased.
        """
        nStamps = 3
        ss = make_stamps(nStamps)
        for i, stamp in enumerate(ss):
            stamp.metadata = {'SRCID': 100 + i, 'DETECTOR': i, 'LABEL': f'stamp{i}'}
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, None)
            self.assertEqual(len(ss2), nStamps)
            for i, stamp in enumerate(ss2):
                # Each stamp must get its own values, not a neighbour's. The
                # extensions of the stamps are interleaved, so this fails if
                # the reader indexes the headers by position rather than by
                # EXTVER.
                self.assertEqual(stamp.metadata['SRCID'], 100 + i)
                self.assertEqual(stamp.metadata['DETECTOR'], i)
                self.assertEqual(stamp.metadata['LABEL'], f'stamp{i}')
                # Metadata shared by every stamp is still available.
                self.assertIn('RA_DEG', stamp.metadata)

    def testStampMetadataWithArchive(self):
        """Test the metadata of a single stamp when an Archive is present.

        The archive adds extensions that carry no EXTVER, interleaved with
        those of the stamps.
        """
        nStamps = 3
        ss = make_stamps(nStamps, use_archive=True)
        for i, stamp in enumerate(ss):
            stamp.metadata = {'SRCID': 200 + i}
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, None)
            self.assertEqual([stamp.metadata['SRCID'] for stamp in ss2],
                             [200 + i for i in range(nStamps)])

    def testStampMetadataMultipleValues(self):
        """Test that a key holding several values keeps all of them.

        Copying the metadata key by key would keep only the last value.
        """
        nStamps = 3
        ss = make_stamps(nStamps)
        for i, stamp in enumerate(ss):
            metadata = PropertyList()
            metadata.set('SRCID', 100 + i)
            metadata.set('APFLUX', [1.5 + i, 2.5 + i, 3.5 + i])
            metadata.setComment('SRCID', f'source {i}')
            metadata.add('HISTORY', 'first')
            metadata.add('HISTORY', 'second')
            stamp.metadata = metadata
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, None)
            for i, stamp in enumerate(ss2):
                self.assertEqual(stamp.metadata.getArray('APFLUX'), [1.5 + i, 2.5 + i, 3.5 + i])
                self.assertEqual(stamp.metadata.getArray('HISTORY'), ['first', 'second'])
                self.assertEqual(stamp.metadata.getComment('SRCID'), f'source {i}')

    def testStampMetadataOverridesSharedValue(self):
        """Test a stamp holding a different value to the shared metadata.

        A value that matches the one in the primary header is not repeated in
        the extension, but one that differs is, so that it can override it.
        """
        sharedMetadata = PropertyList()
        sharedMetadata['DETECTOR'] = 1
        sharedMetadata['SURVEY'] = 'wide'
        ss = stamps.Stamps([], metadata=sharedMetadata)
        for i in range(3):
            metadata = PropertyList()
            metadata['SRCID'] = 100 + i
            metadata['SURVEY'] = 'wide'  # same as the shared value
            if i == 1:
                metadata['DETECTOR'] = 99  # differs from the shared value
            ss.append(stamps.Stamp(stamp_im=afwImage.MaskedImageF(10, 10), metadata=metadata))
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            headers = image_hdu_metadata(f.name)
            # Only the differing value is repeated in an extension.
            self.assertEqual([('DETECTOR' in md) for md in headers], [False, True, False])
            self.assertEqual([('SURVEY' in md) for md in headers], [False, False, False])
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, None)
            self.assertEqual([stamp.metadata['DETECTOR'] for stamp in ss2], [1, 99, 1])
            self.assertEqual([stamp.metadata['SURVEY'] for stamp in ss2], ['wide'] * 3)

    def testStampMetadataDoesNotOverrideSharedList(self):
        """Test that a stamp cannot override a shared key holding a list.

        Those lists hold one entry per stamp and are indexed by stamp number,
        so replacing one with a single stamp's stale copy would misalign the
        others.
        """
        ss = make_stamps(3)
        # Every stamp shares the metadata holding the RA_DEG list, so appending
        # leaves each stamp holding a list that is one entry short.
        ss.append(ss[-1])
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            for md in image_hdu_metadata(f.name):
                self.assertNotIn('RA_DEG', md)
                self.assertNotIn('DEC_DEG', md)
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, None)
            self.assertEqual(len(ss2), 4)
            for stamp1, stamp2 in zip(ss, ss2):
                self.assertAlmostEqual(stamp1.position.getRa().asDegrees(),
                                       stamp2.position.getRa().asDegrees())

    def testStampMetadataRewrite(self):
        """Test that rewriting a file does not make its headers grow.

        Reading merges the primary header into the metadata of each stamp, so
        writing it back out must not repeat those keys, nor the ones the writer
        produces for each extension by itself.
        """
        nStamps = 3
        ss = make_stamps(nStamps)
        for i, stamp in enumerate(ss):
            stamp.metadata = {'SRCID': 100 + i}
        with tempfile.NamedTemporaryFile() as f1, tempfile.NamedTemporaryFile() as f2:
            ss.writeFits(f1.name)
            ss2 = stamps.Stamps.readFitsWithOptions(f1.name, None)
            ss2.writeFits(f2.name)
            ss3 = stamps.Stamps.readFitsWithOptions(f2.name, None)
            self.assertEqual([stamp.metadata['SRCID'] for stamp in ss3],
                             [100 + i for i in range(nStamps)])
            for md1, md2 in zip(image_hdu_metadata(f1.name), image_hdu_metadata(f2.name)):
                self.assertEqual(sorted(md1.names()), sorted(md2.names()))
                # A duplicated card shows up as a key holding several values.
                for key in md2.names():
                    self.assertEqual(len(md2.getArray(key)), len(md1.getArray(key)))

    def roundtrip(self, ss):
        """Round trip a Stamps object to disk and check values
        """
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            options = PropertyList()
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, options)
            self.assertEqual(len(ss), len(ss2))
            for s1, s2 in zip(ss, ss2):
                self.assertMaskedImagesAlmostEqual(s1.stamp_im, s2.stamp_im)
                self.assertAlmostEqual(s1.position.getRa().asDegrees(),
                                       s2.position.getRa().asDegrees())
                self.assertAlmostEqual(s1.position.getDec().asDegrees(),
                                       s2.position.getDec().asDegrees())

                for k, v in s1.metadata.items():
                    self.assertIn(k, s2.metadata)
                    self.assertAlmostEqual(v, s2.metadata[k])

    def roundtripWithArchive(self, ss):
        """Round trip a Stamps object, including Archive elements, and check values
        """
        transformTest = TransformTestBaseClass()
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            options = PropertyList()
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, options)
            self.assertEqual(len(ss), len(ss2))
            for s1, s2 in zip(ss, ss2):
                self.assertMaskedImagesAlmostEqual(s1.stamp_im, s2.stamp_im)
                self.assertAlmostEqual(s1.position.getRa().asDegrees(),
                                       s2.position.getRa().asDegrees())
                self.assertAlmostEqual(s1.position.getDec().asDegrees(),
                                       s2.position.getDec().asDegrees())
                transformTest.assertTransformsEqual(s1.archive_element, s2.archive_element)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
