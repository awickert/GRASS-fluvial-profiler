"""
GRASS-inclusive tests for v.fluvial.profiler (rebuilt on rivernetworkx).

Like the rivernetworkx GRASS suite, these need a live GRASS session and run
locally (not in CI). They require v.stream.network AND v.fluvial.profiler on the
GRASS addon path, and rivernetworkx importable. Run with, e.g.:

    MPLBACKEND=Agg PYTHONPATH=. GRASS_ADDON_BASE=/tmp/rnx-addon \
        grass --tmp-location XY --exec python \
        testsuite/test_vfluvialprofiler_grass.py
"""

import json
import os
import tempfile

from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestProfiler(TestCase):
    """Full path: DEM -> streams -> v.stream.network -> v.fluvial.profiler."""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=100, s=0, e=100, w=0, res=1)
        cls.runModule('r.mapcalc',
                      expression='dem = 0.02*(x()-50)^2 + (100 - y())',
                      overwrite=True)
        cls.runModule('r.watershed', elevation='dem', accumulation='accum',
                      flags='s', overwrite=True)
        cls.runModule('r.stream.extract', elevation='dem', accumulation='accum',
                      stream_vector='streams', direction='draindir',
                      threshold=100, d8cut=0, overwrite=True)
        cls.runModule('v.stream.network', map='streams')

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster',
                      name='dem,accum,draindir')
        cls.runModule('g.remove', flags='f', type='vector', name='streams,raw')
        cls.del_temp_region()

    def test_downstream_profile_and_json(self):
        tmp = tempfile.mkdtemp()
        txt = os.path.join(tmp, 'prof.txt')
        js = os.path.join(tmp, 'prof.json')
        # no accumulation: the synthetic catchment drains off-map (negative
        # accumulation), covered by test_incomplete_catchment_errors.
        self.assertModule('v.fluvial.profiler', cat=1, streams='streams',
                          elevation='dem', dx_target=5, outfile=txt, json=js)
        with open(txt) as f:
            lines = f.read().splitlines()
        self.assertGreater(len(lines), 2)            # header + data
        self.assertIn('s', lines[0].split())         # labelled columns
        d = json.load(open(js))
        ids = [n['id'] for n in d['nodes']]
        self.assertIn(0, ids)                         # off-map outlet present

    def test_incomplete_catchment_errors(self):
        # sampling accumulation on this off-map-draining catchment must fail
        # (negative flow accumulation; issue #9)
        self.assertModuleFail('v.fluvial.profiler', cat=1, streams='streams',
                              elevation='dem', accumulation='accum')

    def test_requires_tostream(self):
        # A raw r.stream.extract vector has no tostream column; the module must
        # fail helpfully rather than proceed.
        self.runModule('r.stream.extract', elevation='dem', accumulation='accum',
                       stream_vector='raw', direction='draindir_raw',
                       threshold=100, d8cut=0, overwrite=True)
        self.assertModuleFail('v.fluvial.profiler', cat=1, streams='raw',
                              elevation='dem')


if __name__ == '__main__':
    test()
