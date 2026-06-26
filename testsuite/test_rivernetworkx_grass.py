"""
GRASS-inclusive tests for rivernetworkx.grass_io.

These require a live GRASS session and run locally (not in CI), like the
r.richdem gunittest suite. Run with, e.g.:

    python -m grass.gunittest.main --location <test-location> --location-type xy \
        testsuite/test_rivernetworkx_grass.py

The build-network test additionally needs v.fluvial.network on the GRASS path.

rivernetworkx must be importable (e.g. `pip install -e .` in the GRASS env).
"""

import numpy as np
import networkx as nx

from grass.gunittest.case import TestCase
from grass.gunittest.main import test

import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster, sample_raster


class TestSampleRaster(TestCase):
    """sample_raster against a raster with known values (no addon needed)."""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=100, s=0, e=100, w=0, res=1)
        # dem_xy = easting + northing at each cell center
        cls.runModule('r.mapcalc', expression='dem_xy = x() + y()',
                      overwrite=True)

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster', name='dem_xy')
        cls.del_temp_region()

    def test_sample_matches_cell_values(self):
        arr, bounds = _read_raster('dem_xy')
        px = [12.5, 47.5, 88.5]
        py = [88.5, 47.5, 12.5]
        vals = sample_raster(arr, px, py, **bounds)
        # nearest cell center value == its (easting + northing)
        for vx, vy, v in zip(px, py, vals):
            self.assertAlmostEqual(v, vx + vy, delta=1.0)

    def test_out_of_region_is_nan(self):
        arr, bounds = _read_raster('dem_xy')
        vals = sample_raster(arr, [-10.0, 200.0], [50.0, 50.0], **bounds)
        self.assertTrue(all(v != v for v in vals))  # NaN != NaN


class TestBuildNetwork(TestCase):
    """Full read path: DEM -> streams -> v.fluvial.network -> build_network.

    Requires v.fluvial.network on the GRASS path.
    """

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=100, s=0, e=100, w=0, res=1)
        # Parabolic valley converging on x=50, sloping down toward the north edge.
        cls.runModule('r.mapcalc',
                      expression='dem = 0.02*(x()-50)^2 + (100 - y())',
                      overwrite=True)
        cls.runModule('r.watershed', elevation='dem', accumulation='accum',
                      flags='s', overwrite=True)
        cls.runModule('r.stream.extract', elevation='dem', accumulation='accum',
                      stream_vector='streams', direction='draindir',
                      threshold=100, d8cut=0, overwrite=True)
        cls.runModule('v.fluvial.network', map='streams')

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster',
                      name='dem,accum,draindir,accum_holes,dem_holes')
        cls.runModule('g.remove', flags='f', type='vector', name='streams')
        cls.del_temp_region()

    def test_build_network(self):
        # no accumulation: this synthetic catchment drains off-map (negative
        # accumulation), which is exercised separately by the error test below.
        G = rnx.build_network('streams', elevation='dem')
        self.assertGreater(G.number_of_nodes(), 1)
        self.assertIn(0, G.nodes)                       # off-map outlet present
        self.assertTrue(nx.is_directed_acyclic_graph(G))
        for n in G.nodes:
            self.assertIn('s', G.nodes[n])              # cumulative distance set
        self.assertEqual(G.nodes[0]['s'][0], 0)         # outlet is the origin
        self.assertGreater(max(G.nodes[n]['s'][0] for n in G.nodes), 0)

    def test_incomplete_catchment_errors(self):
        # off-map contributing area shows up as negative accumulation; sampling
        # it is currently unsupported and must fail loudly (gscript.fatal).
        with self.assertRaises(SystemExit):
            rnx.build_network('streams', accumulation='accum')

    def test_assume_complete_allows_offmap(self):
        # asserting a complete basin skips the off-map error and recovers
        # positive drainage area (the negative head is treated as a boundary
        # artifact)
        G = rnx.build_network('streams', elevation='dem', accumulation='accum',
                              assume_complete=True)
        self.assertGreater(G.number_of_nodes(), 1)
        A = np.concatenate([np.asarray(d['A'], float)
                            for _, _, d in G.edges(data=True) if 'A' in d])
        A = A[np.isfinite(A)]
        self.assertTrue((A >= 0).all())

    def test_nodata_accumulation_errors(self):
        # NULL accumulation under part of the network is a coverage error, even
        # with assume_complete (the NaN-coverage check is unconditional)
        self.runModule('r.mapcalc',
                       expression='accum_holes = if(y() > 50, null(), accum)',
                       overwrite=True)
        with self.assertRaises(SystemExit):
            rnx.build_network('streams', elevation='dem',
                              accumulation='accum_holes', assume_complete=True)

    def test_nodata_elevation_warns_not_errors(self):
        # NULL elevation under part of the network is a warning, not a fatal
        # (unlike accumulation): it builds, leaving NaN at the missing cells.
        self.runModule('r.mapcalc',
                       expression='dem_holes = if(y() > 50, null(), dem)',
                       overwrite=True)
        G = rnx.build_network('streams', elevation='dem_holes')
        self.assertGreater(G.number_of_nodes(), 1)

    def test_elevation_sampled(self):
        G = rnx.build_network('streams', elevation='dem')
        # at least some edge carries finite sampled elevations
        import numpy as np
        finite = [np.isfinite(np.asarray(d['z'], float)).any()
                  for _, _, d in G.edges(data=True) if 'z' in d]
        self.assertTrue(any(finite))


class TestCompleteCatchment(TestCase):
    """A complete catchment (central peak draining radially outward) keeps
    positive accumulation at the channel heads; negatives appear only on the
    outlet boundary cells, which are tolerated, so building WITH accumulation
    succeeds (unlike the off-map-draining catchment above)."""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=100, s=0, e=100, w=0, res=1)
        cls.runModule('r.mapcalc',
                      expression='demp = 100 - sqrt((x()-50)^2 + (y()-50)^2)',
                      overwrite=True)
        cls.runModule('r.watershed', elevation='demp', accumulation='accp',
                      flags='s', overwrite=True)
        cls.runModule('r.stream.extract', elevation='demp', accumulation='accp',
                      stream_vector='streamsp', direction='ddp',
                      threshold=10, d8cut=0, overwrite=True)
        cls.runModule('v.fluvial.network', map='streamsp')

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster', name='demp,accp,ddp')
        cls.runModule('g.remove', flags='f', type='vector', name='streamsp')
        cls.del_temp_region()

    def test_build_with_accumulation_succeeds(self):
        # negatives only at the outlet boundary -> no head is negative -> builds
        G = rnx.build_network('streamsp', elevation='demp', accumulation='accp')
        self.assertGreater(G.number_of_nodes(), 1)
        self.assertIn(0, G.nodes)
        # drainage area is returned as a positive magnitude (outlet boundary
        # cells, negative in r.watershed, are recovered via abs)
        edge_A = np.concatenate([np.asarray(d['A'], float)
                                 for _, _, d in G.edges(data=True) if 'A' in d])
        edge_A = edge_A[np.isfinite(edge_A)]
        self.assertTrue((edge_A >= 0).all())
        self.assertGreater(edge_A.max(), 0)


if __name__ == '__main__':
    test()
