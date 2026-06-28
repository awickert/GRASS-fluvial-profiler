"""
GRASS-inclusive tests for r.fluvial.channelheads.

Need a live GRASS session; run locally (not in CI). Require r.fluvial.channelheads on
the GRASS addon path and rivernetworkx importable. Note this module does NOT
need v.stream.network first (it reads the raw r.stream.extract network without
topology). Run with, e.g.:

    MPLBACKEND=Agg PYTHONPATH=. GRASS_ADDON_BASE=/tmp/rnx-addon \
        grass --tmp-location XY --exec python \
        testsuite/test_rfluvialchannelheads_grass.py
"""

from grass import script as gscript
from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestFluvialHollow(TestCase):
    """DEM -> r.watershed -> r.stream.extract -> r.fluvial.channelheads (no v.stream.network)."""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=200, s=0, e=200, w=0, res=1)
        # a gentle slope (drains off the south edge) plus random roughness gives
        # a real dendritic network (a smooth analytic DEM does not converge to
        # channels). Seeded so the test is reproducible.
        cls.runModule('r.mapcalc', seed=42,
                      expression='dem = 0.4*y() + rand(0.0, 6.0)', overwrite=True)
        cls.runModule('r.watershed', elevation='dem', accumulation='accum',
                      flags='s', overwrite=True)
        cls.runModule('r.stream.extract', elevation='dem', accumulation='accum',
                      stream_vector='streams', direction='draindir',
                      threshold=100, d8cut=0, overwrite=True)
        # deliberately NO v.stream.network: this module needs no topology

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster',
                      name='dem,accum,draindir')
        cls.runModule('g.remove', flags='f', type='vector',
                      name='streams,transition')
        cls.del_temp_region()

    def test_finds_transition_without_topology(self):
        # -c: the synthetic catchment drains off-map, so treat negative
        # accumulation at hollow heads as a boundary artifact. window=0 = none.
        self.assertModule('r.fluvial.channelheads', flags='c', elevation='dem',
                          accumulation='accum', streams='streams',
                          points='transition', window=0, min_slope=1e-9,
                          overwrite=True)
        self.assertVectorExists('transition')
        topo = gscript.vector_info_topo('transition')
        self.assertGreater(int(topo['points']), 0)       # produced points
        cols = gscript.vector_columns('transition')
        self.assertIn('source_cat', cols)                # provenance carried
        self.assertIn('a_star', cols)                    # the fitted break area

    def test_offmap_without_flag_fails(self):
        # without -c, off-map inflow (negative accumulation at heads) must error
        # rather than build a misleading result (issue #9).
        self.assertModuleFail('r.fluvial.channelheads', elevation='dem',
                              accumulation='accum', streams='streams',
                              points='transition_bad', window=0, min_slope=1e-9)


class TestChannelHeadsDrEICH(TestCase):
    """method=dreich (DrEICH): a fastscape-evolved landscape -> chi-z channel heads.
    Uses r.fluvial.fastscape to build coherent terrain (also a cross-module check)."""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=1000, s=0, e=1000, w=0, res=10)
        cls.runModule('r.mapcalc', seed=3, expression='rough = rand(0.0, 2.0)',
                      overwrite=True)
        cls.runModule('r.fluvial.fastscape', input='rough', output='dem',
                      k=3e-5, m=0.45, n=1.0, uplift=3e-3, dt=1500.0, nsteps=250,
                      overwrite=True)

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster', name='rough,dem,dir,netr')
        cls.runModule('g.remove', flags='f', type='vector', name='heads,net')
        cls.del_temp_region()

    def test_finds_dreich_heads(self):
        self.assertModule('r.fluvial.channelheads', method='dreich', elevation='dem',
                          points='heads', threshold=15, window_radius=30,
                          tan_curv_threshold=0.005, min_segment_length=5,
                          overwrite=True)
        self.assertVectorExists('heads')
        topo = gscript.vector_info_topo('heads')
        self.assertGreater(int(topo['points']), 2)        # several channel heads

    def test_emits_fluvial_network(self):
        # also request the downstream fluvial network as a linked directed graph.
        self.assertModule('r.fluvial.channelheads', method='dreich', elevation='dem',
                          points='heads', network='net', threshold=15,
                          window_radius=30, tan_curv_threshold=0.005,
                          min_segment_length=5, overwrite=True)
        self.assertVectorExists('net')
        topo = gscript.vector_info_topo('net')
        self.assertGreater(int(topo['lines']), 2)         # several stream links
        # the network carries v.stream.network topology: a 'tostream' column whose
        # non-zero values all point to real link cats (a converging directed graph).
        cols = gscript.vector_columns('net')
        self.assertIn('tostream', cols)
        recs = gscript.read_command('v.db.select', map='net', columns='cat,tostream',
                                    flags='c').strip().split('\n')
        cats, tos = set(), []
        for line in recs:
            c, t = line.split('|')
            cats.add(int(c)); tos.append(int(t))
        self.assertIn(0, tos)                             # at least one outlet
        for t in tos:
            self.assertTrue(t == 0 or t in cats)         # no dangling downstream link

    def test_external_routing_and_raster_network(self):
        # routing supplied externally (r.watershed -s), plus the raster network.
        self.runModule('r.watershed', flags='s', elevation='dem', drainage='dir',
                       overwrite=True)
        self.assertModule('r.fluvial.channelheads', method='dreich', elevation='dem',
                          direction='dir', network='net', raster_network='netr',
                          threshold=15, window_radius=30, tan_curv_threshold=0.005,
                          min_segment_length=5, overwrite=True)
        # raster stream network is CELL, with link cats matching the vector.
        self.assertRasterExists('netr')
        info = gscript.raster_info('netr')
        self.assertEqual(info['datatype'], 'CELL')
        vcats = set(int(r) for r in gscript.read_command(
            'v.db.select', map='net', columns='cat', flags='c').split())
        self.assertEqual(int(info['max']), max(vcats))    # same cat range

    def test_requires_at_least_one_output(self):
        # method=dreich with no output target must fail, not silently do nothing.
        self.assertModuleFail('r.fluvial.channelheads', method='dreich',
                              elevation='dem', threshold=15, window_radius=30,
                              tan_curv_threshold=0.005, min_segment_length=5)


if __name__ == '__main__':
    test()
