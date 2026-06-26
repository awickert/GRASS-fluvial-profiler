"""
GRASS-inclusive tests for r.stream.hollow.

Need a live GRASS session; run locally (not in CI). Require r.stream.hollow on
the GRASS addon path and rivernetworkx importable. Note this module does NOT
need v.stream.network first (it reads the raw r.stream.extract network without
topology). Run with, e.g.:

    MPLBACKEND=Agg PYTHONPATH=. GRASS_ADDON_BASE=/tmp/rnx-addon \
        grass --tmp-location XY --exec python \
        testsuite/test_rstreamhollow_grass.py
"""

from grass import script as gscript
from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestStreamHollow(TestCase):
    """DEM -> r.watershed -> r.stream.extract -> r.stream.hollow (no v.stream.network)."""

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
        self.assertModule('r.stream.hollow', flags='c', elevation='dem',
                          accumulation='accum', streams='streams',
                          output='transition', window=0, min_slope=1e-9,
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
        self.assertModuleFail('r.stream.hollow', elevation='dem',
                              accumulation='accum', streams='streams',
                              output='transition_bad', window=0, min_slope=1e-9)


if __name__ == '__main__':
    test()
