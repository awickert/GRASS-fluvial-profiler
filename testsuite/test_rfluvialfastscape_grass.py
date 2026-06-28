"""
GRASS-inclusive tests for r.fluvial.fastscape.

Need a live GRASS session; run locally (not in CI). Require r.fluvial.fastscape on
the GRASS addon path and rivernetworkx importable. Run with, e.g.:

    MPLBACKEND=Agg PYTHONPATH=. GRASS_ADDON_BASE=/tmp/rnx-addon \
        grass --tmp-location XY --exec python \
        testsuite/test_rfluvialfastscape_grass.py
"""

from grass import script as gscript
from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestFluvialFastscape(TestCase):
    """A noisy low-relief block, uplifted and incised, should build relief and a
    concave channel network (steady-state slope-area)."""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=300, s=0, e=300, w=0, res=10)
        # small initial roughness; edges become base level inside the module
        cls.runModule('r.mapcalc', seed=7, expression='dem0 = rand(0.0, 1.0)',
                      overwrite=True)

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster',
                      name='dem0,dem_evolved', quiet=True)
        cls.del_temp_region()

    def test_builds_relief(self):
        self.assertModule('r.fluvial.fastscape', input='dem0', output='dem_evolved',
                          k=1e-5, m=0.5, n=1.0, uplift=1e-3, dt=2000.0, nsteps=200,
                          overwrite=True)
        self.assertRasterExists('dem_evolved')
        # r.in.bin does not populate range metadata, so query stats directly
        stats = gscript.parse_command('r.univar', map='dem_evolved', flags='g')
        zmax = float(stats['max'])
        # uplift x time vs steady erosion -> finite, bounded relief well above the
        # initial 1 m roughness, and not a blow-up
        self.assertGreater(zmax, 5.0)
        self.assertLess(zmax, 1.0e4)


if __name__ == '__main__':
    test()
