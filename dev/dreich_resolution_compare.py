"""Compare DrEICH channel heads across resolutions against the 1 m "truth"
(Mid Bailey Run). Run INSIDE the MidBaileyRun/dreich_restest GRASS session:

    PROJ_DATA=/usr/share/proj PYTHONPATH=<repo> \
      grass /home/awickert/grassdata/MidBaileyRun/dreich_restest \
        --exec python3 dev/dreich_resolution_compare.py [std|adapt]

Truth = heads_std_1m. For each resolution it reports head count, and head-recovery
vs truth: recall (fraction of 1 m heads with a detected head within tol),
precision (fraction of detected heads within tol of a 1 m head), and the median
nearest-neighbour distance each way. tol = max(50 m, 2*res) by default.
"""
import sys
import numpy as np
from grass.script import read_command

RES = [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 30]


def coords(vmap):
    """(N,2) easting/northing array of a vector points map, or None if absent."""
    try:
        txt = read_command('v.out.ascii', input=vmap, format='point',
                           separator='comma', quiet=True).strip()
    except Exception:
        return None
    if not txt:
        return np.empty((0, 2))
    pts = []
    for line in txt.split('\n'):
        f = line.split(',')
        pts.append((float(f[0]), float(f[1])))
    return np.array(pts)


def nn_dist(a, b):
    """For each point in a, Euclidean distance to nearest point in b."""
    if len(a) == 0 or len(b) == 0:
        return np.full(len(a), np.inf)
    d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
    return d.min(1)


def main():
    scheme = sys.argv[1] if len(sys.argv) > 1 else 'adapt'
    truth = coords('heads_std_1m')
    nt = len(truth)
    print('truth (1 m): %d heads. Radius for coincidence = 30 m (Grieve et al. 2016).' % nt)
    print('  (b) count_kept = N/Ntruth ; (a) fit = median distance of a detected '
          'head to nearest 1 m head')
    print('res  heads  count_kept  fit_med(m)  sensitivity  reliability')
    for N in RES:
        vmap = 'heads_std_1m' if N == 1 else 'heads_%s_%dm' % (scheme, N)
        det = coords(vmap)
        if det is None:
            print('%3d    --   (map %s absent)' % (N, vmap))
            continue
        tol = 30.0                                    # paper's fixed coincidence radius
        d_t = nn_dist(truth, det)                     # each truth head -> nearest detected
        d_d = nn_dist(det, truth)                     # each detected head -> nearest truth
        sens = float(np.mean(d_t <= tol)) if nt else float('nan')   # = recall
        rel = float(np.mean(d_d <= tol)) if len(det) else float('nan')  # = precision
        print('%3d  %5d    %6.3f     %7.1f      %.2f         %.2f'
              % (N, len(det), len(det) / nt, np.median(d_d) if len(det) else -1, sens, rel))


if __name__ == '__main__':
    main()
