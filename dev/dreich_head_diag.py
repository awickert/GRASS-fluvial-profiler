"""Diagnose the 7 per-junction channel-head mismatches vs _CH."""
import numpy as np
import dev.lsdtt_dreich as L

NR, NC, ND, NDF = L.NR, L.NC, L.ND, L.NDF
z = np.fromfile(f'{L.ALG}/bailey_run_dem_fill.flt', dtype='<f4').reshape(NR, NC)
tcurv = np.fromfile(f'{L.ALG}/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(NR, NC)
fi = L.build_flowinfo(z); L.contributing_area(fi)
src = L.get_sources_index_threshold(fi, 100)
L.junction_network(fi, src)
vj = L.find_valleys(fi, tcurv, src)
L.build_svector(fi); L.distance_from_outlet(fi)

row_of, col_of, recv, JV = fi['row_of'], fi['col_of'], fi['recv'], fi['JunctionVector']
NodeIndex = fi['NodeIndex']

# per-junction profiles
nz = np.where(vj != ND)[0]
order = np.lexsort((col_of[nz], row_of[nz]))
junction_list = vj[nz[order]]
prof = {}
for jn in junction_list:
    dn = int(JV[jn])
    hilltop = L.find_farthest_upslope_node(fi, dn)
    nodeseq, chi, elev = L.build_channel_chi(fi, z, hilltop, dn, 1000.0, 0.525)
    head = L.calculate_channel_heads(nodeseq, chi, elev, 10)
    prof[int(jn)] = dict(dn=dn, hilltop=hilltop, nodeseq=nodeseq, chi=chi, elev=elev, head=head)

ref = np.fromfile(f'{L.ALG}/bailey_run_dem_CH.flt', dtype='<f4').reshape(NR, NC) != NDF
mine = np.load('/tmp/dreich_ch_mine.npy')
ro_r, ro_c = np.where(ref & ~mine)

def junction_of(node):
    cur = node
    for _ in range(100000):
        if vj[cur] != ND:
            return int(vj[cur])
        r = int(recv[cur])
        if r == cur:
            return None
        cur = r
    return None

for k in range(len(ro_r)):
    R = int(NodeIndex[ro_r[k], ro_c[k]])
    jn = junction_of(R)
    p = prof.get(jn)
    print('=== ref-only head node %d at (%d,%d) ; junction %s ===' % (R, ro_r[k], ro_c[k], jn))
    if p is None:
        print('   could not map to a junction'); continue
    ns = p['nodeseq']
    in_prof = R in ns
    myhead = p['head']
    mh_pos = ns.index(myhead) if myhead in ns else -1
    print('   hilltop=%d  profile_len=%d  my_head=%d (pos %d)  R_in_profile=%s'
          % (p['hilltop'], len(ns), myhead, mh_pos, in_prof))
    if in_prof:
        rpos = ns.index(R)
        print('   R pos in profile = %d  (my head pos = %d)  -> SPLIT difference' % (rpos, mh_pos))
        # show test landscape near both
        chi, elev = p['chi'], p['elev']; end = len(chi)
        for hill in sorted(set([mh_pos, rpos, mh_pos-1, mh_pos+1, rpos-1, rpos+1])):
            if 10 <= hill <= end-10:
                r2c, _ = L._reg32(chi[hill:], elev[hill:])
                _, dwh = L._reg32(chi[:hill], elev[:hill])
                test = np.float32(r2c - np.float32((dwh-np.float32(2))/np.float32(2)))
                print('      hill=%d test=%.9f (R2c=%.7f DWh=%.7f)' % (hill, test, r2c, dwh))
    else:
        print('   R NOT on my profile -> different HILLTOP. R drainage area=%d, my hilltop area=%d'
              % (fi['ncontrib'][R], fi['ncontrib'][p['hilltop']]))
