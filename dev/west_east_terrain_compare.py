import numpy as np
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command, region as _region
run_command('g.region', raster='dem', quiet=True)
r=_region(); cn=(float(r['n'])+float(r['s']))/2; ce=(float(r['e'])+float(r['w']))/2; ns=float(r['nsres']); ew=float(r['ewres'])
run_command('g.region', n=cn+600*ns, s=cn-600*ns, e=ce+600*ew, w=ce-600*ew, align='dem', quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_w', format='percent', overwrite=True, quiet=True)
slp,b=_read_raster('slp_w'); dem,_=_read_raster('dem'); acc,_=_read_raster('accum')
nr,nc=slp.shape
S=(slp/100.).ravel(); E=dem.ravel(); A=np.abs(acc).ravel()
cols=np.tile(np.arange(nc),nr); X=b['west']+(cols+0.5)*b['ewres']
fin=np.isfinite(S)&np.isfinite(E)&np.isfinite(A)&(A>0)
XMID=402400.0
print('cut~197; matched elevation band 210-300 m, bluff cells:')
for name,m in [('WEST',(X<XMID)),('EAST',(X>=XMID))]:
    mm=fin&m&(E>197)
    band=mm&(E>=210)&(E<300)
    print('  %s: n=%d medElev=%.0f medSlope(all bluff)=%.3f  | band210-300: n=%d medSlope=%.3f IQR[%.3f,%.3f]'
          %(name, mm.sum(), np.median(E[mm]), np.median(S[mm]), band.sum(),
            np.median(S[band]), np.percentile(S[band],25), np.percentile(S[band],75)))
# binned slope-area, aggregate, west vs east (bluff)
for name,m in [('WEST',(X<XMID)),('EAST',(X>=XMID))]:
    mm=fin&m&(E>197)&(S>1e-4)
    la=np.log10(A[mm]); ls=np.log10(S[mm]); ed=np.linspace(la.min(),la.max(),12); idx=np.digitize(la,ed)
    s=''.join(' (%.1f,%+.2f)'%(np.median(la[idx==i]),np.median(ls[idx==i])) for i in range(1,12) if (idx==i).sum()>50)
    print('  %s binned SA:%s'%(name,s))
