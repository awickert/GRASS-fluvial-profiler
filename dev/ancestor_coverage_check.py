import numpy as np, networkx as nx
from rivernetworkx.core import build_graph
from rivernetworkx.grass_io import _read_raster, sample_raster
from grass.pygrass.vector import VectorTopo
from grass.script import run_command, region as _region
run_command('g.region', raster='dem', quiet=True)
r=_region(); cn=(float(r['n'])+float(r['s']))/2; ce=(float(r['e'])+float(r['w']))/2; ns=float(r['nsres']); ew=float(r['ewres'])
run_command('g.region', n=cn+600*ns, s=cn-600*ns, e=ce+600*ew, w=ce-600*ew, align='dem', quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_h', drainage='ddir_h', overwrite=True, quiet=True)
run_command('r.stream.extract', elevation='dem', accumulation='acc_h', threshold=133, stream_vector='str_h', stream_raster='srast_h', direction='sdir_h', d8cut=0, overwrite=True, quiet=True)
run_command('r.stream.basins', direction='sdir_h', stream_rast='srast_h', basins='bas_h', overwrite=True, quiet=True)
acc,b=_read_raster('acc_h'); basin_r=_read_raster('bas_h')[0].ravel()
vt=VectorTopo('str_h'); vt.open('r'); recs=[]
for ln in vt.viter('lines'):
    if ln.cat is None: continue
    en=ln.to_array(); x,y=en[:,0],en[:,1]
    keep=np.concatenate(([True],(np.diff(x)!=0)|(np.diff(y)!=0))); x,y=x[keep],y[keep]
    if len(x)<2: continue
    Av=np.abs(sample_raster(acc,x,y,**b))
    if Av[0]>Av[-1]: x,y,Av=x[::-1],y[::-1],Av[::-1]
    recs.append({'cat':int(ln.cat),'x':x,'y':y,'A':Av})
vt.close()
def key(x,y): return (round(float(x),2),round(float(y),2))
upi={key(r['x'][0],r['y'][0]):r['cat'] for r in recs}
for r in recs:
    t=upi.get(key(r['x'][-1],r['y'][-1]),0); r['tostream']=0 if t==r['cat'] else t
rbc={r['cat']:r for r in recs}; G=build_graph(recs)
o=np.argsort(basin_r,kind='stable'); bs=basin_r[o]; uq,st=np.unique(bs,return_index=True); cells_of={}
for i,u in enumerate(uq):
    if np.isfinite(u):
        e=st[i+1] if i+1<len(st) else len(o); cells_of[int(u)]=o[st[i]:e]
print('cat   maxA   nAncestors  cells_of[anc](unmasked)   ratio cells/maxA')
for c in sorted(rbc,key=lambda c: rbc[c]['A'].max(),reverse=True)[:6] + sorted(rbc,key=lambda c: rbc[c]['A'].max())[len(rbc)//2:len(rbc)//2+4]:
    anc=nx.ancestors(G,c)|{c}; tot=sum(len(cells_of[a]) for a in anc if a in cells_of)
    print('  %5d  %7.0f   %6d    %8d            %.2f'%(c, rbc[c]['A'].max(), len(anc), tot, tot/max(1,rbc[c]['A'].max())))
