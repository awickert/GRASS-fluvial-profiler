options = {}
options['window'] = 0
options['units'] = 'km2'
options['streams'] = 'streams_all'
options['direction'] = 'upstream'
options['accum_mult'] = 1
options['plots'] = ''
#options['cat'] = 113
options['cat'] = 103
options['elevation'] = 'DEM'


##################
# IMPORT MODULES #
##################
# CUSTOM
import RiverNetwork as rn
# PYTHON
import numpy as np
from matplotlib import pyplot as plt
import sys
# GRASS
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.gis import region
from grass.pygrass import vector # Change to "v"?
from grass.script import vector_db_select
from grass.pygrass.vector import Vector, VectorTopo
from grass.pygrass.raster import RasterRow
from grass.pygrass import utils
from grass import script as gscript
from grass.pygrass.vector.geometry import Point
import warnings
from multiprocessing import Pool
from grass.script import array as garray
from scipy.interpolate import RegularGridInterpolator

# Parsing
window = float(options['window'])
accum_mult = float(options['accum_mult'])
if options['units'] == 'm2':
    accum_label = 'Drainage area [m$^2$]'
elif options['units'] == 'km2':
    accum_label = 'Drainage area [km$^2$]'
elif options['units'] == 'cumecs':
    accum_label = 'Water discharge [m$^3$ s$^{-1}$]'
elif options['units'] == 'cfs':
    accum_label = 'Water discharge [cfs]'
else:
    accum_label = 'Flow accumulation [$-$]'
plots = options['plots'].split(',')

# Attributes of streams
colNames = np.array(vector_db_select(options['streams'])['columns'])
colValues = np.array(vector_db_select(options['streams'])['values'].values())
warnings.warn('tostream is not generalized')
tostream = colValues[:,colNames == 'tostream'].astype(int).squeeze()
cats = colValues[:,colNames == 'cat'].astype(int).squeeze() # = "fromstream"

# We can loop over this list to get the shape of the full river network.
selected_cats = []
segment = int(options['cat'])
selected_cats.append(segment)
x = []
z = []

# ParallelTest

# Import data by segment



# NECESSARY

from grass.pygrass.gis import region
from grass.pygrass.vector.basic import Bbox 

class BoundingBox(object):
    """
    Easily define a bounding box around your data source, padded to include
    the raster grid cells (if these are important)
    """
    def __init__(self, points_xy=None, align_to_region=True, xmin=None, 
                 xmax=None, ymin=None, ymax=None):
        if points_xy is not None:
            points = np.array(points_xy)
            self.xmin = np.min(points[:,0])
            self.xmax = np.max(points[:,0])
            self.ymin = np.min(points[:,1])
            self.ymax = np.max(points[:,1])
        else:
            self.xmin = xmin
            self.ymin = ymin
            self.xmax = xmax
            self.ymax = ymax
        if align_to_region is not None:
            reg = region.Region()
            self.xmin = np.floor( (self.xmin - reg.get_bbox().west) / 
                                   reg.ewres ) * \
                                   reg.ewres + reg.get_bbox().west
            self.ymin = np.floor( (self.ymin - reg.get_bbox().south ) / 
                                   reg.nsres ) * \
                                   reg.nsres + reg.get_bbox().south
            self.xmax = np.ceil( (self.xmax - reg.get_bbox().east ) / 
                                  reg.ewres ) * \
                                  reg.ewres + reg.get_bbox().east
            self.ymax = np.ceil( (self.ymax - reg.get_bbox().north ) / 
                                  reg.nsres ) * \
                                  reg.nsres + reg.get_bbox().north
        self.bbox = Bbox()
        self.bbox.north = self.ymax
        self.bbox.south = self.ymin
        self.bbox.west = self.xmin
        self.bbox.east = self.xmax



# START TESTS
"""
# DOWNSTREAM
selected_cats = []
segment = int(options['cat'])
selected_cats.append(segment)
while selected_cats[-1] != 0:
    selected_cats.append(int(tostream[cats == selected_cats[-1]]))
if selected_cats[-1] == 0:
    selected_cats = selected_cats[:-1] # remove 0 at end if flow is offmap
"""

# Extract x points in network
data = vector.VectorTopo(options['streams']) # Create a VectorTopo object
data.open('r') # Open this object for reading
segments = []
for cat in selected_cats:
    points_with_cat = data.cat(cat_id=cat, vtype='lines')[0]
    subcoords = []
    for point in points_with_cat:
        subcoords.append([point.x, point.y])
    segments.append( rn.Segment(_id=cat, to_ids=tostream[cats == cat]) )
    segments[-1].set_EastingNorthing(ENarray=subcoords)
    segments[-1].calc_x_from_EastingNorthing()
data.close()

net = rn.Network(segments)


bbox = BoundingBox(points_xy=net.segments_xy_flattened())
reg = region.Region()
reg.set_bbox(bbox.bbox)
reg.write()

DEM = garray.array()
DEM.read(options['elevation'])
DEM = np.flipud(DEM)

# nearest or linear?
x = np.arange(reg.west + reg.ewres/2., reg.east, reg.ewres)
y = np.arange(reg.south + reg.nsres/2., reg.north, reg.nsres)
itp = RegularGridInterpolator( (x, y), DEM.transpose(), method='nearest')

for segment in net.segment_list:
    segment.set_z( itp(segment.EastingNorthing) )
    
net.compute_x_in_network()

# END TESTS


for segment in net.segment_list:
    plt.plot(segment.x, segment.z)

#itp = RegularGridInterpolator( (y, x), DEM, method='nearest')
#res = itp((coords[0][:,1], coords[0][:,0]))



# DOES NOT PARALLELIZE IN A FRIENDLY WAY

def readInSegment(raster, seg):
    pass
    
def readPoint(raster, point):
    return raster.get_value(point)

def readPoint(xy):
    DEM = RasterRow(options['elevation'])
    DEM.open('r')
    point = Point(xy[0], xy[1])
    DEM.close()
    return DEM.get_value(point)


coordsList = []
for point in coords[-1]:
    coordsList.append([point.x, point.y])

p = Pool(3)
p.map(readPoint, coordsList)


DEM = RasterRow(options['elevation'])
DEM.open('r')
z = []
for __i in range(len(netcats)):
    cat = netcats[__i]
    zsub = []
    for _j in range(len(E[__i])):
        zsub.append(DEM.get_value(Point(E[__i][_j], N[__i][_j])))
    z.append(zsub)
    
