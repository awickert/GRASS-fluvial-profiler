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
#from grass.pygrass.gis import region
from joblib import Parallel, delayed


z = RasterRow('DEM')
z.open('r')

xy = [[529733.000000, 4924574.000000]]

def get_values(coords_list):
    for coords in coords_list:
        print z.get_value(Point(coords[0], coords[1]))
    
def get_value(coords):
    print z.get_value(Point(coords[0], coords[1]))


Parallel(n_jobs=2)(delayed(get_value)(coords for coords in xy*4))
