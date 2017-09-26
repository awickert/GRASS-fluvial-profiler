#!/usr/bin/env python
############################################################################
#
# MODULE:       v.gsflow.reaches
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Builds grid for the MODFLOW component of GSFLOW
#
# COPYRIGHT:    (c) 2016-2017 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  uses inputs from r.stream.extract
 
# More information
# Started December 2016

#%module
#% description: Builds grid for the MODFLOW component of GSFLOW
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#%end

#%option G_OPT_V_INPUT
#%  key: basin
#%  label: Study basin, over which to build a MODFLOW grid
#%  required: yes
#%end

#%option
#%  key: dx
#%  label: Cell size (x / E / zonal), in map units
#%  required: yes
#%end

#%option
#%  key: dy
#%  label: Cell size (y / N / meridional), in map units
#%  required: yes
#%end

#%option G_OPT_V_OUTPUT
#%  key: output
#%  label: MODFLOW grid
#%  required: yes
#%end

##################
# IMPORT MODULES #
##################
# PYTHON
import numpy as np
from matplotlib import pyplot as plt
import sys
import warnings
# GRASS
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import miscellaneous as m
from grass.pygrass.gis import region
from grass.pygrass import vector # Change to "v"?
from grass.script import vector_db_select
from grass.pygrass.vector import Vector, VectorTopo
from grass.pygrass.raster import RasterRow
from grass.pygrass import utils
from grass import script as gscript

###############
# MAIN MODULE #
###############

def main():
    """
    Builds a grid for the MODFLOW component of the USGS hydrologic model,
    GSFLOW.
    """

    options, flags = gscript.parser()
    basin = options['basin']
    dx = options['dx']
    dy = options['dy']
    grid = options['output']
    
    # Create grid
    gscript.use_temp_region()
    g.region(vector=basin, ewres=dx, nsres=dy)
    v.mkgrid(map=grid, overwrite=gscript.overwrite())

    # Cell numbers (row, column, continuous ID)
    v.db_addcolumn(map=grid, columns='id int', quiet=True)
    colNames = np.array(gscript.vector_db_select(grid, layer=1)['columns'])
    colValues = np.array(gscript.vector_db_select(grid, layer=1)['values'].values())
    cats = colValues[:,colNames == 'cat'].astype(int).squeeze()
    rows = colValues[:,colNames == 'row'].astype(int).squeeze()
    cols = colValues[:,colNames == 'col'].astype(int).squeeze()
    nrows = np.max(rows)
    ncols = np.max(cols)
    _id = ncols*(rows-1) + cols
    _id_cat = []
    for i in range(len(_id)):
      _id_cat.append( (_id[i], cats[i]) )
    gridTopo = VectorTopo(grid)
    gridTopo.open('rw')
    cur = gridTopo.table.conn.cursor()
    cur.executemany("update "+grid+" set id=? where cat=?", _id_cat)
    gridTopo.table.conn.commit()
    gridTopo.close()

if __name__ == "__main__":
    main()
