#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.networkx
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Convert a river network into a NetworkX object
#
# COPYRIGHT:    (c) 2025 Andrew Wickert
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
# Started 14 October 2016, 07 December 2025
#%module
#% description: Create a NetworkX river-network object with cumulative downstream distances
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#%end
#%option G_OPT_V_INPUT
#%  key: streams
#%  label: Vector input of stream network created by r.stream.extract
#%  required: yes
#%end
#%option G_OPT_V_OUTPUT
#%  key: outstream
#%  label: Vector output stream
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: elevation
#%  label: Topography (DEM)
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: accumulation
#%  label: Flow accumulation raster
#%  required: no
#%end
#%option
#%  key: units
#%  type: string
#%  label: Flow accumulation units
#%  options: m2, km2, cumecs, cfs
#%  required: no
#%end
#%option
#%  key: accum_mult
#%  type: double
#%  label: Multiplier to convert flow accumulation to your chosen unit
#%  answer: 1
#%  required: no
#%end
#%option
#%  key: outfile
#%  type: string
#%  label: Output file for NetworkX / Pandas
#%  required: no
#%end

##################
# IMPORT MODULES #
##################
# CUSTOM
# temporary patch
import sys
sys.path.insert(0, '/home/awickert/dataanalysis/GRASS-fluvial-profiler/v.stream.profiler/')
#import RiverNetwork as rn # Now using NetworkX instead of my custom network code
# PYTHON
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import networkx as nx
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
from grass.script import core as gcore
from grass.pygrass.vector.geometry import Point
import warnings
from grass.script import array as garray
from grass.raster import read_raster
from scipy.interpolate import RegularGridInterpolator
from grass.pygrass.gis import region
from grass.pygrass.vector.basic import Bbox

###################
# PARSER - GLOBAL #
###################

options, flags = gscript.parser()

_cat = int(options['cat'])
overwrite_flag = gscript.overwrite()
elevation = options['elevation']
if elevation == '': elevation = None    
accumulation = options['accumulation']
if accumulation == '': accumulation = None
streams = options['streams']
if streams == '': streams = None
outstream = options['outstream']
if outstream == '': outstream = None
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

###################
# UTILITY MODULES #
###################

def values_from_raster( cats, rastname ):
    vt = vector.VectorTopo(streams) # Create a VectorTopo object
    vt.open('r') # Open this object for reading
    rast_lol = []
    for _cat in cats:
        print(_cat)
        coords = vt.cat(cat_id=_cat, vtype='lines')[0]
        rast_list = []
        for _c in coords:
            x, y = _c.x, _c.y  # map coordinates
            with RasterRow(rastname, mode="r") as rast:
                _reg = region.Region()
                rast_list.append(rast.get_value((x, y), _reg))
        rast_lol.append(rast_list)
    vt.close()
    return rast_lol

# Move node elevations from segments to nodes

def drop_downstream_edge_array_values(G, attr_names):
    """
    Remove downstream-most entries in arrays, which duplicate those that are
    or will be found on the nodes at tributary junctions.
    """
    for parent, child, data in G.edges(data=True):
        if child == 0:
            # Leave downstream-most elevations on segments leaving the map.
            # For valid reasons, these might not all have the same elevation.
            continue
        for attr in attr_names:
            if attr not in data:
                continue

            arr = data[attr]

            # Handle Python list or numpy array
            if arr is None or len(arr) == 0:
                continue

            # ---- remove from edge array ----
            # If it's a list, mutate in place:
            if isinstance(arr, list):
                arr.pop(-1)
            else:
                # e.g. numpy array: create a sliced copy
                data[attr] = arr[:-1]


def pull_first_from_edges_to_parents(G, attr_names):
    """
    Assign the upstream-most entry in arrays of attributes held on the edges
    to the upstream parent nodes, and then remove it from the edges.
    
    Therefore, nodes will share information at tributary junctions.  
    
    For each directed edge parent -> child and each attribute in attr_names:
    - take the first element of the edge's list attribute
    - append it to a list attribute on the parent node
    - remove that first element from the edge list
    """
    for parent, child, data in G.edges(data=True):
        for attr in attr_names:
            if attr not in data:
                continue

            arr = data[attr]

            # Handle Python list or numpy array
            if arr is None or len(arr) == 0:
                continue

            first = arr[0]

            # ---- move to parent node ----
            # store as a list on the parent node (accumulate values)
            node_attr = G.nodes[parent].get(attr, [])
            # copy if it's not a list yet
            if not isinstance(node_attr, list):
                node_attr = [node_attr]
            node_attr.append(first)
            G.nodes[parent][attr] = node_attr

            # ---- remove from edge array ----
            # If it's a list, mutate in place:
            if isinstance(arr, list):
                arr.pop(0)
            else:
                # e.g. numpy array: create a sliced copy
                data[attr] = arr[1:]

def bfs_upward(G, start):
    """
    Breadth-first search going upwards.
    Use this to update overall distance upstream from all outlets
    in a single sweep through the network, relying on those network components
    closer to the outlet to be updated first
    """
    R = G.reverse(copy=False)  # just a view, no data duplication
    for node in nx.bfs_tree(R, start):
        yield node



#########################
# PLOTTING. MOVE LATER. #
#########################

#plt.ion()
plt.figure()
for n in bfs_upward(G, 0):
    edges = G.in_edges(n)
    for parent, child in edges:
        plt.plot( G.edges[parent,child]['s'], G.edges[parent,child]['z'], 'k-', linewidth=3, alpha=1 )
        plt.plot( G.nodes[parent]['s'], G.nodes[parent]['z'], 'ko', alpha=1 )
plt.xlabel('Upstream distance')
plt.ylabel('Elevation') # hard-coded for now
plt.show()




# Now try with a farther upstream starting point

# Iterate in BFS through all: test and print
for n in bfs_upward(G, 1258):
    print(n)
    edges = G.in_edges(n, data=True)
    for parent, child, data in edges:
        print(parent, child)

#plt.ion()
plt.figure()
for n in bfs_upward(G, 1258):
    edges = G.in_edges(n)
    for parent, child in edges:
        plt.plot( G.edges[parent,child]['s'], G.edges[parent,child]['z'], 'k-', linewidth=3, alpha=1 )
        plt.plot( G.nodes[parent]['s'], G.nodes[parent]['z'], 'ko', alpha=1 )
plt.xlabel('Upstream distance')
plt.ylabel('Elevation') # hard-coded for now
plt.show()

        
###############
# MAIN MODULE #
###############

def main():
    """
    Links each river segment to the next downstream segment in a tributary 
    network by referencing its category (cat) number in a new column. "0"
    means that the river exits the map.
    """

    # Attributes of streams
    colNames = vector_db_select(streams)['columns']
    colValues = vector_db_select(streams)['values'].values()
    df_edges = pd.DataFrame( data=colValues, columns=colNames )
    cats = list( df_edges['cat'].astype(int) ) # = "fromstream"
    
    # Vector topology
    vt = vector.VectorTopo(streams) # Create a VectorTopo object
    vt.open('r') # Open this object for reading

    # Add distance and position information
    _x = [] # easting
    _y = [] # northing
    _su = [] # upstream-directed along-stream distance
    _sd = [] # downstream-directed along-stream distance
    for _cat in cats:
        print(_cat)
        # Extract and calculate E, N, and along-stream distance
        coords = vt.cat(cat_id=_cat, vtype='lines')[0]
        EN = coords.to_array()
        _diffs = np.diff(EN, axis=0)
        ds_downstream = ( (_diffs**2).sum(axis=1) )**.5
        s_downstream = np.concatenate( [[0], np.cumsum(ds_downstream)] )
        s_upstream = s_downstream[-1] - s_downstream
        # Insert results into DataFrame
        _x.append(EN[:,0])
        _y.append(EN[:,1])
        _su.append(s_upstream)
        _sd.append(s_downstream)

    vt.close()

    df_edges['s_upstream'] = _su
    df_edges['s_downstream'] = _sd
    df_edges['x'] = _x
    df_edges['y'] = _y

    # Get all attributes
    df_edges['z'] = values_from_raster( cats, 'dem' )
    df_edges['A'] = values_from_raster( cats, 'accumulation' )

    df_edges['cat'] = df_edges['cat'].astype(int)
    df_edges['tostream'] = df_edges['tostream'].astype(int)


# from setupDomain.py (modified)


# Generate network structure with data on edges
G = nx.from_pandas_edgelist(df_edges, source='cat', target='tostream', edge_key='cat', edge_attr=True, create_using=nx.DiGraph)

# Move data from edges to nodes
#G.add_nodes_from((n, dict(d)) for n, d in df_nodes.iterrows())
attr_names = ['x', 'y', 's_upstream', 's_downstream', 'z', 'A']

gcore.message("Moving upstream-most data points from streams to nodes above.")
pull_first_from_edges_to_parents(G, attr_names)
gcore.message("Processing started…")
drop_downstream_edge_array_values(G, attr_names)

# Drop x1 and y1 and x2 and y2 later !!!!!!!!!!!!!!! ######################


# Define all values as nan or 0 for offmap
G.nodes[0]['x'] = [np.nan]
G.nodes[0]['y'] = [np.nan]
G.nodes[0]['s_upstream'] = [0]
G.nodes[0]['s_downstream'] = [np.nan]
G.nodes[0]['z'] = [np.nan]
G.nodes[0]['A'] = [np.nan]
G.nodes[0]['s'] = [0] # Total distance upstream of outlet

# Overall downstream distance

# Iterate in BFS through all: test and print
for n in bfs_upward(G, 0):
    print(n)
    edges = G.in_edges(n, data=True)
    for parent, child, data in edges:
        print(parent, child)

# Iterate in BFS through all, and update values.
# "s" will just be total distance upstream of outlet.
for n in bfs_upward(G, 0):
    print(n)
    edges = G.in_edges(n)
    for parent, child in edges:
        print(parent, child)
        # Update node
        G.nodes[parent]['s'] = [G.nodes[child]['s'][0] + G.nodes[parent]['s_upstream'][0]]
        # Update edge
        G.edges[parent,child]['s'] = G.nodes[child]['s'][0] + G.edges[parent,child]['s_upstream']





    coords = []
    _i = 0
    for i in range(len(data)):
        if type(data.read(i+1)) is vector.geometry.Line:
            if data.read(i+1).cat in selected_cats:
                coords.append(data.read(i+1).to_array())
                gcore.percent(_i, len(selected_cats), 100./len(selected_cats))
                _i += 1
    gcore.percent(1, 1, 1)
    coords = np.vstack(np.array(coords))
    
    _dx = np.diff(coords[:,0])
    _dy = np.diff(coords[:,1])
    x_downstream_0 = np.hstack((0, np.cumsum((_dx**2 + _dy**2)**.5)))
    x_downstream = x_downstream_0.copy()

if __name__ == "__main__":
    colNames = vector_db_select(streams)['columns']
    colValues = vector_db_select(streams)['values'].values()
    dfnet = pd.DataFrame( data=colValues, columns=colNames )

    main()

