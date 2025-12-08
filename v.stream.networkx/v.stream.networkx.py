#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.profiler
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
# Started 14 October 2016
#%module
#% description: Build a linked stream network: each link knows its downstream link
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#%end
#%option
#%  key: cat
#%  label: Starting line segment category
#%  required: yes
#%  guidependency: layer,column
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
#%option
#%  key: direction
#%  type: string
#%  label: Which directon to march: up or down
#%  options: upstream,downstream
#%  answer: downstream
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
#%option G_OPT_R_INPUT
#%  key: slope
#%  label: Map of slope created by r.slope.area
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
#%  key: dx_target
#%  type: double
#%  label: Target distance between output stream points [map units]
#%  required: no
#%end
#%option
#%  key: window
#%  label: Averaging distance [map units]
#%  required: no
#%end
#%option
#%  key: plots
#%  type: string
#%  label: Plots to generate
#%  options: LongProfile,SlopeAccum,SlopeDistance,AccumDistance
#%  required: no
#%  multiple: yes
#%end
#%option
#%  key: outfile
#%  type: string
#%  label: output file
#%  required: no
#%end

##################
# IMPORT MODULES #
##################
# CUSTOM
# temporary patch
import sys
sys.path.insert(0, '/home/awickert/dataanalysis/GRASS-fluvial-profiler/v.stream.profiler/')
import RiverNetwork as rn
# PYTHON
import numpy as np
import pandas as pd
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
slope = options['slope']
if slope == '': slope = None    
accumulation = options['accumulation']
if accumulation == '': accumulation = None
direction = options['direction']
if direction == '': direction = None
streams = options['streams']
if streams == '': streams = None
outstream = options['outstream']
if outstream == '': outstream = None
try:
    window = float(options['window'])
except:
    window = None
try:
    dx_target = float(options['dx_target'])
except:
    dx_target = None
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

"""
# Test where it isn't working -- just needed to be re-run?
_cat = 149
print(_cat)
coords = vt.cat(cat_id=_cat, vtype='lines')[0]
rast_list = []
for _c in coords:
    x, y = _c.x, _c.y  # map coordinates
    with RasterRow(rastname, mode="r") as rast:
        _reg = region.Region()
        rast_list.append(rast.get_value((x, y), _reg))
rast_lol.append(rast_list)
"""

# Move node elevations from segments to nodes

def drop_downstream_edge_array_values(G, attr_names):
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

"""
# Test
parent = 1258
child = 1259
data = G.get_edge_data(1258, 1259)

# Test again

parent = 149
child = 150

#df_edges.loc[ df_edges['cat'] == 149, 'z' ] = values_from_raster( [parent], 'dem' )[0]
#df_edges.loc[parent-1]['A'] = values_from_raster( [parent], 'accumulation' )

df_edges.at[ parent-1, 'z' ] = values_from_raster( [parent], 'dem' )[0]
df_edges.at[ parent-1, 'A' ] = values_from_raster( [parent], 'accumulation' )[0]

data = G.get_edge_data(parent, child)

# Check out upstream node
G.nodes[parent]
# Check out link (edge): empty as it should be!
G.edges[parent, child]

attr = 'x'
arr = data[attr]

first = arr[0]

# ---- move to parent node ----
# store as a list on the parent node (accumulate values)
node_attr = G.nodes[parent].get(attr, [])
# copy if it's not a list yet
if not isinstance(node_attr, list):
    node_attr = [node_attr]
node_attr.append(first)
G.nodes[parent][attr] = node_attr
"""


### ALL THIS STUFF PURE COPY/PASTE
colNames = vector_db_select(streams)['columns']
colValues = vector_db_select(streams)['values'].values()
df_edges = pd.DataFrame( data=colValues, columns=colNames )
cats = list( df_edges['cat'].astype(int) ) # = "fromstream"

vt = vector.VectorTopo(streams) # Create a VectorTopo object
vt.open('r') # Open this object for reading

_x = []
_y = []
_su = []
_sd = []
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
### END

# Get all attributes
df_edges['z'] = values_from_raster( cats, 'dem' )
df_edges['A'] = values_from_raster( cats, 'accumulation' )

df_edges['cat'] = df_edges['cat'].astype(int)
df_edges['tostream'] = df_edges['tostream'].astype(int)


# from setupDomain.py (modified)
import networkx as nx

# Generate network structure with data on edges
G = nx.from_pandas_edgelist(df_edges, source='cat', target='tostream', edge_key='cat', edge_attr=True, create_using=nx.DiGraph)

# Move data from edges to nodes
#G.add_nodes_from((n, dict(d)) for n, d in df_nodes.iterrows())
attr_names = ['x', 'y', 's_upstream', 's_downstream', 'z', 'A']
pull_first_from_edges_to_parents(G, attr_names)
drop_downstream_edge_array_values(G, attr_names)

# Drop x1 and y1 and x2 and y2 later !!!!!!!!!!!!!!! ######################


"""
# Now, to see if it worked: try to plot!

# To do this, assume we will have an acyclic binary tree and walk up it
# We start by defining the number of the mouth node; this could later
# set externally
mouth_node = 1258
s = nx.bfs_tree(G, mouth_node, reverse=True)
# We start by placing distance = 0 at the mouth node
# !!!!!!!!!!!!! THIS IS WHERE I CREATED THE INITIAL PROBLEM !!!!!!!!!!!
G.nodes[mouth_node]['s'] = 0

# Find upstream and downstream nodes and links
children = list(G.successors(mouth_node)) # Downstream
parents = list(G.predecessors(mouth_node)) # Upstream
"""

"""
# This is where our 1258 problem was!
# And our offsets
mouth_node = 0
# Loop to get all upstream distances
for node in s.nodes:
    # First, update the distance upstream from the node
    if node == mouth_node:
        # I never like having a single-case "if". seems like a big waste.
        # Ways around this?
        G.nodes[node]['s_upstream'] = 0. # MAKE THIS A LIST OF 1 ITEM
        # Could also check if out_edges is nonexistent...
    else:
        # Otherwise, look at the next link downstream
        edges = G.out_edges(node)
        if len(edges) != 1:
            # Declare an error if we don't have exactly 1 downstream edge
            # This is built for only directed, convergent, acyclic graphs
            sys.exit() # error
        else:
            edge = next(iter(edges))
        # First item = upstream
        # Try first on edge
        if len( G.edges[edge]['x'] ) > 0:
            _dx = G.nodes[node]['x'] - G.edges[edge]['x'][0]
            _dy = G.nodes[node]['y'] - G.edges[edge]['y'][0]
            ds = ( _dx**2 + _dy**2 )**.5
            # Update the node
            G.nodes[node]['s_upstream'] = ds + G.edges[edge]['s_upstream'][0]
        # If no items on edge, then it just connects two nodes with no cell
        # between them (i.e., tributary junctions on subsequent cells)
        else:
            _dx = G.nodes[node]['x'] - G.nodes[edge[-1]]['x'][0]
            _dy = G.nodes[node]['y'] - G.nodes[edge[-1]]['y'][0]
            ds = ( _dx**2 + _dy**2 )**.5
            # Update the node
            G.nodes[node]['s_upstream'] = ds + G.nodes[edge[-1]]['s_upstream'][0]
    # Then walk upstream to the edges touching this node and iterate
    # over them.
    # For a tributary network, this should be 1 for the mouth node, 2
    # for all middle nodes, and 0 for the upstream-most-end nodes
    edges_to_node = G.in_edges(node)
    for edge in edges_to_node:
        # Find distances upstream from next downstream node
        # I have set the convention that x,y,z are arranged
        # upstream to downstream
        # MODIFICATION: Nodes are already 1-element lists
        _x = np.hstack(( G.edges[edge]['x'], G.nodes[node]['x'] )) 
        _y = np.hstack(( G.edges[edge]['y'], G.nodes[node]['y'] )) 
        # Continue without assuming that we've already calculated dx, dy, ds
        ds = ( np.diff(_x)**2 + np.diff(_y)**2 )**.5
        # Increasing (positive) distance with distance from the river mouth 
        G.edges[edge]['s_upstream']= np.cumsum(ds[::-1])[::-1] + G.nodes[node]['s_upstream']
"""
    


# Define all values as nan or 0 for offmap
G.nodes[0]['x'] = [np.nan]
G.nodes[0]['y'] = [np.nan]
G.nodes[0]['s_upstream'] = [0]
G.nodes[0]['s_downstream'] = [np.nan]
G.nodes[0]['z'] = [np.nan]
G.nodes[0]['A'] = [np.nan]
G.nodes[0]['s'] = [0] # Total distance upstream of outlet

# Overall downstream distance
def bfs_upward(G, start):
    R = G.reverse(copy=False)  # just a view, no data duplication
    for node in nx.bfs_tree(R, start):
        yield node
"""
for n in bfs_upward(G, 0):
    print(n)
    

for i in range(len(edges)):
        edge = next(iter( edges ))
        print(edge)
"""

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

"""
# It's 1258! Probably just not a list because I was using it as an example. 
# ABOVE: Because I had it as the mouth_node in the now-commented section.
# This code section is shorter and better.
        try:
            G.nodes[parent]['s'] = [G.nodes[child]['s'][0] + G.nodes[parent]['s_upstream'][0]]
        except:
            # Find and fix the reason this is a float later
            # Probably because of these two-unit tribs
            # THIS WILL CAUSE PROBLEMS: FIX FIX FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            G.nodes[parent]['s'] = [np.nan]#[G.nodes[child]['s'][0] + G.nodes[parent]['s_upstream']]
"""

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







# Try plotting
# Color each of the segments, up to and including the nodes above them.
# Use gray for the downstream connectors

# Because of this plan, we can start by simply looping over nodes and seeing
# what is downstream of each of them
s_on_segment_list = []
var_on_segment_list = []
varname = 'z' # <-- should be able to alter how we obtain this
for node in G.nodes:
    # By this point, just assume binary tree
    edges = G.out_edges(node)
    if len(edges) == 0:
        break # we are at the mouth -- this should happen only once.
    elif len(edges) > 1:
        sys.exit() # something bad has happened
    else:
        edge = next(iter( edges ))
    _s = np.hstack(( G.nodes[node]['s_upstream'], G.edges[edge]['s_upstream'] ))
    _var = np.hstack(( G.nodes[node][varname], G.edges[edge][varname] ))
    s_on_segment_list.append(_s)
    var_on_segment_list.append(_var)

# Then we can go through nodes and find what is upstream of them
s_on_node_list = []
var_on_node_list = []
varname = 'z' # <-- should be able to alter how we obtain this
for node in G.nodes:
    # Let's see what's coming in
    edges = G.in_edges(node)
    if len(edges) == 0:
        continue # upstream-most nodes; these are already included in segments
    else:
        # Here, should have multiple connections
        for edge in edges:
            if len( G.edges[edge]['s_upstream'] ) > 0:
                _s = np.hstack(( [G.edges[edge]['s_upstream'][-1]], G.nodes[node]['s_upstream'] ))
            else:
                _s = np.hstack(( G.nodes[edge[0]]['s_upstream'], G.nodes[node]['s_upstream'] ))                
            if len( G.edges[edge][varname] ) > 0:
                _var = np.hstack(( [G.edges[edge][varname][-1]],
                                    G.nodes[node][varname] ))
            else:
                _var = np.hstack(( G.nodes[edge[0]][varname],
                                    G.nodes[node][varname] ))
            s_on_node_list.append(_s)
            var_on_node_list.append(_var)

#plt.ion()
plt.figure()
for i in range(len(s_on_node_list)):
    plt.plot( s_on_node_list[i], var_on_node_list[i], 'k-', linewidth=3,
                alpha=0.3 )
for i in range(len(s_on_segment_list)):
    plt.plot( s_on_segment_list[i], var_on_segment_list[i], 'k-',
                linewidth=3, alpha=1 )

plt.xlabel('Upstream distance')
plt.ylabel('Elevation') # hard-coded for now

plt.show()


for i in range(len(dfe)):
    print( np.min(dfe.loc[i]['s_upstream']) )

for i in range(len(s_on_segment_list)):
    print( np.min(s_on_segment_list[i]) )

for node in G.nodes:
    print( np.min(G.edges[edge]['s_upstream']) )

# Add overall upstream distance to the attributes
s_on_segment_list = []
var_on_segment_list = []
varname = 'z' # <-- should be able to alter how we obtain this
for node in G.nodes:
    # By this point, just assume binary tree
    edges = G.out_edges(node)
    if len(edges) == 0:
        continue # we are at the mouth
    elif len(edges) > 1:
        sys.exit() # something bad has happened
    else:
        edge = next(iter( edges ))
    _s = np.hstack(( G.nodes[node]['s_upstream'], G.edges[edge]['s_upstream'] ))
    _var = np.hstack(( G.nodes[node][varname], G.edges[edge][varname] ))
    s_on_segment_list.append(_s)
    var_on_segment_list.append(_var)


# Do this for one sub-watershed

parents_recursive = nx.ancestors(G, 1248)

plt.figure()
for i in range(len(parents_recursive)):
    plt.plot( parents_recursive[i], parents_recursive[i], 'k-', linewidth=3,
                alpha=0.3 )
for i in range(len(s_on_segment_list)):
    plt.plot( s_on_segment_list[i], var_on_segment_list[i], 'k-',
                linewidth=3, alpha=1 )

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

    # Parsing inside function
    _cat = int(options['cat'])
    overwrite_flag = gscript.overwrite()
    elevation = options['elevation']
    if elevation == '': elevation = None    
    slope = options['slope']
    if slope == '': slope = None    
    accumulation = options['accumulation']
    if accumulation == '': accumulation = None
    direction = options['direction']
    if direction == '': direction = None
    streams = options['streams']
    if streams == '': streams = None
    outstream = options['outstream']
    if outstream == '': outstream = None
    outfile = options['outfile']
    if outfile == '': outfile = None
    # !!!!!!!!!!!!!!!!!
    # ADD SWITCHES TO INDIVIDUALLY SMOOTH SLOPE, ACCUM, ETC.
    # !!!!!!!!!!!!!!!!!
    try:
        window = float(options['window'])
    except:
        window = None
    try:
        dx_target = float(options['dx_target'])
    except:
        dx_target = None
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
    colNames = vector_db_select(streams)['columns']
    colValues = vector_db_select(streams)['values'].values()
    warnings.warn('tostream is not generalized')
    dfnet = pd.DataFrame( data=colValues, columns=colNames )
    tostream = dfnet['tostream'].astype(int)
    cats = dfnet['cat'].astype(int) # = "fromstream"

    # Get all cats in network
    data = vector.VectorTopo(streams) # Create a VectorTopo object
    data.open('r') # Open this object for reading

    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # UPDATE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """
    ##### FIND RIGHT SPOT TO ADD CLASS STUFF HERE/BELOW ####
    
    # Extract x points in network
    data = vector.VectorTopo(streams) # Create a VectorTopo object
    data.open('r') # Open this object for reading
    
    coords = []
    _i = 0
    for i in range(len(data)):
        if type(data.read(i+1)) is vector.geometry.Line:
            if data.read(i+1).cat in selected_cats:
                coords.append(data.read(i+1).to_array())
                gscript.core.percent(_i, len(selected_cats), 100./len(selected_cats))
                _i += 1
    gscript.core.percent(1, 1, 1)
    coords = np.vstack(np.array(coords))
    
    _dx = np.diff(coords[:,0])
    _dy = np.diff(coords[:,1])
    x_downstream_0 = np.hstack((0, np.cumsum((_dx**2 + _dy**2)**.5)))
    x_downstream = x_downstream_0.copy()
    
    data.close()
    """

if __name__ == "__main__":
    colNames = vector_db_select(streams)['columns']
    colValues = vector_db_select(streams)['values'].values()
    dfnet = pd.DataFrame( data=colValues, columns=colNames )

    main()

