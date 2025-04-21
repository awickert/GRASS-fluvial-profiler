#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.network
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Attach IDs of upstream and downstream nodes as well as the
#               category value of the next downstream stream segment
#               (0 if the stream exits the map)
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
# Started 14 October 2016

#%module
#% description: Build a linked stream network: each link knows its downstream link
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#%end

#%option G_OPT_V_INPUT
#%  key: map
#%  label: Vector stream network from r.stream.extract
#%  required: yes
#%  guidependency: layer,column
#%end

#%option
#%  key: upstream_easting_column
#%  type: string
#%  description: Upstream easting (or x or lon) column name
#%  answer: x1
#%  required : no
#%end

#%option
#%  key: upstream_northing_column
#%  type: string
#%  description: Upstream northing (or y or lat) column name
#%  answer: y1
#%  required : no
#%end

#%option
#%  key: downstream_easting_column
#%  type: string
#%  description: Downstream easting (or x or lon) column name
#%  answer: x2
#%  required : no
#%end

#%option
#%  key: downstream_northing_column
#%  type: string
#%  description: Downstream northing (or y or lat) column name
#%  answer: y2
#%  required : no
#%end

#%option
#%  key: tostream_cat_column
#%  type: string
#%  description: Adjacent downstream segment cat (0 = offmap flow)
#%  answer: tostream
#%  required : no
#%end

##################
# IMPORT MODULES #
##################
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

###############
# MAIN MODULE #
###############

def main():
    """
    Links each river segment to the next downstream segment in a tributary 
    network by referencing its category (cat) number in a new column. "0"
    means that the river exits the map.
    """

    options, flags = gscript.parser()
    streams = options['map']
    x1 = options['upstream_easting_column']
    y1 = options['upstream_northing_column']
    x2 = options['downstream_easting_column']
    y2 = options['downstream_northing_column']

    streamsTopo = VectorTopo(streams)
    #streamsTopo.build()

    # 1. Get vectorTopo
    streamsTopo.open(mode='rw')
    """
    points_in_streams = []
    cat_of_line_segment = []

    # 2. Get coordinates
    for row in streamsTopo:
        cat_of_line_segment.append(row.cat)
        if type(row) == vector.geometry.Line:
            points_in_streams.append(row)
    """
    
    # 3. Coordinates of points: 1 = start, 2 = end
    try:
        streamsTopo.table.columns.add(x1,'double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add(y1,'double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add(x2,'double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add(y2,'double precision')
    except:
        pass
    try:
        streamsTopo.table.columns.add('tostream','int')
    except:
        pass
    streamsTopo.table.conn.commit()

    # Is this faster than v.to.db?
    """
    cur = streamsTopo.table.conn.cursor()
    for i in range(len(points_in_streams)):
        cur.execute("update streams set x1="+str(points_in_streams[i][0].x)+" where cat="+str(cat_of_line_segment[i]))
        cur.execute("update streams set y1="+str(points_in_streams[i][0].y)+" where cat="+str(cat_of_line_segment[i]))
        cur.execute("update streams set x2="+str(points_in_streams[i][-1].x)+" where cat="+str(cat_of_line_segment[i]))
        cur.execute("update streams set y2="+str(points_in_streams[i][-1].y)+" where cat="+str(cat_of_line_segment[i]))
    streamsTopo.table.conn.commit()
    streamsTopo.build()
    """
    # v.to.db Works more consistently, at least
    # NO MORE BUILDING COLUMNS FIRST???
    streamsTopo.close()
    v.to_db(map=streams, option='start', columns=x1+','+y1, overwrite=True)
    v.to_db(map=streams, option='end', columns=x2+','+y2, overwrite=True)

    # 4. Read in and save the start and end coordinate points
    colNames = vector_db_select(streams)['columns']
    colValues = vector_db_select(streams)['values'].values()
    dfnet = pd.DataFrame( data=colValues, columns=colNames )

    # 5. Build river network
    x1y1_all = dfnet[[x1,y1]]
    for i in range( len(dfnet) ):
        fr_idx = i # Index is the one that flow comes from
        row = dfnet.loc[i]
        fr_cat = row['cat']
        catx = dfnet[row[x2] == dfnet[x1]].cat
        caty = dfnet[row[y2] == dfnet[y1]].cat
        cats = list(set(catx).intersection(set(caty)))
        if ( len(cats) == 0):
            dfnet.loc[fr_idx, 'to_cat'] = -1
            continue
        elif len(cats) > 1:
            print("Diverging graph!")
            # Should probably exit condition here
            print(i)
            print(row)
            print(catx)
            print(caty)
            break
        else:
            # Must be just one cat
            cat = cats[0]
            to_idx = catx[catx == cat].index[0]
            to_cat = dfnet.loc[to_idx]['cat']
            dfnet.loc[fr_idx, 'to_cat'] = to_cat
            dfnet.loc[to_idx, 'fr_cat'] = fr_cat

    # This gives us a set of downstream-facing adjacencies.
    # We will update the database with it.
    streamsTopo.build()
    streamsTopo.open('rw')
    cur = streamsTopo.table.conn.cursor()
    # Default to 0 if no stream flows to it
    cur.execute("update "+streams+" set tostream=0")
    for i in range(len(dfnet)):
        cat = dfnet.loc[i, 'cat']
        to_cat = dfnet.loc[i, 'to_cat']
        # Might be able to chunk it all straight in, but keeping the old
        # way, which isn't too slow.
        cur.execute("update "+streams+" set tostream="+str(to_cat)+" where cat="+cat)
    streamsTopo.table.conn.commit()
    #streamsTopo.build()
    streamsTopo.close()

    gscript.message('Drainage topology built. Check "tostream" column for the downstream cat.')
    gscript.message('A cat value of 0 indicates the downstream-most segment.')

if __name__ == "__main__":
    main()

