#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.profiler
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Build long profiles and slope--accumulation (e.g., 
#               slope--area) diagrams of a river network
#
# COPYRIGHT:    (c) 2016-2018 Andrew Wickert
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

def moving_average(x, y, window):
    """
    Create a moving average every <window/2> points with an averaging
    distance of <window>, but including the the first point + window/2
    and the last point - window/2
    (so distance to last point could be irregular)
    """
    x = np.array(x)
    y = np.array(y)
    out_x = np.arange(x[0]+window/2., x[-1]-window/2., window)
    out_x = np.hstack((out_x, x[-1]-window/2.))
    out_y = []
    for _x in out_x:
        out_y.append( np.mean(y[ (x < _x + window/2.) * 
                                 (x > _x - window/2.) ]))
    return out_x, out_y

def get_xEN(cat, x0=0., streams=streams):
    """
    """
    data = vector.VectorTopo(streams) # Create a VectorTopo object
    data.open('r') # Open this object for reading

    coords = data.cat(cat_id=cat, vtype='lines')[0]

    x_downstream = [x0]
    E = [ coords[0].x ]
    N = [ coords[0].y ]
    for _i in range(1, len(coords)):
        x_downstream.append(coords[_i-1].distance(coords[_i]))
        E.append(coords[_i].x)
        N.append(coords[_i].y)
    x_downstream = np.cumsum(x_downstream)[::-1]
    
    data.close()
    
    return x_downstream, E, N

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
                                   reg.ewres + reg.get_bbox().west \
                                   - reg.ewres
            self.ymin = np.floor( (self.ymin - reg.get_bbox().south ) / 
                                   reg.nsres ) * \
                                   reg.nsres + reg.get_bbox().south \
                                   - reg.nsres
            self.xmax = np.ceil( (self.xmax - reg.get_bbox().east ) / 
                                  reg.ewres ) * \
                                  reg.ewres + reg.get_bbox().east \
                                   + reg.ewres
            self.ymax = np.ceil( (self.ymax - reg.get_bbox().north ) / 
                                  reg.nsres ) * \
                                  reg.nsres + reg.get_bbox().north \
                                   + reg.nsres
        self.bbox = Bbox()
        self.bbox.north = self.ymax
        self.bbox.south = self.ymin
        self.bbox.west = self.xmin
        self.bbox.east = self.xmax


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
    colNames = np.array(vector_db_select(streams)['columns'])
    colValues = np.array(vector_db_select(streams)['values'].values())
    warnings.warn('tostream is not generalized')
    tostream = colValues[:,colNames == 'tostream'].astype(int).squeeze()
    cats = colValues[:,colNames == 'cat'].astype(int).squeeze() # = "fromstream"

    # We can loop over this list to get the shape of the full river network.
    selected_cats = []
    segment = _cat
    selected_cats.append(segment)

    # Get all cats in network
    data = vector.VectorTopo(streams) # Create a VectorTopo object
    data.open('r') # Open this object for reading

    if direction == 'downstream':
        gscript.message("Extracting drainage pathway...",)
        # Get network
        while selected_cats[-1] != 0:
            selected_cats.append(int(tostream[cats == selected_cats[-1]]))
        #x.append(selected_cats[-1])
        selected_cats = selected_cats[:-1] # remove 0 at end
        gscript.message("Done.")
        
        
    elif direction == 'upstream':
        gscript.message("Extracting drainage network...",)
        # GENERALIZE COLUMN NAME!!!!!!!!
        tostream_col = np.where(np.array(data.table.columns.names())
                                == 'tostream')[0][0]
        terminalCats = [_cat]
        terminal_x_values = [0]
        netcats = []
        net_tocats = []
        while len(terminalCats) > 0:
            for cat in terminalCats:
                netcats.append(cat)
                # ALSO UNADVISABLE NAME -- NEED TO GET TOSTREAM, GENERALIZED
                #print data.table_to_dict()
                colnum = np.where( np.array(data.table.columns.names()) 
                                   == 'tostream')[0][0]
                net_tocats.append(data.table_to_dict()[cat][colnum])
            oldcats = terminalCats
            terminalCats = []
            for cat in oldcats:
                terminalCats += list(cats[tostream == cat])
        #data.close()
        netcats = np.array(netcats)
        net_tocats = np.array(net_tocats)
        
        selected_cats = netcats
        gscript.message("Done.")
        
    segments = []
    for cat in selected_cats:
        points_with_cat = data.cat(cat_id=cat, vtype='lines')[0]
        subcoords = []
        for point in points_with_cat:
            subcoords.append([point.x, point.y])
        segments.append( rn.Segment(_id=cat, to_ids=tostream[cats == cat]) )
        segments[-1].set_EastingNorthing(ENarray=subcoords)
        segments[-1].calc_x_from_EastingNorthing()
        # x grid spacing
        #print segments[-1].Easting[-1], segments[-1].Northing[-1]
        #print segments[-1].EastingNorthing[-1]
        #print ""
        if dx_target is not None:
            dx_target = float(dx_target)
            segments[-1].set_target_dx_downstream(dx_target)
            segments[-1].densify_x_E_N()
    data.close()
    
    net = rn.Network(segments)

    bbox = BoundingBox(points_xy=net.segments_xy_flattened())
    reg_to_revert = region.Region()
    reg = region.Region() # to limit region for computational efficiency
    reg.set_bbox(bbox.bbox)
    reg.write()
    
    # Network extraction
    if outstream:
        selected_cats_str = list(np.array(selected_cats).astype(str))
        selected_cats_csv = ','.join(selected_cats_str)
        v.extract( input=streams, output=outstream, \
                   cats=selected_cats_csv, overwrite=overwrite_flag )
    
    
    # All coordinates
    coords = net.segments_xy_flattened()
    #x_downstream = 
    
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
  
    
    # TEMPORARY!!!!
    #x_downstream = get_xEN()
    #x_downstream_0 = x_downstream[0]

    # Analysis

    # Downstream distances -- 0 at mouth
    net.compute_x_in_network()

    # Elevation
    if elevation:
        gscript.message("Elevation")
        _include_z = True
        # Load DEM
        griddata = garray.array()
        griddata.read(elevation)
        griddata = np.flipud(griddata)
        # Interpolate: nearest or linear?
        x = np.arange(reg.west + reg.ewres/2., reg.east, reg.ewres)
        y = np.arange(reg.south + reg.nsres/2., reg.north, reg.nsres)
        itp = RegularGridInterpolator( (x, y), griddata.transpose(), 
                                       method='nearest')
        _i = 0
        _lasti = 0
        _nexti = 0
        for segment in net.segment_list:
            try:
                segment.set_z( itp(segment.EastingNorthing) )
            except:
                print segment.EastingNorthing
                print np.vstack((segment.Easting_original, segment.Northing_original)).transpose()
                sys.exit()
            if _i > _nexti:
                gscript.core.percent( _i, len(net.segment_list), np.floor(_i - _lasti))
                _nexti = float(_nexti) + len(net.segment_list)/10.
                if _nexti > len(net.segment_list):
                    _nexti = len(net.segment_list) - 1
            _lasti = _i
            _i += 1
        gscript.core.percent(1, 1, 1)
        del griddata
        #warnings.warn('Need to handle window in network')
        #gscript.core.percent(1, 1, 1)
    else:
        _include_z = False

    # Slope
    if slope:
        gscript.message("Slope")
        _include_S = True
        _slope = RasterRow(slope)
        _slope.open('r')
        _i = 0
        _lasti = 0
        _nexti = 0
        for segment in net.segment_list:
            sen = segment.EastingNorthing # all E,N
            S = []
            for row in sen:
                #try:
                S.append(_slope.get_value(Point(row[0], row[1])))
                #except:
                #    print "ERROR"
                if _i > _nexti:
                    gscript.core.percent(_i, len(coords), np.floor(_i - _lasti))
                    _nexti = float(_nexti) + len(coords)/10.
                    if _nexti > len(coords):
                        _nexti = len(coords) - 1
                _lasti = _i
                _i += 1
            # MAKE SETTER FOR THIS!!!!
            segment.channel_slope = np.array(S)
        if window is not None:
            pass
            #net.smooth_window()
            #_x_downstream, _S = moving_average(x_downstream_0, S, window)
        _slope.close()
        S = np.array(S)
        S_0 = S.copy()
        gscript.core.percent(1, 1, 1)
    else:
        _include_S = False

    # Accumulation / drainage area
    if accumulation:
        gscript.message("Accumulation")
        _include_A = True
        accumulation = RasterRow(accumulation)
        accumulation.open('r')
        _i = 0
        _lasti = 0
        _nexti = 0
        for segment in net.segment_list:
            A = []
            sen = segment.EastingNorthing # all E,N
            for row in sen:
                A.append(accumulation.get_value(Point(row[0], row[1])) 
                                                          * accum_mult)
                if _i > _nexti:
                    gscript.core.percent(_i, len(coords), np.floor(_i - _lasti))
                    _nexti = float(_nexti) + len(coords)/10.
                    if _nexti > len(coords):
                        _nexti = len(coords) - 1
                _lasti = _i
                _i += 1
            # MAKE SETTER FOR THIS!!!!
            segment.channel_flow_accumulation = np.array(A)
        accumulation.close()
        A = np.array(A)
        A_0 = A.copy()
        """
        if window is not None:
            _x_downstream, A = moving_average(x_downstream_0, A, window)
        """
        gscript.core.percent(1, 1, 1)
    else:
        _include_A = False

    # Revert to original region
    reg_to_revert

    # Smoothing
    if window is not None:
        net.smooth_window(window)

    # Plotting
    if 'LongProfile' in plots:
        plt.figure()
        if window:
            for segment in net.segment_list:
                plt.plot(segment.x/1000., segment.z_smoothed, 'k-', linewidth=2)
        else:
            for segment in net.segment_list:
                plt.plot(segment.x/1000., segment.z, 'k-', linewidth=2)
        #plt.plot(x_downstream/1000., z, 'k-', linewidth=2)
        plt.xlabel('Distance from mouth [km]', fontsize=16)
        plt.ylabel('Elevation [m]', fontsize=20)
        plt.tight_layout()
    if 'SlopeAccum' in plots:
        plt.figure()
        if window:
            for segment in net.segment_list:
                _y_points = segment.channel_slope_smoothed[
                                segment.channel_flow_accumulation_smoothed > 0
                                ]
                _x_points = segment.channel_flow_accumulation_smoothed[
                                segment.channel_flow_accumulation_smoothed > 0
                                ]
        else:
            for segment in net.segment_list:
                _y_points = segment.channel_slope[
                                    segment.channel_flow_accumulation > 0
                                    ]
                _x_points = segment.channel_flow_accumulation[
                                    segment.channel_flow_accumulation > 0
                                    ]
        plt.loglog(_x_points, _y_points, 'k.', alpha=.5)
        plt.xlabel(accum_label, fontsize=20)
        plt.ylabel('Slope [$-$]', fontsize=20)
        plt.tight_layout()
    if 'SlopeDistance' in plots:
        plt.figure()
        for segment in net.segment_list:
            plt.plot(segment.x/1000., segment.channel_slope, 'k-', linewidth=2)
        plt.xlabel('Distance downstream [km]', fontsize=16)
        plt.ylabel('Slope [$-$]', fontsize=20)
        plt.tight_layout()
    if 'AccumDistance' in plots:
        plt.figure()
        for segment in net.segment_list:
            _x_points = segment.x[segment.channel_flow_accumulation > 0]
            _y_points = segment.channel_flow_accumulation[
                                         segment.channel_flow_accumulation > 0
                                         ]
            plt.plot(_x_points/1000., _y_points, 'k.', alpha=.5)
        plt.xlabel('Distance downstream [km]', fontsize=16)
        plt.ylabel(accum_label, fontsize=20)
        plt.tight_layout()
    plt.show()
    
    # Saving data -- will need to update for more complex data structures!
    if outfile:
        net.compute_profile_from_starting_segment()
        _outfile = np.vstack((net.long_profile_header, net.long_profile_output))
        np.savetxt(outfile, _outfile, '%s')
    else:
        pass
        
    #print net.accum_from_headwaters[1] - net.slope_from_headwaters[1]

    """
    for segment in net.segment_list:
        print segment.channel_flow_accumulation_smoothed
        print segment.channel_slope_smoothed
        print segment.channel_flow_accumulation_smoothed - \
              segment.channel_slope_smoothed
    """
    
    """
    if options['outfile_original'] is not '':
        header = ['x_downstream', 'E', 'N']
        outfile = np.hstack((np.expand_dims(x_downstream_0, axis=1), coords))
        if _include_S:
            header.append('slope')
            outfile = np.hstack((outfile, np.expand_dims(S_0, axis=1)))
        if _include_A:
            if (options['units'] == 'm2') or (options['units'] == 'km2'):
                header.append('drainage_area_'+options['units'])
            elif (options['units'] == 'cumecs') or (options['units'] == 'cfs'):
                header.append('water_discharge_'+options['units'])
            else:
                header.append('flow_accumulation_arbitrary_units')
            outfile = np.hstack((outfile, np.expand_dims(A_0, axis=1)))
        header = np.array(header)
        outfile = np.vstack((header, outfile))
        np.savetxt(options['outfile_original'], outfile, '%s')
    if options['outfile_smoothed'] is not '':
        header = ['x_downstream', 'E', 'N']
        # E, N on smoothed grid
        x_downstream, E = moving_average(x_downstream_0, coords[:,0], window)
        x_downstream, N = moving_average(x_downstream_0, coords[:,1], window)
        # Back to output
        outfile = np.hstack((np.expand_dims(x_downstream, axis=1),
                             np.expand_dims(E, axis=1),
                             np.expand_dims(N, axis=1)))
        if _include_S:
            header.append('slope')
            outfile = np.hstack((outfile, np.expand_dims(S, axis=1)))
        if _include_A:
            if (options['units'] == 'm2') or (options['units'] == 'km2'):
                header.append('drainage_area_'+options['units'])
            elif (options['units'] == 'cumecs') or (options['units'] == 'cfs'):
                header.append('water_discharge_'+options['units'])
            else:
                header.append('flow_accumulation_arbitrary_units')
            outfile = np.hstack((outfile, np.expand_dims(A, axis=1)))
        header = np.array(header)
        outfile = np.vstack((header, outfile))
        np.savetxt(options['outfile_smoothed'], outfile, '%s')
    """
    
    #print net.segment_list[0].x net.segment_list[0].E, net.segment_list[0].N
    
if __name__ == "__main__":
    main()

