#! /usr/local/env python

# 06 May 2018
# ADW

import numpy as np
from scipy.interpolate import interp1d
import warnings

class Segment(object):

    def __init__(self, _id, to_ids=None, segment_type='stream', 
                 direction='to_downstream'):
        """
        A segment of a stream, between tributary junctions;
        ID of the stream should be an integer >0
        """
        self.id = int(_id)
        self.to_ids = np.array([to_ids]).astype(int)
        self.x = None
        self.z = None
        self.XP = None
        self.B = None
        #self.E = None
        #self.N = None
        self.target_dx_downstream = None
        self.segment_type = segment_type
        self.allowed_segment_types = ['stream',
                                      'lake']
        self.direction = direction
        self.from_ids = []
        
        if self.id == 0:
            warnings.warn('ID of 0 by default indicates offmap flow')
        
    ###########
    # SETTERS #
    ###########

    def set_x(self, x):
        """
        Set distances downstream
        x can be overwritten with the network
        x_local stays local
        """
        self.x = np.array(x)
        self.x_local = self.x.copy()
        
    def set_z(self, z):
        """
        Set elevation
        """
        self.z = np.array(z)
        
    def set_XP(self, XP):
        """
        Set cross-sectional elevation profile; from this, the at-channel
        valley width can always be determined.
        """
        self.XP = XP
        
    def set_B(self, B):
        """
        Set valley width at the elevation of the channel.
        Just a scalar value at each downstream position
        """
        self.B = B
        
    def set_B_from_XP(self):
        """
        Use the elevation and cross-sectional profile to set the width
        of the channel
        """
        # Something with interpolation
        # Look at TerraPIN for how this may work
        pass
        
    def set_type(self, segment_type):
        if segment_type in self.allowed_segment_types:
            self.segment_type = segment_type
        else:
            warnings.warn('Segment type not valid: please choose one of'+
                         ', '.join(a))
        
    def set_Q(self, Q):
        """
        Set discharge with distance downstream
        """
        self.Q = Q

    def set_A(self, A):
        """
        Set drainage area with distance downstream
        """
        self.A = A
        
    def set_target_dx_downstream(self, target_dx_downstream):
        """
        Set spacing between points in a downstream direction
        This will almost certainly not divide evenly into the length of a
        stream segment, and will therefore be used as the target distance 
        between nodes, which will be adjusted for an integer number of 
        divisible distances using np.ceil
        """
        self.target_dx_downstream = target_dx_downstream
        
    def set_fromids(self, from_id):
        """
        Set the ID of the stream from which 
        """
        self.from_ids = from_id

    def append_fromids(self, from_id):
        """
        Set the ID of the stream from which 
        """
        self.from_ids.append(from_id)
        
    #def set_E(self, E):
    #    """
    #    Set easting of points
    #    """
    #    self.E = E
    
    #def set_N(self, N):
    #    """
    #    Set northing of points
    #    """
    #    self.N = N
        
    def set_EastingNorthing(self, ENarray=None, Easting=None, Northing=None):
        """
        Set Northing and Easting at once
        """
        if ENarray is not None:
            self.EastingNorthing = np.array(ENarray)
            self.Easting = self.EastingNorthing[:,0]
            self.Northing = self.EastingNorthing[:,1]
        elif (Easting is not None) and (Northing is not None):
            self.Easting = Easting
            self.Northing = Northing
        else:
            warnings.warn("Improper inputs")            
        
    def calc_x_from_EastingNorthing(self):
        """
        Calculate downstream distances from Easting and Northing, relative
        only to the individual segment
        """
        if (self.Easting is not None) and (self.Northing is not None):
            dx = ( np.diff(self.Easting)**2 + np.diff(self.Northing)**2 )**.5
            self.x = np.hstack((0., np.cumsum(dx)))
            self.x_local = self.x.copy()
        else:
            warnings.warn("Both Easting and Northing must be defined")
            
    def densify_x_E_N(self):
        """
        Especially if E, N, are taken only at vertices, there can be a 
        signficant lack of resolution in the network. This seeks to densify 
        the network based on a target downstream dx.
        
        This will not preserve points at the coordinates of all of the corners
        of the fluvial network. It will create backups of the original grid
        to do this. It could be reprogrammed to optionally preserve these 
        corners.
        
        Run this after calc_x_from_EastingNorthing and set_target_dx_downstream
        """
        if self.target_dx_downstream is None:
            warnings.warn("Nothing will be done; run set_target_dx_downstream")
        elif (self.Easting is None) or (self.Northing is None):
            warnings.warn("Nothing will be done; Set Easting and Northing")
        elif self.x_local is None:
            warnings.warn("Nothing will be done; Set x (and x_local)")
        else:
            # Backup original
            self.Easting_original = self.Easting.copy()
            self.Northing_original = self.Northing.copy()
            self.x_local_original = self.x_local.copy()
            ixE = interp1d(self.x_local_original, self.Easting_original)
            ixN = interp1d(self.x_local_original, self.Northing_original)
            #print self.Easting[-1], self.Northing[-1]
            # +1 because it includes the bookends
            # Ceil instead of round for short segments
            nx = np.ceil( ( np.max(self.x_local) - np.min(self.x_local) )
                            / self.target_dx_downstream ) + 1.
            self.x_local = np.linspace( np.min(self.x_local_original), 
                                        np.max(self.x_local_original), 
                                        nx )
            self.Easting = ixE(self.x_local)
            self.Northing = ixN(self.x_local)
            self.EastingNorthing = \
                 np.vstack((self.Easting, self.Northing)).transpose()
                    
class BoundingBox(object):
    """
    Easily define a bounding box around your data source, padded to include
    the raster grid cells (if these are important)
    """
    def __init__(self, points_xy=None):
        points = np.array(points_xy)
        self.west = np.min(points[:,0])
        self.east = np.max(points[:,0])
        self.south = np.min(points[:,1])
        self.north = np.max(points[:,1])

class Network(object):

    def __init__(self, segment_list, downstream_x=0., offmap_id=0):
        """
        Tools to build a river network
        Downstream_x = the x-value or values for the offmap cell(s);
        they default to 0 and will be expanded to a list in this case
        offmap_id is the ID given to offmap segments (to_id = offmap_id if a 
        segment flows offmap). This defaults to 0 following GRASS GIS 
        convention.
        """
        self.segment_list = np.array(segment_list)
        # Convenient info to have in lists
        # identifiers and connections
        self.ids = []
        self.to_ids = []
        for segment in self.segment_list:
            self.ids.append(segment.id)
            self.to_ids.append(segment.to_ids)
        self.ids = np.array(self.ids)
        #self.to_ids = np.array(self.to_ids)
        self.offmap_id = offmap_id # 0 default following GRASS GIS
        # Get ID of downstream-most segment
        # Is downstream-most if "tosegment" is not in the list
        for segment in self.segment_list:
            if np.sum(self.ids == segment.to_ids[0][0]) == 0:
                self.downstream_most_id = segment.id
                break
        self.downstream_x = downstream_x
        
    def compute_fromstreams(self):
        """
        Compute IDs of flow contributors to each segment.
        This will be redundant with to_id, from which it is derived, but
        will also be convenient.
        """
        for segment in self.segment_list:
            from_ids = []
            for _i in range(len(self.to_ids)):
                if (self.to_ids[_i] == segment.id).any():
                    from_ids.append(self.ids[_i])
            segment.set_fromids(from_ids)
            
    def segments_xy_flattened(self):
        """
        Return all points in segments, flattened into a single Numpy array
        """
        allpoints = []
        for segment in self.segment_list:
            allpoints.append(segment.EastingNorthing)
        return np.vstack(allpoints)
        
    def get_bounding_box(self):
        """
        Produces a bounding box in 
        """
        allpoints = segments_xy_flattened(self)
        self.BoundingBox = BoundingBox(allpoints)
        
    def get_segment_ids(self):
        """
        Create a list of all segment IDs
        """
        self.ids = []
        for segment in self.segment_list:
            self.ids.append(segment.id)
        self.ids = np.array(self.ids)

    def get_to_segment_ids(self):
        """
        Create a list of all segment IDs
        """
        self.to_ids = []
        for segment in self.segment_list:
            self.to_ids.append(segment.to_ids)
        self.to_ids_flattened = np.squeeze(np.hstack(self.to_ids))
    
    def compute_offmap_segments(self):
        """
        get list of all segments that flow offmap
        """
        self.get_segment_ids()
        self.offmap_segments = []
        # With the below, we can probbaly just ditch this top part of code
        for segment in self.segment_list:
            if segment.to_ids == 0:
                self.offmap_segments.append(segment.id)
        # More robust -- see which ones do not go to an on-map segment
        for segment in self.segment_list:
            if len([i for i in segment.to_ids if i in self.ids]) == 0:
                self.offmap_segments.append(segment.id)
                
        
    def compute_headwater_segments(self):
        """
        get list of segments that are at the headwaters of a drainage network
        (or the one that is the upstream-most in the case in which an 
        upstream-most segment was defined)
        """
        self.get_to_segment_ids()
        self.headwater_segments = []
        for segment in self.segment_list:
            if segment.id not in self.to_ids_flattened:
                self.headwater_segments.append(segment.id)
    
    def set_hierarchy_order(self):
        """
        Get IDs in order going upstream.
        Helps with calculating distance downstream
        """
        self.get_to_segment_ids()
        self.hierarchy_ids = []
        _ids = [self.downstream_most_id]
        while(len(_ids) > 0):
            self.hierarchy_ids.append(_ids)
            __ids = []
            for _id in _ids:
                __ids.append(self.ids[self.to_ids_flattened == _id])
            __ids = np.hstack(__ids)
            _ids = list(__ids)
    
    def compute_x_in_network(self, overlapping_termini=True):
        """
        Get distance along each stream path
        * Straightforward for downstream at map resolution
        * For upstream, need to march through tributary network to 
          obtain connectivity, and then turn this into a set of distances
        * For reduced resolution (for compute time): need to keep track of 
          endpoints from initial survey
        "Ovelapping termini" means that the nodes between each segment are
        shared among adjacent segments
        
        Can move from bottom up in all cases.
        """
        

        self.overlapping_termini = overlapping_termini
        self.compute_offmap_segments()
        self.set_hierarchy_order()

        # This is ordered from downstream to upstream
        for _id_list in self.hierarchy_ids:
            for _id in _id_list:
                segment = self.segment_list[self.ids == _id][0]
                # Will break in a branching network
                downstream_segment_id = int(segment.to_ids)
                downstream_segment = self.segment_list[self.ids == 
                                                       downstream_segment_id]
                if len(downstream_segment) == 0:
                    x_downstream = self.downstream_x
                elif len(downstream_segment) > 1:
                    sys.exit("Downstream-branching network: not supported!")
                else:
                    downstream_segment = downstream_segment[0]
                    x_downstream = downstream_segment.x[0]
                segment.x = segment.x_local - segment.x_local[-1] \
                                + x_downstream

                
    #def interpolate_xEN_in_network(self):
        
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

    def smooth_window(self):
        """
        Smoothes using a moving window
        """
        for segment in self.segment_list:
            pass
    
    def compute_profile_from_starting_segment(self):
        """
        Compute a single long profile from a starting river segment
        This code currently assumes that you have built the network from 
        just a starting segment, and so doesn't walk downstream, but rather
        takes the segments in order. This could be easily fixed by checking 
        the indices of the next river downstream.
        """
        outlist = []
        for segment in self.segment_list:
            for i in range(len(segment.z)):
                print [segment.x[i], segment.Easting[i], segment.Northing[i], segment.z[i]]
                outlist.append([segment.x[i], segment.Easting[i], segment.Northing[i], segment.z[i]])
        self.long_profile_output = outlist
        
        
