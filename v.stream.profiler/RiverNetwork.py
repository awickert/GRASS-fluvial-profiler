#! /usr/local/env python

# 06 May 2018
# ADW

import numpy as np
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
        self.E = None
        self.N = None
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
            self.Easting = self.EastingNorthing[:,1]
            self.Northing = self.EastingNorthing[:,1]
        elif (Easting is not None) and (Northing is not None):
            self.Easting = self.Easting
            self.Northing = self.Northing
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
        else:
            warnings.warn("Both Easting and Northing must be defined")
        
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
        segment flows offmap). Thid defaults to 0 following GRASS GIS 
        convention.
        """
        self.segment_list = segment_list
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
        self.to_ids_flattened = np.hstack(self.to_ids)
    
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
        headwater_segments = []
        for segment in self.segment_list:
            if segment.id not in self.to_ids_flattened:
                headwater_segments.append(segment.id)
    
    def compute_x_in_network(self, overlapping_termini=True):
        """
        Get distance along each stream path
        "Ovelapping termini" means that the nodes between each segment are
        shared by 
        """
        # Flag so we know
        # e.g., for numerical modeling on the network
        self.overlapping_termini = overlapping_termini
        # Find all offmap segments
        self.compute_offmap_segments()
        # Find segments that flow to these, etc. (see other code)
        for terminus in self.offmap_segments:
            x_downstream = [0.] # start
            terminal_ids = [terminus]
            while len(terminal_ids) > 0:
                x_downstream_next = []
                terminal_ids_next = []
                _j = 0
                # March upstream, replacing lists of ids and associated
                # downstream x-values as you go
                for _id in terminal_ids:
                    # Update x
                    segment = np.array(self.segment_list)[self.ids == _id][0]
                    segment.x = segment.x - segment.x[-1] + x_downstream[_j]
                    # Updates for upcoming round
                    for _i in range(len(self.to_ids)):
                        # Next upstream segment IDs
                        if (self.to_ids[_i] == _id).any():
                            terminal_ids_next.append(self.ids[_i])
                        # Update the list of x values at the end of the next set
                        if self.overlapping_termini:
                            # To produce the same count
                            if (self.to_ids[_i] == _id).any():
                                x_downstream_next.append(segment.x[0])
                        else:
                            sys.exit("I don't know how to deal with this case")
                    # segment.set_fromids(from_ids) # or append? unnecessary...
                    _j += 1
                # Replace main arrays with temp arrays
                #print terminal_ids
                #print x_downstream
                terminal_ids = terminal_ids_next
                x_downstream = x_downstream_next
        
    def compute_profile_from_starting_segment(self):
        """
        Compute a single long profile
        Note that this may not be easy if there is a fork in a river!
        Note to self: consider simplifying network geometry to just
        what I will need (tributaries joining)
        Although I would prefer to keep it general, philosophically...
        """
        pass
