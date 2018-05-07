#! /usr/local/env python

# 06 May 2018
# ADW

import numpy as np
import warnings

class Segment(object):

    def __init__(self, _id, to_id=None, segment_type='stream', 
                 direction='to_downstream'):
        """
        A segment of a stream, between tributary junctions
        """
        self.id = _id
        self.to_id = to_id
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
        
    ###########
    # SETTERS #
    ###########

    def set_x(self):
        """
        Set distances downstream
        """
        self.x = x
        
    def set_z(self, z):
        """
        Set elevation
        """
        self.z = z
        
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
        
    def set_fromid(self, from_id):
        """
        Set the ID of the stream from which 
        """
        self.from_ids = from_id

    def append_fromid(self, from_id):
        """
        Set the ID of the stream from which 
        """
        self.from_ids.append(from_id)
        
    def set_E(self, E)
        """
        Set easting of points
        """
        self.E = E
    
    def set_N(self, N)
        """
        Set northing of points
        """
        self.N = N
        
    def calc_x_from_E_N(self, _x):
        """
        Calculate downstream distances from Easting and Northing, relative
        only to the individual segment
        """
        if self.E and self.N:
            dx = ( np.diff(self.E)**2 + np.diff(self.N)**2 )**.5
            self.x = np.hstack((0., np.cumsum(dx)))
        else:
            warnings.warn("Both Easting and Notrthing must be defined")
        
class Network(object):

    def __init__(self, segment_list):
        """
        Tools to build a river network
        """
        self.segment_list = segment_list
        # Convenient info to have in lists
        # identifiers and connections
        self.ids = []
        self.to_ids = []
        for segment in self.segment_list:
            self.ids.append(segment.id)
            self.to_ids.append(segment.to_id)
        self.ids = np.array(self.ids)
        self.to_ids = np.array(self.to_ids)
        
    def compute_fromstreams(self):
        """
        Compute IDs of flow contributors to each segment.
        This will be redundant with to_id, from which it is derived, but
        will also be convenient.
        """
        for segment in self.segment_list:
            self.from_ids = list(self.ids[self.to_ids == segment.id])
            segment.set_fromid(self.from_ids)
            
        
