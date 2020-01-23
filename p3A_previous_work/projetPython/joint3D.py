# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:24:28 2017

@author: Manon
"""
import numpy as np

class Joint3D:
    
    def __init__(self, position, quaternion, name=None):
        assert (position.shape == (3,))
        
        self.position = position
        self.abs_pos = None
        self.quaternion = quaternion
        self.name = name

    def transformation(self, m):
        p = self.position
        assert(m.shape==(p.shape,p.shape))
        self.position = m.dot(p)
        return 
    
    
