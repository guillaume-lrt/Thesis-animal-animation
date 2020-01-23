# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:57:20 2017

@author: Manon
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import quaternion
from plot_utils import plot_cylinder
from skeleton2D import Skeleton2D

class SkeletonTree: 
    
    def __init__(self, j, children=[]):
        self.root = j
        self.children = children
    
    def add_child(self, t):
        self.children.append(t)
        
    def update_absolute_pos(self, t=[], q=[]): 
        assert(len(t)==len(q))
        absolute_position = self.root.position
        for i in range(len(t)-1, -1, -1):
            absolute_position = t[i] + quaternion.as_rotation_matrix(q[i]).dot(absolute_position)
        self.root.abs_pos = absolute_position
        #print(self.root.name,": ", absolute_position)
        new_t = t[:]
        new_t.append(self.root.position)
        new_q = q[:]
        new_q.append(self.root.quaternion)
        for c in self.children:
            c.update_absolute_pos(new_t, new_q)
        return
    
    def list_abs_position(self):
        L = []
        L.append(self.root.abs_pos)
        for c in self.children:
            L.extend(c.list_abs_position())
        return L
    
    def plot_absolute(self, ax): 
        root_abs = self.root.abs_pos
        #print(self.root.name,": ",root_abs)
        for c in self.children:
            c_abs = c.root.abs_pos
            ax.plot([root_abs[0], c_abs[0]], [root_abs[1], c_abs[1]], [root_abs[2], c_abs[2]])
            plot_cylinder(ax, root_abs, c_abs)
            c.plot_absolute(ax)


    def plot(self):
        self.update_absolute_pos()
        fig = plt.figure()
        #fig.axis('equal')
        ax = fig.add_subplot(111, projection='3d')
        self.plot_absolute(ax)
        plt.show()
        
    def project(self, p):
        L = []
        for c in self.children:
            c_2D = c.project(p)
            L.append(c_2D)
        return Skeleton2D(p(self.root.abs_pos), L)

        
