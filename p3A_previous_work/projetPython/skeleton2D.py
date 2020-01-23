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

class Skeleton2D:

    def __init__(self, p, children=[]):
        self.pos = p
        self.children = children

    def add_child(self, t):
        self.children.append(t)

    def transformation(self, m):
        """Transformation in homogenous coordinates"""
        self.pos = m.dot(self.pos)
        for c in self.children:
            c.transformation(m)
        return

    def list_abs_position(self):
        L = [self.pos]
        for c in self.children:
            L.extend(c.list_abs_position())
            return L

    def plot_aux(self):
        root_abs = self.pos
        t = False
        for c in self.children:
            c_abs = c.pos
            if (root_abs[2] != 0 and c_abs[2]):
                print(root_abs[2],c_abs[2])
                plt.plot([root_abs[0]/root_abs[2], c_abs[0]/c_abs[2]],  [root_abs[1]/root_abs[2], c_abs[1]/c_abs[2]])
                c.plot_aux()
                t = True
        return t

    def plot(self):
        fig = plt.figure()
        if self.plot_aux():
            plt.show()
