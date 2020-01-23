# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 18:08:21 2017

@author: Manon
"""
import cv2
import numpy as np
from joint3D import Joint3D
import quaternion
from skeletonTree import SkeletonTree
from quat_utils import *
from projection_utils import *
from fit import fit

def main():
    un = np.quaternion(1, 0, 0, 0)
    r = quat_from_alpha_vec(np.pi/4, [1, 0, 0])
    head = Joint3D(np.array([0, -0.5, 1]), r, "head")
    nose = Joint3D(np.array([0, -1, 0]), un, "nose")
    neck = Joint3D(np.array([0, -5, 0]), un, "neck")
    leftfoot = Joint3D(np.array([-0.5, -1, -1]), un, "leftfoot")
    rightfoot = Joint3D(np.array([0.5, -1, -1]), un, "rightfoot")
    hip = Joint3D(np.array([0, 0, 0]), un, "hip")
    leftbfoot = Joint3D(np.array([-0.5, 1, -1]), un, "leftbfoot")
    rightbfoot = Joint3D(np.array([0.5, 1, -1]), un, "rightbfoot")
    tail = Joint3D(np.array([0,1,2]), un, "tail")


    noset = SkeletonTree(nose)
    headt = SkeletonTree(head, [noset])

    lfn = SkeletonTree(leftfoot)
    rfn = SkeletonTree(rightfoot)
    neckt = SkeletonTree(neck, [headt, lfn, rfn])

    tailn = SkeletonTree(tail)
    lbfn = SkeletonTree(leftbfoot)
    rbfn = SkeletonTree(rightbfoot)
    hipt = SkeletonTree(hip, [neckt, tailn, lbfn, rbfn])

    hipt.plot()
    d = np.array([0, 0, 1]) #np.random.randint(1,10,3)
    d = d.astype(float)/np.linalg.norm(d)
    p = projection_on_plane(d)
    s2D = hipt.project(p)

    s2D.plot()

    # fit(p,d)

if __name__=="__main__":
    print("Hello")
main()
