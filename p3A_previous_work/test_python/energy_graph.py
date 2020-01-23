#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 10:43:42 2018

@author: coraliedavid
"""

import numpy as np
import os
import matplotlib.pyplot as plt

colours = ['r', 'b', 'g', 'gray', 'pink', 'orange']


fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Energy graph")
ax1.set_xlabel('Iterations')
ax1.set_ylabel('Energy')
#with open("energy.txt") as f:
#    lines = f.readlines()
#    for line in lines:
#        try:
#            l = line.split(";")
#            plt.plot(l[:-1])
#        except:
#           print(line)

energy = []
T = []
with open("energy.txt") as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        try:
            l = line.split(";")
            T.append(float(l[2]))
            energy.append(l[1])
            if l[0]=="0":
               plt.plot(i, 0, 'r+')
            if l[0]=="1":
               plt.plot(i, 0.1, 'b+')
        except:
           print(line)
           
plt.plot(energy)
plt.xlim([00,3000])
#leg = ax1.legend()
plt.figure()
plt.show()