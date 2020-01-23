#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 21:31:41 2018

@author: coraliedavid
"""

import matplotlib.pyplot as plt
import numpy as np

X = list(range(2, 25))
#for n in [1,10, 100, 1000, 1e4, 1e6, 1e10, 1e30]:
#    Y = [1-(1-1./np.math.factorial(x-2))**n for x in X]
#    Y2 = [Y[i-1]-Y[i] for i in range(len(Y))]
#    plt.plot(X, Y, label = str(n))
   
#plt.xlim(0, 25)
#plt.ylim(-.01, 1.01)
#plt.legend(title="Number of iterations")
#plt.title("Convergence")
#plt.xlabel("Number of nodes")
#plt.ylabel("Probability of reaching the optimal path")
#plt.axhline(0.99, ls='--')
#plt.axhline(0.95, ls='--')

for a in [0.95, 0.99]:
    Y = [int(np.log(1-a)/np.log(1-1./np.math.factorial(x-2))) for x in X]
    #Y2 = [Y[i-1]-Y[i] for i in range(len(Y))]
    plt.plot(X, Y, label = str(a))

plt.show()
