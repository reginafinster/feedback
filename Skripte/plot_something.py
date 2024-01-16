#!/usr/bin/python3
# coding: utf8

import numpy as np
import os
import math as m
import subprocess
#from matplotlib import pyplot as plt


X = np.arange(1.0,2.1,0.1)
Y = np.array(X/((X+1)*(X+1)));


resultsdir = os.path.abspath("results")
if not os.path.exists(resultsdir): os.mkdir(resultsdir)
path_results = os.path.join(resultsdir, "gamma.dat")


with open(path_results, 'w+') as f:
    for m in range(len(X)):
        f.write("%f\t\t" % X[m])
        f.write("%f" % Y[m])
        f.write("\n")
f.close()


#plt.plot(X,Y);
#plt.show()




















