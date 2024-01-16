#!/usr/bin/python3
# coding: utf8

import numpy as np
import os
import math as m
import subprocess
#from matplotlib import pyplot as plt


#X = np.arange(0.0,2.1,0.1)
#Y = np.array(1 - (X/((X+1)*(X+1))));
Gamma1 = 1.0
Gamma2Min = 0.5
Gamma2Max = 2.0
X = np.arange(Gamma2Min,Gamma2Max,0.1)
Y = np.array(0.5*6.28*(1 - ((4.0*Gamma1*X)/((Gamma1+X)*(Gamma1+X)))))

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




















