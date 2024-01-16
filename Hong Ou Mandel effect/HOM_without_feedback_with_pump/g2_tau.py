#!/usr/bin/python3
# coding: utf8

import numpy as np
import os
import math as m
import subprocess
#from matplotlib import pyplot as plt


Gamma1 = 1.0
Gamma2 = 0.5
td = 0.0
prefactor = 0.25*Gamma1*Gamma1*Gamma2*Gamma2*np.exp(-2.0*Gamma1*Gamma2*td)
X = np.arange(0.0,100.0,0.1)
Y = np.array( prefactor*(np.exp(-2.0*Gamma1*X) + np.exp(-2.0*Gamma2*X) - 2.0*np.exp(-(Gamma1+Gamma2)*X)) )
t_part = np.array(prefactor*(np.exp(-2.0*Gamma1*X)))
r_part = np.array(prefactor*(np.exp(-2.0*Gamma2*X)))
rt_part = np.array(prefactor*(1.0*np.exp(-(Gamma1+Gamma2)*X)))

resultsdir = os.path.abspath("results")
if not os.path.exists(resultsdir): os.mkdir(resultsdir)
path_results = os.path.join(resultsdir, "gamma_tau.dat")


with open(path_results, 'w+') as f:
    for m in range(len(X)):
        f.write("%f\t" % X[m])
        f.write("%f\t" % Y[m])
        f.write("%f\t" % t_part[m])
        f.write("%f\t" % r_part[m])
        f.write("%f" % rt_part[m])
        f.write("\n")
f.close()


#plt.plot(X,Y);
#plt.show()




















