#!/usr/bin/python3
# coding: utf8

import sys
import numpy as np
import math as m
import os
from os import listdir
from os.path import isfile, join
from operator import itemgetter

import matplotlib
import pylab as plt

np.set_printoptions(formatter={"float_kind": lambda x: "%g" % x})

directory = os.getcwd()
resultsdir = directory
if not os.path.exists(resultsdir): os.mkdir(resultsdir)

#read in steadystate file in directory
try:
    steadyfile = [f for f in listdir(directory) if (isfile(join(directory, f)) and f.endswith("steadylist.dat"))]
except:
    print(("check input dir {}. The following error occured").format(directory))
    raise

# read in file 
valuelist = []
datafile = open(os.path.join(directory, steadyfile[0]), "r")
for i, line in enumerate(datafile):
    if line != '\n':
        if (i == 0):
            header = line
        else:
            readdata = line.split('\t\t')
            readdata = readdata[:-1]
            valuelist.append(list(map(lambda x: float(x), readdata)))
datafile.close()
   
# sort list of list first by phi and next by tau
sorted_valuelist = sorted(valuelist, key=itemgetter(3, 2))

# get list of relevant params for 2d plot
phi_values = []
phi_values.append(sorted_valuelist[1][3])
detector_list = []
check_list = []
tau_list = []
detector_list_list = []
check_list_list = []
tau_list_list = []
total_excit = sorted_valuelist[0][-3]
for m in range(len(sorted_valuelist)):
    detector_list.append(sorted_valuelist[m][-5])
    check_list.append(sorted_valuelist[m][-2])
    tau_list.append(sorted_valuelist[m][2])
    if m != (len(sorted_valuelist)-1):
        if sorted_valuelist[m][3] != sorted_valuelist[m+1][3]:
            phi_values.append(sorted_valuelist[m+1][3])
            detector_list_list.append(detector_list)
            check_list_list.append(check_list)
            tau_list_list.append(tau_list)
            detector_list = []
            check_list = []
            tau_list = []

# get ranges of phi and tau
max_phi = max([sublist[3] for sublist in valuelist])
min_phi = min([sublist[3] for sublist in valuelist])

max_tau = max([sublist[2] for sublist in valuelist])
min_tau = min([sublist[2] for sublist in valuelist])

# write results to file 
path_resultsfile = os.path.join(resultsdir,("tau_lines.dat"))
with open(path_resultsfile, 'w+') as f:
    f.write("# tau\t\t")
    for value in phi_values:
        f.write("phi=%.2f\t\t" %value)
    f.write("\n")
    # m ist Anzahl der Spalten
    for m in range(len(detector_list_list[0])):
        f.write("%f\t\t" % tau_list_list[0][m])
        for n in range(len(phi_values)-1):
            f.write("%f\t\t" % detector_list_list[n][m])
        f.write("\n")
f.close()
