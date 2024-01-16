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


directory = os.getcwd()
resultsdir = os.path.abspath("results")
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
   
# sort list of list first by tau and next by phi
sorted_valuelist = sorted(valuelist, key=itemgetter(2, 3))

# get list of relevant params for 2d plot
tau_values = []
tau_values.append(sorted_valuelist[1][2])
detector_list = []
check_list = []
phi_list = []
detector_list_list = []
check_list_list = []
phi_list_list = []
total_excit = sorted_valuelist[0][-3]
for m in range(len(sorted_valuelist)):
    detector_list.append(sorted_valuelist[m][-5])
    check_list.append(sorted_valuelist[m][-2])
    phi_list.append(sorted_valuelist[m][3])
    if m != (len(sorted_valuelist)-1):
        if sorted_valuelist[m][2] != sorted_valuelist[m+1][2]:
            tau_values.append(sorted_valuelist[m+1][2])
            detector_list_list.append(detector_list)
            check_list_list.append(check_list)
            phi_list_list.append(phi_list)
            detector_list = []
            check_list = []
            phi_list = []

# get ranges of phi and tau
max_phi = max([sublist[3] for sublist in valuelist])
min_phi = min([sublist[3] for sublist in valuelist])

max_tau = max([sublist[2] for sublist in valuelist])
min_tau = min([sublist[2] for sublist in valuelist])

# write results to file and check dimensions for plotting
dimY = 0
path_resultsfile = os.path.join(resultsdir,("tau_lines.dat"))
with open(path_resultsfile, 'w+') as f:
    f.write("phi\t\t")
    for value in tau_values:
        f.write("signal at tau=%f\t\t" %value)
    f.write("\n")
    # m ist Anzahl der Spalten
    for m in range(len(detector_list_list[0])):
        f.write("%f\t\t" % phi_list_list[0][m])
        for n in range(len(tau_values)-2):
            f.write("%f\t\t" % detector_list_list[n][m])
        f.write("\n")
f.close()
