#!/usr/bin/python3
# coding: utf8

import sys
import numpy as np
import math as m
import os

dataname1 = sys.argv[1]
dataname2 = sys.argv[2]
datafile1 = open(dataname1, "r")
datafile2 = open(dataname2, "r")

valuelist1 = []
valuelist2 = []

for i, line in enumerate(datafile1):
    if line != '\n':
        readdata       = line.split('\t\t')
        time1          = float(readdata[0])
        value1         = float(readdata[1])
        valuelist1.append((time1,value1))
datafile1.close()

for i, line in enumerate(datafile2):
    if line != '\n':
        readdata       = line.split('\t\t')
        time2          = float(readdata[0])
        value2         = float(readdata[1])
        valuelist2.append((time2,value2))
datafile2.close()

deviation_list = []
i=0;
while (i < 40000):
    # make sure we deal with the values for the same timestep:
    if (valuelist1[i][0] == valuelist2[i][0]):
        # calculate difference from reference file (file1) and append it to ouput list with the respective timestep:
        deviation_list.append([valuelist1[i][0],(valuelist1[i][1] - valuelist2[i][1])])
        i+=1
    else:
        printf("Timesteps are not chosen equally in both files, cannot make calculation")
        break

currentfilename = os.path.abspath(dataname2[:-4])

path_resultsfile = os.path.join(str(currentfilename + "_deviation.dat"))
with open(path_resultsfile, 'w+') as f:
    for m in range(len(deviation_list)):
        for n in range(len(deviation_list[m])):
            f.write("%f\t\t" % deviation_list[m][n])
        f.write("\n")
f.close()





