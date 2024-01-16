#!/usr/bin/python3
# coding: utf8

import sys
import numpy as np
import math as m
import os
from os import listdir
from os.path import isfile, join

dataname1 = sys.argv[1]
directory = sys.argv[2]

try:
    allfiles = [f for f in listdir(directory) if (isfile(join(directory, f)) and f.endswith(".dat"))]
except:
    print(("check input dir {}. The following error occured").format(directory))
    raise

datafile1 = open(dataname1, "r")
valuelist1 = []
for i, line in enumerate(datafile1):
    if line != '\n':
        readdata       = line.split('\t\t')
        time1          = float(readdata[0])
        value1         = float(readdata[1])
        valuelist1.append((time1,value1))
datafile1.close()

phi_sat_list = []

for onefile in allfiles:
    datafile2 = open(os.path.join(directory,onefile), "r")
    valuelist2 = []
    for i, line in enumerate(datafile2):
        if line != '\n':
            readdata       = line.split('\t\t')
            time2          = float(readdata[0])
            value2         = float(readdata[1])
            valuelist2.append((time2,value2))
    datafile2.close()

    deviation_list = []
    i=0;
    while (i < 40000): # WARNING hier habe ich auf T=100 abgeschnitten
        if (valuelist1[i][0] == valuelist2[i][0]):
            deviation_list.append([valuelist1[i][0],(valuelist1[i][1] - valuelist2[i][1])])
            i+=1
        else:
            printf("Timesteps are not chosen equally in both files, cannot make calculation")
            break
    print([onefile[-6:-4],deviation_list[-1][1]])
    phi_sat_list.append([float(onefile[-6:-4]),deviation_list[-1][1]])
    print(phi_sat_list)

    currentfilename = os.path.abspath(onefile[:-4]) # WARNING hier schneide ich die Endung .dat ab
    path_resultsfile = os.path.join(str(currentfilename + "_deviation.dat"))
    with open(path_resultsfile, 'w+') as f:
        for m in range(len(deviation_list)):
            for n in range(len(deviation_list[m])):
                f.write("%f\t\t" % deviation_list[m][n])
            f.write("\n")
    f.close()

phi_sat_filename = os.path.join(os.getcwd(),"phi_sat.dat")
with open(phi_sat_filename, 'w+') as f:
    for m in range(len(phi_sat_list)):
        for n in range(len(phi_sat_list[m])):
            f.write("%f\t\t" % phi_sat_list[m][n])
        f.write("\n")
f.close()



