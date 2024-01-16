#!/usr/bin/python3
# coding: utf8

import sys
import numpy as np
import math as m
import os
from os import listdir
from os.path import isfile, join


directory = os.path.abspath(sys.argv[1])
resultsdir = os.path.abspath("averaged_results")
if not os.path.exists(resultsdir): os.mkdir(resultsdir)

try:
    allfiles = [f for f in listdir(directory) if (isfile(join(directory, f)) and f.endswith(".dat"))]
    #allfiles = [f for f in allfiles if f.endswith(".dat")]
except:
    print(("check input dir {}. The following error occured").format(directory))
    raise

print(len(allfiles))
#valuematrix = [0 for x in range(len(allfiles))]
valuelist = []
firstfile = open(os.path.join(directory,allfiles[0]), "r")
for i, line in enumerate(firstfile):
    if line != '\n':
            if (i in range(0,17)):
                print(i)
                continue
            else:
                #values_temp_list = list(map(lambda x: float(x), line.split('\t\t')))
                readdata = line.split("\t")[:-1]
                valuelist.append(list(map(lambda x: float(x), readdata)))
firstfile.close()
for k in range(1,len(allfiles)):
    datafile = open(os.path.join(directory,allfiles[k]), "r")
    for i, line in enumerate(datafile):
        if line != '\n':
            if (i in range(0,17)):
                continue
            else:
                readdata = line.split('\t')[:-1]
                for j in range(len(readdata)):
                    value = float(readdata[j])
                    valuelist[i-17][j] += value
    datafile.close()
valuearray = np.array(valuelist)
valuearray /= len(allfiles)

average_valuelist = [0.0] *len(valuearray)
for i in range(len(valuearray)):
    average_valuelist[i] += sum(valuearray[i][1:])
totalvaluearray = np.array(average_valuelist)
print(totalvaluearray)
print(valuearray[1][1:])
totalvaluearray /= len(valuelist[1][1:])

path_resultsfile = os.path.join(resultsdir, ("av_res_" + sys.argv[1] + ".dat"))
path_averagefile = os.path.join(resultsdir, ("av_res_allsites_" + sys.argv[1] + ".dat"))

with open(path_resultsfile, 'w+') as f:
    for m in range(len(valuearray)):
        for n in range(len(valuearray[m])):
            f.write("%f\t\t" % valuearray[m][n])
        f.write("\n")
f.close()

#with open(path_averagefile, 'w+') as f:
    #for m in range(len(valuearray)):
        #f.write("%f\t\t" % valuearray[m][0])
        #f.write("%f\t\t" % totalvaluearray[m])
        #f.write("\n")
#f.close()





