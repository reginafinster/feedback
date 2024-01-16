#!/usr/bin/python3
# coding: utf8

import sys
import numpy as np
import math as m
import os
from os import listdir
from os.path import isfile, join


directory = os.getcwd()
currdirectory = os.path.join(directory,"currcontainer")
resultsdir = os.path.abspath("summary_stst")
if not os.path.exists(resultsdir): os.mkdir(resultsdir)

#first we read in the steadystate values and parameters and sum them on the valuelist
try:
    allfiles = [f for f in listdir(directory) if (isfile(join(directory, f)) and f.endswith("steady.dat"))]
    #allfiles = [f for f in allfiles if f.endswith(".dat")]
except:
    print(("check input dir {}. The following error occured").format(directory))
    raise

print(len(allfiles))
#start read
valuelist = []
numberofcols = 0
firstfile = open(os.path.join(directory,allfiles[0]), "r")
for i, line in enumerate(firstfile):
    if line != '\n':
        if (i == 0):
            continue
        else:
            readdata = line.split("\t\t")[:-1]
            numberofcols = len(readdata)
            valuelist.append(list(map(lambda x: float(x), readdata)))
firstfile.close()
print(numberofcols)
numberofrows = len(allfiles)
valuearray = [[0 for x in range(numberofcols)] for y in range(numberofrows)]
print(valuelist)
print(valuearray[0])
for i in range(numberofcols):
    valuearray[0][i] = valuelist[0][i]
for k in range(1,len(allfiles)):
    datafile = open(os.path.join(directory,allfiles[k]), "r")
    for i, line in enumerate(datafile):
        if line != '\n':
            if (i == 0):
                continue
            else:
                readdata = line.split('\t\t')[:-1]
                for j in range(len(readdata)):
                    value = float(readdata[j])
                    print(value)
                    print("k%i" % k)
                    valuearray[k][j] += value
                    print(valuearray)
    datafile.close()
#start write    
path_resultsfile = os.path.join(resultsdir,("steadylist.dat"))
with open(path_resultsfile, 'w+') as f:
    f.write("# N \t\t Gamma \t\t tau \t\t phi \t\t init \t\t init_pos \t\t shifted")
    N = len(valuearray[0])-11
    for i in range(N):
        f.write(" \t\t occupation%i" % (N-i))
    f.write(" \t\t detector \t\t trapped-bins \t\t det+sumocc \t\t steady_check")
    f.write("\n")
    for m in range(len(valuearray)):
        for n in range(len(valuearray[m])):
            f.write("%f\t\t" % valuearray[m][n])
        f.write("\n")
f.close()
