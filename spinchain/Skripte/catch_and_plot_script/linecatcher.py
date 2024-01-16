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
from matplotlib import cm
#matplotlib.use('Qt5Agg')
import pylab as plt
from collections import OrderedDict

cmaps = OrderedDict()


directory = os.getcwd()
currdirectory = os.path.join(directory,"currcontainer")
resultsdir = os.path.abspath("summary_stst")
if not os.path.exists(resultsdir): os.mkdir(resultsdir)

#read all steadystate files in directory
try:
    allfiles = [f for f in listdir(directory) if (isfile(join(directory, f)) and f.endswith("steady.dat"))]
except:
    print(("check input dir {}. The following error occured").format(directory))
    raise


valuelist = []
for file in allfiles:
    datafile = open(os.path.join(directory, file), "r")
    for i, line in enumerate(datafile):
        if line != '\n':
            if (i == 0):
                header = line
            else:
                readdata = line.split('\t\t')
                numberofcols = len(readdata)
                valuelist.append(list(map(lambda x: float(x), readdata)))
    datafile.close()
   
# sort list of list first by tau and next by phi
sorted_valuelist = sorted(valuelist, key=itemgetter(2, 3))

# get list of relevant params for 2d plot
detector_list = []
steady_list_list = []
check_list = []
steady_list = []
total_excit = sorted_valuelist[0][-3]
for m in range(len(sorted_valuelist)):
    detector_list.append(sorted_valuelist[m][-5])
    check_list.append(sorted_valuelist[m][-2])

# get ranges of phi and tau
max_phi = max([sublist[3] for sublist in valuelist])
min_phi = min([sublist[3] for sublist in valuelist])

max_tau = max([sublist[2] for sublist in valuelist])
min_tau = min([sublist[2] for sublist in valuelist])

# write results to file and check dimensions for plotting
dimY = 0
path_resultsfile = os.path.join(resultsdir,("steadylist.dat"))
with open(path_resultsfile, 'w+') as f:
    f.write(header)
    for m in range(len(sorted_valuelist)):
        steady_list.append(sorted_valuelist[m][7:-5])
        for n in range(len(sorted_valuelist[m])):
            f.write("%f\t\t" % sorted_valuelist[m][n])
        if m != (len(sorted_valuelist)-1):
            if sorted_valuelist[m][2] != sorted_valuelist[m+1][2]:
                steady_list_list.append([steady_list,sorted_valuelist[m][2]])
                f.write("\n")
                dimY += 1
                if dimY==1:
                    dimX = m+1
                steady_list = []
        else:
            steady_list_list.append([steady_list,sorted_valuelist[m][2]])
        f.write("\n")
            
dimY += 1
f.close()

# make 2d plots
detector_array = np.array(detector_list)
detector_array = np.reshape(detector_array, (dimY, dimX))
detector_array = np.rot90(detector_array)

normalized_detector_array = detector_array/total_excit

check_array = np.array(check_list)
check_array = np.reshape(check_array, (dimY, dimX))
check_array = np.rot90(check_array)

#plt.imshow(detector_array, interpolation='none', aspect='auto', extent=[min_tau,max_tau,min_phi,max_phi], vmin=0, vmax=total_excit)
#plt.legend()
#plt.xlabel('delay time J$\\tau$')
#plt.ylabel('feedback phase $\phi$')
#plt.title('integrated detector signal')
#cb = plt.colorbar()
#cb.set_label('value of integrated signal')
#plt.savefig("detector.pdf", dpi = 300,  bbox_inches='tight', pad_inches=0)

plt.imshow(normalized_detector_array, interpolation='none', aspect='auto', extent=[min_tau,max_tau,min_phi,max_phi], vmin=0, vmax=1, cmap = cm.gist_heat)
plt.legend()
plt.xlabel('delay time J$\\tau$')
plt.ylabel('feedback phase $\phi$')
plt.title('integrated detector signal')
#cb.remove()
plt.xticks([5,10,15,20,25,30,35,40,45,50],[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
cb = plt.colorbar()
cb.set_label('value of integrated signal')
plt.savefig("normalized_detector.pdf", dpi = 300,  bbox_inches='tight', pad_inches=0)

plt.imshow(check_array, interpolation='none', aspect='auto', extent=[min_tau,max_tau,min_phi,max_phi])
plt.legend()
plt.xlabel('delay time J$\\tau$')
plt.ylabel('feedback phase $\phi$')
plt.title('steady state signal at detector')
cb.remove()
cb = plt.colorbar()
cb.set_label('value of integrated signal')
plt.savefig("check.pdf", dpi = 300,  bbox_inches='tight', pad_inches=0)

#for steadylist in steady_list_list:
    #steady_array = np.array(steadylist[0])
    #plt.imshow(steady_array, interpolation='none', aspect='auto', extent=[1,len(steadylist[0][0]),min_phi,max_phi])
    #plt.legend()
    #plt.xlabel('site position')
    #plt.ylabel('feedback phase $\phi$')
    #plt.title('steady state magnetization')
    #cb.remove()
    #cb = plt.colorbar()
    #cb.set_label('magnetization')
    #filename = "magprofile_tau_" + str(steadylist[1]) + ".pdf"
    #plt.savefig(filename, dpi = 300,  bbox_inches='tight', pad_inches=0)




