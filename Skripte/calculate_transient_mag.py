#!/usr/bin/python3
# coding: utf8

import sys
import numpy as np
import math as m

dataname1 = sys.argv[1]
line1 = int(sys.argv[2])
dataname2 = sys.argv[3]
line2 = int(sys.argv[4])
datafile1 = open(dataname1, "r")
datafile2 = open(dataname2, "r")

valuelist1 = []
valuelist2 = []

for i, line in enumerate(datafile1):
    if line != '\n':
        if i == 0:
            continue
        else:
        #read data
            readdata       = line.split('\t\t')
            time1          = float(readdata[0])
            value1         = float(readdata[line1])
            
            #write data in lists
            valuelist1.append((time1,value1))
datafile1.close()

for i, line in enumerate(datafile2):
    if line != '\n':
        if i == 0:
            continue
        else:
        #read data
            readdata       = line.split('\t\t')
            time2          = float(readdata[0])
            value2         = float(readdata[line2])
            
            #write data in lists
            valuelist2.append((time2,value2))
datafile2.close()

# calculate square errors
percentage_deviation_list = []
valuelist = []
i=0;
j=0;
while (j < len(valuelist2)):
    if (valuelist1[i][0] == valuelist2[j][0]):
        valuelist.append((valuelist1[i][1], valuelist2[j][1]))
        j+=1
    i += 1
for val in valuelist:
    error = m.sqrt((val[0] - val[1])**2)
    if (min(m.fabs(val[0]),m.fabs(val[1])) == 0.0):
        percentage_deviation = 0.0
    else:
        percentage_deviation = error/(min(m.fabs(val[0]),m.fabs(val[1])))
    percentage_deviation_list.append(percentage_deviation)

# calculate mean and maximal square error
percentage_deviation_list.remove(max(percentage_deviation_list))
percentage_deviation_array = np.array(percentage_deviation_list)



mean_percentage_deviation = sum(percentage_deviation_array)/len(percentage_deviation_array)
max_percentage_deviation = max(percentage_deviation_array)

print("mean deviation: %.1f percent, max deviation: %.1f percent" % (mean_percentage_deviation*100.0, max_percentage_deviation*100.0))





