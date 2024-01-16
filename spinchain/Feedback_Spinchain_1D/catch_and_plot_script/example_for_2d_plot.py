#!/usr/bin/python3
# coding: utf8

import numpy as np
import math as m
import subprocess
import matplotlib
from matplotlib import cm
matplotlib.use('Qt5Agg')
import pylab as plt

#this file is for plotting a bifurcation diagram for the phase. bifurcation parameter is either coupling strength or pump current

tauList           = []
philist           = []
detectorlist      = []
    
datafile = open("sortiert.dat", "r")   
for i, line in enumerate(datafile):
    if line != '\n':
        if i == 0:
            continue
        #read number of steps (needed for reshaping matrix)
        elif i == 1:
            readdata       = line.split('    ')
            kappa_number   = int(readdata[0])
            C_number       = int(readdata[1])
            pump_number    = int(readdata[2])
        
        else:
        #read data
            readdata       = line.split('    ')
            C              = float(readdata[0])
            kappa          = float(readdata[1])
            pump           = float(readdata[2])
            cosinus        = float(readdata[7])
            
            #write data in lists
            Clist.append(C)
            kappalist.append(kappa)
            pumplist.append(pump)
            if cosinus == 2.0:
                cosinus = np.nan
            cosinuslist.append(cosinus)
datafile.close()

#get minimal and maximal value of pumpcurrent

Cmin = min(Clist)
Cmax = max(Clist)

pumpmin = pumplist[0]
pumpmax = pumplist[-1]

kappa = kappalist[0]

#convert list to array for plotting and reshape:
cosinusarray = np.array(cosinuslist)
cosinusarray = np.reshape(cosinusarray, (pump_number, (C_number-1)))
cosinusarray = np.flipud(cosinusarray)

#plot with pumpcurrent as bif par, plot only part aroung Hopf-bifurcation:
font = {'size' : 14}
plt.rc('font', **font)
plt.figure(figsize = (8,4))
plt.imshow(cosinusarray, interpolation='none', aspect='auto', extent=[Cmin,Cmax,pumpmin,pumpmax], cmap = cm.Set1)
plt.legend()
plt.xlabel('coupling phase $C$')
plt.ylabel('pump current')
plt.title('Cosine of $\Delta \phi$\n feedback strength = %.2f'%kappa)
plt.xticks([Cmin,(Cmax-Cmin)/4,(Cmax-Cmin)/2,(3*(Cmax-Cmin))/4,Cmax],[r'$0$',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2\pi$'])
cb = plt.colorbar()
cb.set_label('Cos($\Delta \phi$)')
plt.savefig("pump_cosine.pdf", dpi = 300,  bbox_inches='tight', pad_inches=0)

#optional:plot complete range
#plt.imshow(cosinusarray, interpolation='none', aspect='auto', extent=[0,m.pi*2,pumpmin,pumpmax])
#plt.legend()
#plt.xlabel('coupling phase')
#plt.ylabel('pumpcurrent')
#plt.title('2-laser-network, value of cosinus of phase')
#plt.xticks([0,m.pi/2,m.pi,3*m.pi/2,2*m.pi],[r'$0$',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2\pi$'])
#plt.colorbar()
#plt.show()




















