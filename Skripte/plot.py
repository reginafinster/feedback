#!/usr/bin/python3
# coding: utf8

import numpy as np
import math as m
import subprocess
from matplotlib import pyplot as plt


X = np.arange(-2,2,0.1)
Y = np.array(np.sqrt(np.exp(2.0*X)))
YY = np.array(np.sqrt(np.exp(-2.0*X)))


plt.plot(X,Y);
plt.plot(X,YY);
plt.show()




















