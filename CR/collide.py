#!/usr/bin/env python3
# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#plt.style.use('dark_background')

rc("font",**{"family":"serif","serif":["Palatino"]})
rc('text', usetex=True)

L = 4
H = 3

b = [0]*L
m = [1]*L
t = [2]*L
