#!/usr/bin/env python3
# coding:utf-8

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc("font",**{"family":"serif","serif":["Palatino"]})
rc('text', usetex=True)

PUx = './Ux/'
PUy = './Uy/'
PP  = './P/'

files = sorted(os.listdir(PP))
imgs = sorted(os.listdir("./anim/"))
for f in files:
    if f.split(".")[0]+".png" not in imgs:
        Ux = np.gradient(np.loadtxt(PUx+f), axis=0)
        Uy = np.gradient(np.loadtxt(PUy+f), axis=1)
        plt.imsave('./anim/'+f.split(".")[0]+".png",Uy-Ux, vmin=-2,vmax=2, cmap='jet')
#    a = np.loadtxt(PP+f)
#    if np.isnan(a).any():
#        break
#    plt.imshow(a, cmap='jet')
#    plt.pause(0.1)
#    plt.cla()
#    plt.clf()
    #print(np.loadtxt(PUx+f).mean())
