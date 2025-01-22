#!/usr/bin/env python3
# coding:utf-8

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc("font",**{"family":"serif","serif":["Palatino"]})
rc('text', usetex=True)

PU = './U/'
PP  = './P/'
files = sorted(os.listdir(PU))
U = np.zeros((len(files),392,490,2))
for i, f in enumerate(files):
    U[i] = np.fromfile(f"{PU}{f}", dtype=np.float64, offset=4).reshape((392, 490, 2))
U = np.sqrt(np.sum(U*U,-1))
vmax = U.max()

for i, u in enumerate(U):
    plt.imsave(f"./anim/{files[i].split(".")[0]}.png",u,  vmax=vmax, vmin=0, cmap='jet')
