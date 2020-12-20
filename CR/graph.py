#!/usr/bin/env python3
# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#plt.style.use('dark_background')

rc("font",**{"family":"serif","serif":["Palatino"]})
rc('text', usetex=True)

#r = [-1.5,1.5]

#plt.axis("square")
##plt.tight_layout()
#plt.xlim(*r)
#plt.ylim(*r)
#plt.xlabel(r"$\delta_x$")
#plt.ylabel(r"$\delta_y$")
#plt.xticks([-1, 0, 1])
#plt.yticks([-1, 0, 1])
#plt.plot(r , [0.5, 0.5], "w")
#plt.plot(r , [-0.5, -0.5], "w")
#plt.plot([0.5, 0.5], r, "w")
#plt.plot([-0.5, -0.5], r, "")
##plt.plot([1,-1,0,0,1,-1,-1,1],[0,0,1,-1,1,1,-1,-1],"ko")

#c = "#ff004a"
#l = .9
#plt.plot([0],[0],color=c, marker="o")
#plt.arrow(0, 0, l, 0, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, 0, l, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, -l, 0, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, 0, -l, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, l, l, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, -l, l, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, -l, -l, head_width=0.05, head_length=0.1, fc=c, ec=c)
#plt.arrow(0, 0, l, -l, head_width=0.05, head_length=0.1, fc=c, ec=c)

#plt.text(0, 0, r"$\vec{e}_1$",ha="center", va="center", size=20)
#plt.text(1, 0, r"$\vec{e}_2$",ha="center", va="center", size=20)
#plt.text(0, 1, r"$\vec{e}_3$",ha="center", va="center", size=20)
#plt.text(-1, 0, r"$\vec{e}_4$",ha="center", va="center", size=20)
#plt.text(0, -1, r"$\vec{e}_5$",ha="center", va="center", size=20)
#plt.text(1, 1, r"$\vec{e}_6$",ha="center", va="center", size=20)
#plt.text(-1, 1, r"$\vec{e}_7$",ha="center", va="center", size=20)
#plt.text(-1, -1, r"$\vec{e}_8$",ha="center", va="center", size=20)
#plt.text(1, -1, r"$\vec{e}_9$",ha="center", va="center", size=20)

#plt.savefig("ei.svg",transparent=True)
#plt.cla()
#plt.clf()

#img = plt.imread("otter.png")
#shape = img.shape
#shape = (shape[0]*2, shape[1], shape[2])
#Bimg = np.zeros(shape)
#Bimg[:shape[0]//2] = img
#Bimg[shape[0]//2:] = img
#Bimg = np.roll(Bimg, shape[0]//4, 0)
#plt.plot(np.linspace(0,shape[1]-1,2), [shape[0]//4]*2, 'k')
#plt.plot(np.linspace(0,shape[1]-1,2), [shape[0] - (shape[0]//4)]*2, 'k')
#plt.imshow(Bimg)
#plt.yticks([(shape[0]//8 + 1)*i for i in range(8)],[((shape[0]//8)*i + shape[0]//4 +1)%(shape[0]//2 - 1) for i in range(8)])
#plt.xlabel(r"$\delta_x$")
#plt.ylabel(r"$\delta_y$")
#plt.tight_layout()
#plt.savefig("periodic.pdf",transparent=True)

img = plt.imread("Queen.png")
shape = img.shape
plt.plot([1]*2, [0,shape[0]], '#550', linewidth=3, label="Speed inlet")
plt.plot([shape[1]-2]*2, [0,shape[0]],'#606', linewidth=3, label="Pressure outlet")
plt.imshow(img)
plt.xlabel(r"$\delta_x$")
plt.ylabel(r"$\delta_y$")
plt.legend()
plt.tight_layout()
plt.savefig("iolet.pdf",transparent=True)

