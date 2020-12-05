import numpy as np
import cv2
import matplotlib.pyplot as plt

plt.rc("font",**{"family":"serif","serif":["Palatino"]})
plt.rc('text', usetex=True)

figure = plt.figure(figsize=(10, 10), tight_layout=True)

ax1 = figure.add_axes([0.06, 0.025, 0.92, 0.5])
ax2 = figure.add_axes([0.06, 0.55, 0.3, 0.43])
ax3 = figure.add_axes([0.43, 0.55, 0.55, 0.215])
ax4 = figure.add_axes([0.43, 0.765, 0.55, 0.215], sharex=ax3)

ax1.set_xlabel(r"$X$[m]")
ax1.set_ylabel(r"$Y$[m]")
ax2.set_xlabel(r"$V_x$[m$\cdot$s$^{-1}$]")
ax2.set_ylabel(r"$Y$[m]")
ax3.set_xlabel(r"$t$[s]")
ax3.set_ylabel(r"$V$[m$^{3}$]")
ax4.set_ylabel(r"$\Delta V$[m$^{3}$]")
#ax3.set_xscale("log")
#ax3.set_yscale("symlog")
#ax4.set_yscale("symlog")
plt.setp(ax4.get_xticklabels(), visible=False)

cap = cv2.VideoCapture('rot.mp4')

fr = 0
namedata = ["./Ux/%05d.csv","./Uy/%05d.csv"]

x1 = []
y1 = []
x2 = []
y2 = []

Poiseuille = lambda Dv, h: Dv*6/h
speedth = lambda Vmax, y, h: Vmax*(y-(y*y)/h)/h

ret, frame = cap.read()
rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

L = rgb.shape
l = [2/L[0],4/L[1]]
xcoord = [int(.75/l[1]),int(2.5/l[1])]

constr = lambda : rgb[:,j].sum(1).nonzero()[0]
j = xcoord[0]
y1coord = constr()
j = xcoord[1]
y2coord = constr()

y1 = l[0]*y1coord
y2 = l[0]*y2coord

extent = [0, 4, 0, 2]
img = ax1.imshow(rgb, extent=extent)

yth = y1[rgb[y1coord,xcoord[0]].sum(1)>3]
h1 = yth[-1] - yth[0]

yth = y2[rgb[y2coord,xcoord[1]].sum(1)>3]
h2 = yth[-1] - yth[0]
def load(fr):
    U  = np.loadtxt(namedata[0]%fr)
    V  = np.loadtxt(namedata[1]%fr)
    Vx1 = U[y1coord,xcoord[0]]*(rgb[y1coord,xcoord[0]].sum(1)>30)
    Vx2 = U[y2coord,xcoord[1]]*(rgb[y2coord,xcoord[1]].sum(1)>30)
    
    return(Vx1[::-1], Vx2[::-1])

Vx1, Vx2 = load(fr)

line1 = ax2.plot(Vx1,y1,'r')[0]
line2 = ax2.plot(Vx2,y2,'b')[0]

ax1.plot(l[0]*xcoord[0]+y1*0,y1, 'r-')
ax1.plot(l[0]*xcoord[1]+y2*0,y2, 'b-')

dt = 5e-3
V1 = [Vx1.sum()*l[0]]
V2 = [Vx2.sum()*l[0]]
DV = [V2[-1]-V1[-1]]
t = [fr*dt]
line3 = ax3.plot(V1,t,'r')[0]
line4 = ax3.plot(V2,t,'b')[0]

line5 = ax4.plot(DV,t,'#aa00aa')[0]

vmax = Poiseuille(V1[-1],h2)
Vxth = speedth(vmax, yth-yth[0], h2)
line6 = ax2.plot(Vxth, yth, 'k--', label="Écoulement théorique")[0]
ax2.legend()
while(cap.isOpened()):
    
    plt.savefig("./anim2/%05d.svg"%fr)
    
    fr += 1
    
    Vx1, Vx2 = load(fr)
    line1.set_xdata(Vx1)
    line2.set_xdata(Vx2)
    vmax = Poiseuille(V1[-1],h2)
    Vxth = speedth(vmax, yth-yth[0], h2)
    line6.set_xdata(Vxth)
    V1.append(Vx1.sum()*l[0])
    V2.append(Vx2.sum()*l[0])
    t.append(fr*dt)
    line3.set_ydata(V1)
    line3.set_xdata(t)
    line4.set_ydata(V2)
    line4.set_xdata(t)
    
    DV.append(V2[-1]-V1[-1])
    line5.set_ydata(DV)
    line5.set_xdata(t)
    
    ax2.relim()
    ax2.autoscale_view()
    ax3.relim()
    ax3.autoscale_view()
    ax4.relim()
    ax4.autoscale_view()
    
    img.set_array(rgb)
    
    ret, frame = cap.read()
    rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

plt.savefig("./anim2/%05d.svg"%fr)

img.set_array(rgb)

cap.release()
cv2.destroyAllWindows()
