import numpy as np
import cv2
import matplotlib.pyplot as plt

plt.rc("font",**{"family":"serif","serif":["Palatino"]})
plt.rc('text', usetex=True)

figure = plt.figure(figsize=(10, 10), tight_layout=True)

ax1 = figure.add_axes([0.06, 0.025, 0.92, 0.5])
ax2 = figure.add_axes([0.08, 0.55, 0.4, 0.43])
ax3 = figure.add_axes([0.56, 0.55, 0.4, 0.43])

ax1.set_xlabel(r"$X$[m]")
ax1.set_ylabel(r"$Y$[m]")
ax2.set_xlabel(r"$V_x$[m$\cdot$s$^{-1}$]")
ax2.set_ylabel(r"$V_y$[m$\cdot$s$^{-1}$]")
ax3.set_xlabel(r"$V_x$[m$\cdot$s$^{-1}$]")
ax3.set_ylabel(r"$V_y$[m$\cdot$s$^{-1}$]")

ax2.axis('equal')
ax3.axis('equal')

cap = cv2.VideoCapture('rot.mp4')

fr = 0
namedata = ["./Ux/%05d.csv","./Uy/%05d.csv"]

x1 = []
y1 = []
x2 = []
y2 = []

def load(x1,y1,x2,y2,fr):
    U  = np.loadtxt(namedata[0]%fr)
    V  = np.loadtxt(namedata[1]%fr)
    x1.append(U[coord[0]])
    y1.append(V[coord[0]])
    x2.append(U[coord[1]])
    y2.append(V[coord[1]])

ret, frame = cap.read()
rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

L = rgb.shape
l = [2/L[0],4/L[1]]
coord = [(int(1/l[0]),int(1/l[1])),(int(1/l[0]),int(3/l[1]))]
extent = [0, 4, 0, 2]
img = ax1.imshow(rgb, extent=extent)
load(x1,y1,x2,y2,fr)

line1 = ax2.plot(x1,y1,'r')[0]
line2 = ax3.plot(x2,y2,'b')[0]

ax1.plot([l[1]*coord[0][1]],[l[0]*coord[0][0]], 'ro')
ax1.plot([l[1]*coord[1][1]],[l[0]*coord[1][0]], 'bo')

while(cap.isOpened()):
    
    plt.savefig("./anim2/%05d.png"%fr)
    
    fr += 1
    
    load(x1,y1,x2,y2,fr)
    line1.set_xdata(x1)
    line1.set_ydata(y1)
    line2.set_xdata(x2)
    line2.set_ydata(y2)
    ax2.relim()
    ax2.autoscale_view()
    ax3.relim()
    ax3.autoscale_view()
    img.set_array(rgb)
    
    if (fr%2 == 0) or (fr%7 == 0):
        x1 = x1[1:]
        y1 = y1[1:]
        x2 = x2[1:]
        y2 = y2[1:]
    
    ret, frame = cap.read()
    rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

plt.savefig("./anim2/%05d.png"%fr)

img.set_array(rgb)

cap.release()
cv2.destroyAllWindows()
