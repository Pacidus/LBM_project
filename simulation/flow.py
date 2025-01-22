import numpy as np
import cv2
import matplotlib.pyplot as plt

plt.rc("font",**{"family":"serif","serif":["Palatino"]})
plt.rc('text', usetex=True)

figure = plt.figure(figsize=(10, 10), tight_layout=True)

ax1 = figure.add_axes([0.06, 0.025, 0.92, 0.5])

ax1.set_xlabel(r"$X$[m]")
ax1.set_ylabel(r"$Y$[m]")

cap = cv2.VideoCapture('rot.mp4')

fr = 0
namedata = ["./Ux/%05d.csv","./Uy/%05d.csv"]

def load(x1,y1,x2,y2,fr):
    U  = np.loadtxt(namedata[0]%fr)
    V  = np.loadtxt(namedata[1]%fr)
    return(U,V)

ret, frame = cap.read()
rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

L = rgb.shape
l = [2/L[0],2.5/L[1]]
coord = [(int(1/l[0]),int(1/l[1])),(int(1/l[0]),int(3/l[1]))]

extent = [0, 2.5, 0, 2]

img = ax1.imshow(rgb, extent=extent)
Ux, Uy = load(fr)
x = [0]*10
y = [i*l[0] for i in range(10)]
line1 = ax1.plot(x,y,"k.")[0]
while(cap.isOpened()):
    
    plt.savefig("./anim2/%05d.png"%fr)
    
    fr += 1
    
    Ux,Vx = load(fr)
    line1.set_xdata(x)
    line1.set_ydata(y)
    img.set_array(rgb)
    
    if(fr%10 == 0):
        x = x1[1:]
        y = y1[1:]
    
    ret, frame = cap.read()
    rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

plt.savefig("./anim2/%05d.png"%fr)

img.set_array(rgb)

cap.release()
cv2.destroyAllWindows()
