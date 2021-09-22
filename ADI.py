import matplotlib as mt
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy.ma.core import shape, sqrt 

#input 
length = float(10)
width = float(10)
rows = int(10) 
cols = int(10)
iterations = 500
alpha = float(1)
dt = float(0.001)
#boundary temperatures
rb = 0
lb = 1
tb = 1
bb = 0
#calculation values
r = [(i,cols - 1) for i in range(rows)]
l = [(i,0) for i in range(rows)]
t = [(0,i) for i in range(cols)]
b = [(rows - 1,i) for i in range(cols)]
dx = float(length / cols)
dy = float(width / rows)

#weights
a = float((dt*alpha)/(2*dx*dy))
ar = a*dy/dx
al = a*dy/dx
ab = a*dx/dy
at = a*dx/dy

#initialization
temp = np.zeros((rows,cols),dtype=float)
copy = temp
iter = 0
z = []

while iter < iterations:
    for n in range(rows*cols):
        r = int(n/rows)
        c = int(n%rows)
        cf = np.zeros((cols,cols),dtype=float)
        ct = np.zeros((cols),dtype=float)
        for i in range(cols): #horizontal
            apn = 1 + al + ar
            apo = 1 - at - ab
            if i == 0:
                #left
                if r == 0:
                    #top
                    ct[i] = (apo - at)*copy[r,i] + ab*copy[r+1,i] + 2*tb*at
                elif r == rows - 1:
                    #bottom
                    ct[i] = (apo - ab)*copy[r,i] + at*copy[r-1,i] + 2*bb*ab
                else:
                    #mid
                    ct[i] = apo*copy[r,i] + at*copy[r-1,i] + ab*copy[r+1,i]
                cf[i,i+1] = -ar
                cf[i,i] = apn + al
                ct[i] += 2*al*lb
            elif i == cols - 1:
                #right
                if r == 0:
                    #top
                    ct[i] = (apo - at)*copy[r,i] + ab*copy[r+1,i] + 2*tb*at
                elif r == rows - 1:
                    #bottom
                    ct[i] = (apo - ab)*copy[r,i] + at*copy[r-1,i] + 2*bb*ab
                else:
                    #mid
                    ct[i] = apo*copy[r,i] + at*copy[r-1,i] + ab*copy[r+1,i]
                cf[i,i] = apn + ar
                cf[i,i-1] = -al
                ct[i] += 2*ar*rb
            else:
                #midcells
                if r == 0:
                    #top
                    ct[i] = (apo - at)*copy[r,i] + ab*copy[r+1,i] + 2*tb*at
                elif r == rows - 1:
                    #bottom
                    ct[i] = (apo - ab)*copy[r,i] + at*copy[r-1,i] + 2*bb*ab
                else:
                    #mid
                    ct[i] = apo*copy[r,i] + at*copy[r-1,i] + ab*copy[r+1,i]
                cf[i,i] = apn
                cf[i,i+1] = -ar
                cf[i,i-1] = -al

        inv_cf = np.linalg.inv(cf)
        rnew = np.dot(inv_cf,ct)
        for i in range(cols):
            temp[r,i] = rnew[i]
        
        for i in range(rows): #vertical
            apn = 1 + at + ab
            apo = 1 - al - ar
            if i == 0:
                #top
                if c == 0:
                    #left
                    ct[i] = (apo - al)*copy[i,c] + ar*copy[i,c+1] + 2*lb*al
                elif c == cols - 1:
                    #right
                    ct[i] = (apo - ar)*copy[i,c] + al*copy[i,c-1] + 2*rb*ar
                else:
                    #mid
                    ct[i] = apo*copy[i,c] + al*copy[i,c-1] + ar*copy[i,c+1]
                cf[i+1,i] = -ab
                cf[i,i] = apn + at
                ct[i] += 2*at*tb
            elif i == cols - 1:
                #bottom
                if c == 0:
                    #left
                    ct[i] = (apo - al)*copy[i,c] + ar*copy[i,c+1] + 2*lb*al
                elif c == cols - 1:
                    #right
                    ct[i] = (apo - ar)*copy[i,c] + al*copy[i,c-1] + 2*rb*ar
                else:
                    #mid
                    ct[i] = apo*copy[i,c] + al*copy[i,c-1] + ar*copy[i,c+1]
                cf[i,i] = apn + ab
                cf[i-1,i] = -at
                ct[i] += 2*ab*bb
            else:
                #midcells
                if c == 0:
                    #left
                    ct[i] = (apo - al)*copy[i,c] + ar*copy[i,c+1] + 2*lb*al
                elif c == cols - 1:
                    #right
                    ct[i] = (apo - ar)*copy[i,c] + al*copy[i,c-1] + 2*rb*ar
                else:
                    #mid
                    ct[i] = apo*copy[i,c] + al*copy[i,c-1] + ar*copy[i,c+1]
                cf[i,i] = apn
                cf[i+1,i] = -ab
                cf[i-1,i] = -at

        inv_cf = np.linalg.inv(cf)
        rnew = np.dot(inv_cf,ct)
        for i in range(cols):
            temp[i,c] = rnew[i]
    copy = temp
    iter += 1
    z.append(temp)
xcells = np.array([(2*i+1)*dx/2 for i in range(rows)])
ycells = np.array([(2*i+1)*dy/2 for i in range(cols)])
plt.contourf(ycells, xcells,temp,30,cmap=plt.cm.jet)
plt.gca().invert_yaxis()
plt.show()





        
