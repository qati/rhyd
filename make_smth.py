#!/bin/python

import numpy as np
import scipy.signal as sig
import sys

FILE=str(sys.argv[1])
T0=float(sys.argv[2])


def get_line(file):
    NUM_OF_ELEMENTS_IN_ROW=6
    s = file.readline();
    if s=="\n":
        return ['q']
    l = s.split(' ', NUM_OF_ELEMENTS_IN_ROW)
    l[NUM_OF_ELEMENTS_IN_ROW-1]=l[NUM_OF_ELEMENTS_IN_ROW-1].replace('\n','')
    #ln = np.array([],dtype='float64')
    ln=[]
    for i in l:
        #ln=np.append(ln, [float(i)])
        ln.append(float(i))
    return ln

def get_d(file_name):
    NUM_OF_ELEMENTS_IN_ROW=6
    f = open(file_name, "r")
    
    x_min = 0
    y_min = 0
    dy    = 0

    first=[]
    ln  = get_line(f)
    lnb = ln
    while ln!=['q']:
        if dy < 0:
            if ln[1]==0:
                dy = abs(dy)
            else:
                dy += abs(ln[1])
                dy  = abs(dy)
        
        if dy==0:
            x_min = ln[0]
            y_min = ln[1]
            dy    = -abs(ln[1])
    
        first.append([ln[i] for i in range(0, NUM_OF_ELEMENTS_IN_ROW)])
        lnb = ln
        ln  = get_line(f)

    if lnb[1]!=abs(y_min):
        print("Error: y1!=y2")
        sys.exit(0)

    N = 2*int(np.round(abs(x_min)/dy))+1

    first=np.array(first, dtype='float64')    
    if first.shape[0]!=N:
        print("Error! first.shape[0](=",first.shape[0],")!=N(=",N,")")
        sys.exit(0)
    
    data_n0 = np.zeros((N,N), dtype='float64')
    data_p0 = data_n0.copy()
    data_vx = data_n0.copy()
    data_vy = data_n0.copy()
    for i in range(0, N):
        data_n0[0][i] = first[i][2]
        data_vx[0][i] = first[i][3]
        data_vy[0][i] = first[i][4]
        data_p0[0][i] = first[i][5]
    

    for i in range(1,N):
        for j in range(0, N):
            ln = get_line(f)
            if ln != ['q']:
                data_n0[i][j] = ln[2]
                data_vx[i][j] = ln[3]
                data_vy[i][j] = ln[4]
                data_p0[i][j] = ln[5]
            else:
                ln = get_line(f)
                if ln == ['q']:
                    print("Something wrong: ln==['q'] => 2 blank lines!!!")
                    sys.exit(1)
                data_n0[i][j] = ln[2]
                data_vx[i][j] = ln[3]
                data_vy[i][j] = ln[4]
                data_p0[i][j] = ln[5]
            
    print("x_min=%s, y_min=%s, dy=%s, N=%d"%(x_min, y_min, dy, N))
    f.close()
    return x_min, y_min, data_n0, data_vx, data_vy, data_p0

def print_d(file_name, m, N0, VX, VY, P0, NORM_N0, NORM_P0, NORM_VX, NORM_VY):
    f = open(file_name, "w")
    x = np.linspace(m[0], abs(m[0]), N0.shape[0])
    y = np.linspace(m[1], abs(m[1]), N0.shape[0])
    x,y = np.meshgrid(x,y)
    for i in range(0, N0.shape[0]):
        for j in range(0, N0.shape[0]):
            N0[i][j] *= NORM_N0
            P0[i][j] *= NORM_P0
            VX[i][j] *= NORM_VX
            VY[i][j] *= NORM_VY
            s=str(y[i][j])+" "+str(x[i][j])+" "+str(N0[i][j])+" "+str(VX[i][j])+" "+str(VY[i][j])+" "+str(P0[i][j])+"\n"
            f.write(s)
        f.write("\n")
    f.close()
    return

x_min, y_min, n0, vx, vy, p0 = get_d(FILE)

x = np.linspace(x_min, abs(x_min), n0.shape[0])
y = np.linspace(y_min, abs(y_min), n0.shape[0])

X,Y = np.meshgrid(x,y)

z1 = np.zeros(X.shape)
z2 = np.zeros(Y.shape)
for i in range(0,X.shape[0]):
    for j in range(0, X.shape[1]):
        if abs(X[i][j])<1:
            z1[i][j] = np.exp(-1/(1-X[i][j]**2))
        if abs(Y[i][j])<1:
            z2[i][j] = np.exp(-1/(1-Y[i][j]**2))
Z = z1*z2
i
smth_n0 = sig.convolve2d(n0,Z, mode='same')
smth_p0 = sig.convolve2d(p0,Z, mode='same')

print("Time: ", T0);
for i in range(vx.shape[0]):
    for j in range(vx.shape[1]):
        if n0[i][j]==0 and vx[i][j]==0 and smth_n0[i][j]!=0:
            vx[i][j]=Y[i][j]/T0
        if n0[i][j]==0 and vy[i][j]==0 and smth_n0[i][j]!=0:
            vy[i][j]=X[i][j]/T0
smth_vx = vx
smth_vy = vy
#smth_vx = sig.convolve2d(vx,Z, mode='same')
#smth_vy = sig.convolve2d(vy,Z, mode='same')

n0_max = n0[0][0]
p0_max = p0[0][0]
vx_max = vx[0][0]
vy_max = vy[0][0]
smth_max1 = smth_n0[0][0]
smth_max2 = smth_p0[0][0]
smth_max3 = smth_vx[0][0]
smth_max4 = smth_vy[0][0]
for i in range(0, n0.shape[0]):
    for j in range(0, n0.shape[1]):
        if n0_max<n0[i][j]:
            n0_max = n0[i][j]
        if p0_max<p0[i][j]:
            p0_max = p0[i][j]
        if vx_max<vx[i][j]:
            vx_max = vx[i][j]
        if vy_max<vy[i][j]:
            vy_max = vy[i][j]
        if smth_max1<smth_n0[i][j]:
            smth_max1 = smth_n0[i][j]
        if smth_max2<smth_p0[i][j]:
            smth_max2 = smth_p0[i][j]
        if smth_max3<smth_vx[i][j]:
            smth_max3 = smth_vx[i][j]
        if smth_max4<smth_vy[i][j]:
            smth_max4 = smth_vy[i][j]

print_d(FILE, (x_min, y_min), smth_n0, smth_vx, smth_vy, smth_p0, n0_max/smth_max1, p0_max/smth_max2, vx_max/smth_max3, vy_max/smth_max4)
