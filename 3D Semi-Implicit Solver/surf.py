import numpy as np
import math
from numpy.linalg import inv

N = 80

s = np.zeros((N,N))
d = np.zeros((N,N))
q = np.zeros((N,N))
d = np.zeros((N,N))
b = np.zeros((N,N))

si = np.zeros(N)
sj = np.zeros(N)
ai = np.zeros(N)
aj = np.zeros(N)
e = np.zeros((N,N))
Mp = np.zeros((N,N))
r = np.zeros((N,N))


epsilon = 0.01
GUESS = 1.6

def surf_solver(surf, deltaZ,Gi,Gj,A, g, dt, dx, dy):
    print('solving surface...')
    print('computing compact terms')
    termZ = np.transpose(deltaZ)@inv(A)@deltaZ
    #termGi = np.transpose(deltaZ)*inv(A)*Gi
    #termGj = np.transpose(deltaZ)*inv(A)*Gj
    si[0] = g*((dt**2)/(dx**2))*termZ
    sj[0] = g*((dt**2)/(dy**2))*termZ
    for i in range(N-1):
        for j in range(N-1):
            termGi_pl = np.transpose(deltaZ)*inv(A)*((Gi[i+1,j,:]).reshape((N,1)))
            termGi_mn = np.transpose(deltaZ)*inv(A)*(Gi[i,j,:]).reshape((N,1))
            termGj_pl = np.transpose(deltaZ)*inv(A)*(Gj[i,j+1,:]).reshape((N,1))
            termGj_mn = np.transpose(deltaZ)*inv(A)*(Gj[i,j,:]).reshape((N,1))
            si[i] = g*((dt**2)/(dx**2))*termZ
            sj[j] = g*((dt**2)/(dy**2))*termZ
            d[i,j] = 1 + si[i] + si[i-1] + sj[j] +sj[j-1]
            #print(d)
            q[i,j] = surf[i,j] - (dt/dx)*(termGi_pl - termGi_mn) - (dt/dy)*(termGj_pl- termGj_mn)
            ai[i-1] = si[i-1]/(math.sqrt(d[i,j]*d[i-1,j]))
            ai[i] = si[i]/(math.sqrt(d[i,j]*d[i,j]))
            aj[j-1] = sj[j-1]/(math.sqrt(d[i,j]*d[i,j-1]))
            aj[j] = sj[j]/(math.sqrt(d[i,j]*d[i,j]))
            b[i,j] = q[i,j]/(math.sqrt(d[i,j]))

    p = np.zeros((N,N))
    print('setting up gradient method...')
    for i in range(N-1):
        for j in range(N-1):
            e[i,j] = GUESS
            r[i,j] = e[i,j] - ai[i]*e[i+1,j] - ai[i-1]*e[i-1,j] - aj[j-1]*e[i,j-1] - b[i,j]
            p[i,j] = r[i,j]
            Mp[i,j] = p[i,j] - ai[i]*p[i+1,j] - ai[i-1]*p[i-1,j] - aj[j]*p[i,j+1] - aj[j-1]*p[i,j-1]
    print('running conjugate method algorithm...')
    enew = np.zeros((N,N))
    while np.dot(r.all(),r.all()) > epsilon:
        alpha = np.dot(r,r)/np.dot(p,Mp)
        enew = e - alpha*p
        rnew = r - alpha*Mp
        beta = np.dot(rnew,rnew)/np.dot(r,r)
        p = r + beta*p
    print(enew)
    surf = enew*np.sqrt(d)
    return surf
