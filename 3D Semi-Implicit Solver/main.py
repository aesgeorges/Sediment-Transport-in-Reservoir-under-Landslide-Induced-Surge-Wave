import numpy as np
from setup import *
from advance import *

N = 10
L = 20
W = 20
H = 20
dx = L/N
dy = W/N
dz = H/N
dt = 60
nu_h = 0.0001
nu_v = 0.0001
f = 0.0000891 # Coriolis Parameter - at San Francisco latitude 37.773972

U = np.zeros((N,N,N))
V = np.zeros((N,N,N))
W = np.zeros((N,N,N))
Gu = np.zeros((N,N,N))
Gv = np.zeros((N,N,N))
deltaZ = dz*np.ones(N)
Cz = 50

surf = np.zeros((N,N)) # Surface - eta in x,z
length = np.arange(0,L,dx)
for i in range(N):
    surf[i,:] = length[i]*0.5 + 15

Fu, Fv = F_op(U, V, f, nu_h, dt, dx, dy, N)
U,V,deltaZ,Gu,Gv,A = set_matrices(U,V,Cz,Fu,Fv,Gu,Gv,deltaZ,N,dt,dz,nu_v)
surf, U, V = advance(surf,N,U,V,deltaZ,Gu,Gv,A,dt,dx,dy)
print(surf)
#print(U)
#print(V)

