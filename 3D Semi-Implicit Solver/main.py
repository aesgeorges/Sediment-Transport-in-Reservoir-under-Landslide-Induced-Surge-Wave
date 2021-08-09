import numpy as np
from setup import *
from advance import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 20
L = 5
W = 5
H = 5
dx = L/N
dy = W/N
dz = H/N
dt = 60
nu_h = 0.0001
nu_v = 0.0001
f = 0.0000891 # Coriolis Parameter - at San Francisco latitude 37.773972

U = np.zeros((N,N,N))
V = np.zeros((N,N,N))
#W = np.zeros((N,N,N))
Gu = np.zeros((N,N,N))
Gv = np.zeros((N,N,N))
deltaZ = dz*np.ones(N)
Cz = 50

surf = np.zeros((N,N)) # Surface - eta in x,z
length = np.arange(0,L,dx)
for i in range(N):
    surf[i,:] = length[i]*0.5 +0.0001

X = np.arange(0, W, dx)
Y = np.arange(0, L, dy)
X, Y = np.meshgrid(X, Y)
fig = plt.figure()
ax = plt.axes(projection='3d')

Fu, Fv = F_op(U, V, f, nu_h, dt, dx, dy, dz, N)
U,V,deltaZ,Gu,Gv,A = set_matrices(U,V,Cz,Fu,Fv,Gu,Gv,deltaZ,N,dt,dz,nu_v)
#surf, U, V = advance(surf,N,U,V,deltaZ,Gu,Gv,A,dt,dx,dy)
#print(surf)
#print(U)
#print(V)

total_time = 900

def animate(i, surf, U, V, Gu, Gv):
    ax.clear()
    Gu,Gv = update_matrices(U,V,f,nu_h,Gu,Gv,deltaZ,N,dt,dx,dy,dz,nu_v)
    surf, U, V = advance(surf,N,U,V,deltaZ,Gu,Gv,A,dt,dx,dy)
    ax.plot_surface(X,Y,surf,rstride=1, cstride=1, cmap='cool', linewidth=0, antialiased=False)

ani = animation.FuncAnimation(fig, animate, frames = total_time, fargs=(surf, U, V, Gu, Gv), interval = dt)
ani.save('Visualization/semiimplicittest0.gif', writer='imagemagick', fps=24)
#plt.show()
