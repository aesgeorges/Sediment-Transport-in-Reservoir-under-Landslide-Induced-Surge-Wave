import numpy as np
from setup import *
from advance import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 10
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
length = np.arange(L,0,-dx)
for i in range(N):
    surf[:,i] = -0.0000005 + (0.000001/L)*length[i]


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

total_time = 500

ax.plot_surface(X,Y,surf,cmap='cool')
ax.set(xlabel='X axis label', ylabel='Y axis label')
plt.show()
def animate(i, surf, U, V, Gu, Gv):
    ax.clear()
    Gu,Gv = update_matrices(U,V,f,nu_h,Gu,Gv,deltaZ,N,dt,dx,dy,dz,nu_v)
    surf, U, V = advance(surf,N,U,V,deltaZ,Gu,Gv,A,dt,dx,dy)
    #ax.set_xlim(0, 15)
    #ax.set_ylim(0, 15)
    ax.set(xlabel='X axis label', ylabel='Y axis label')

    #ax.set_zlim(-1e-13, 1e-13)
    ax.plot_surface(X,Y,surf,rstride=1, cstride=1, cmap='cool', linewidth=0, antialiased=False)
ani = animation.FuncAnimation(fig, animate, frames = total_time, fargs=(surf, U, V, Gu, Gv), interval = dt)
ani.save('Visualization/testd3.gif', writer='imagemagick', fps=24)
#plt.show()
