import numpy as np
import math
from numpy.lib.shape_base import apply_along_axis
import matplotlib.animation as animation
from scipy import sparse
from surf import surf_solver
from scipy.sparse import diags

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Mesh setup
N = 16
L = 20
W = 20
H = 20
dx = L/N
dy = W/N
dz = H/N
dt = 60
nu = 0.0001 # Horizontal Eddy Viscosity
v_nu = 0.0001*np.ones(N) # Vertical Eddy Viscosity
f = 0.0000891 # Coriolis Parameter - at San Francisco latitude 37.773972
g = 9.81



u = np.zeros((N, N, N)) # Velocity in x
v = np.zeros((N, N, N))# Velocity in y
w = np.zeros((N, N, N)) # Velocity in z
surf = np.zeros((N,N)) # Surface - eta in x,z
Fu = np.zeros((N, N, N))
Fv = np.zeros((N, N, N))

m = np.zeros((N,N)) # k-index of bottom
M = (N/2 )*np.ones((N,N)) # k-index of top

length = np.arange(0,L,dx)
for i in range(N):
    surf[i,:] = length[i]*0.5 + 15
print(surf)
#### Solve for F operator

# Interpolator
def interpolator(x, i, j, k):
    u_ini = 0.4
    v_ini = 0
    w_ini = 0
    a = 0.5 #u_ini*(dt/dx)
    b = 0.5 #v_ini*(dt/dy)
    c = 0.5 #w_ini*(dt/dz)
    l = int(math.modf(a)[1])
    m = int(math.modf(b)[1])
    n = int(math.modf(c)[1])
    p = int(math.modf(a)[0])
    q = int(math.modf(b)[0])
    r = int(math.modf(c)[0])
    c_interp = (1-r)*((1-p)*((1-q)*x[i-l,j-m,k-n] + q*x[i-l,j-m-1,k-n]) + p*((1-q)*x[i-l-1,j-m,k-n] + q*x[i-l-1,j-m-1,k-n])) + r*((1-p)*((1-q)*(x[i-l,j-m,k-n-1] + q*x[i-l,j-m-1,k-n-1])) + p*((1-q)*x[i-l-1,j-m,k-n-1] + q*x[i-l-1,j-m-1,k-n-1]))
    return c_interp

# Solve F operator
def F_op(u, v):
    print("Computing F operator")
    for i in range(1, N-1):
        for j in range(1, N-1):
            for k in range(N):
                u_interp = interpolator(u, i, j, k)
                v_interp = interpolator(v, i, j, k)

                u_interpP = interpolator(u, i+1, j, k)
                u_interpM = interpolator(u, i-1, j, k)
                u_interp_jP = interpolator(u, i, j+1, k)
                u_interp_jM = interpolator(u, i, j-1, k)
                Fu[i,j,k] = u_interp + nu*(dt/(dx**2))*(u_interpP - 2*u_interp + u_interpM) + (1/(dy**2))*(u_interp_jP - 2*u_interp + u_interp_jM) + f*dt*v_interp

                v_interpP = interpolator(v, i+1, j, k)
                v_interpM = interpolator(v, i-1, j, k)
                v_interp_jP = interpolator(v, i, j+1, k)
                v_interp_jM = interpolator(v, i, j-1, k)
                Fv[i,j,k] = v_interp + nu*(dt/(dx**2))*(v_interpP - 2*v_interp + v_interpM) + (1/(dy**2))*(v_interp_jP - 2*v_interp + v_interp_jM) + f*dt*u_interp

    return Fu, Fv

def tridiag_maker(a,b,c,N):
    data = [a.tolist(), b.tolist(), c.tolist()]
    pos = [0,1,-1]
    tridiag = diags(data, pos, (N,N))
    return tridiag

def set_matrices(N):
    print('setting up matrices and vectors...')
    U = np.zeros((N,N,N))
    V = np.zeros((N,N,N))
    stressx = np.zeros((N,N))
    stressy = np.zeros((N,N))
    deltaZ = dz*np.ones(N)
    Gi = np.zeros((N,N,N))
    Gj = np.zeros((N,N,N))
    A = np.zeros((N,N))
    for i in range(N-1):
        for j in range(N-1):
            # Setting Boundary stress conditions
            stressx = (u[i,j,int(M[i+1,j+1])]-u[i,j,int(M[i,j])])/(dz)
            stressy = (v[i,j,int(M[i+1,j+1])]-v[i,j,int(M[i,j])])/(dz)

    for i in range(N):
        for j in range(N-1):
            for k in range(int(m[i,j]),int(M[i,j])):
                U[i,j,k] = u[i,j,k]
                V[i,j,k] = v[i,j,k]
                Gi[i,j,k] = dz*Fu[i,j,k]
                Gj[i,j,k] = dz*Fv[i,j,k]
    Gi[:,:,0] = Gi[:,:,0] + dt*stressx
    Gj[:,:,0] = Gj[:,:,0] + dt*stressy
    #Gi = Gi[:,None]
    #Gj = Gj[:,None]
    a = dz + v_nu*(dt/dz)
    b = -v_nu[:-1]*(dt/dz)
    c = -v_nu[:-1]*(dt/dz)
    for i in range(1,N):
        a[i] = dz +  v_nu[i]*(dt/dz) + v_nu[i-1]*(dt/dz)*(-v_nu[i-1])*(dt/dz)
    A = tridiag_maker(a,b,c,N)
    return U,V,deltaZ,Gi,Gj,A

def solve_u(newsurf, deltaZ,Gi,A,g,dt,dx):
    U = np.zeros((N,N,N))
    for i in range(N-1):
        for j in range(N):
            #U[i,j,:] = ((Gi[i,j] - g*(dt/dx)*(newsurf[i+1,j]-newsurf[i,j])*deltaZ)/A[i,j])
            b = Gi[i,j] - (g*(dt/dx)*(newsurf[i+1,j]-newsurf[i,j])*deltaZ).reshape(N)
            U[i,j] = np.linalg.lstsq(A, b, rcond=None)[0]
                #bloc = (v_nu[k+1]/dz)*(u[i,j,k+1] - u[i,j,k]) - (v_nu[k]/dz)*(u[i,j,k] - u[i,j,k-1])
                #u[i,j,k] = Fu[i,j,k] - g*(dt/dx)*(surf[i+1, j] - surf[i,j]) + (dt/dz)*bloc
    return U

def solve_v(newsurf,deltaZ,Gj,A,g,dt,dy):
    V = np.zeros((N,N,N))
    for i in range(N-1):
        for j in range(N-1):
            b = Gj[i,j] - (g*(dt/dx)*(newsurf[i+1,j]-newsurf[i,j])*deltaZ).reshape(N)
            V[i,j] = np.linalg.lstsq(A, b, rcond=None)[0]
                #bloc = (v_nu[k+1]/dz)*(v[i,j,k+1] - v[i,j,k]) - (v_nu[k]/dz)*(v[i,j,k] - v[i,j,k-1])
                #v[i,j,k] = Fv[i,j,k] - g*(dt/dy)*(surf[i, j+1] - surf[i,j]) + (dt/dz)*bloc
    return V

# Boundary Conditions
#u[M] = 
Fu, Fv = F_op(u, v)
U,V,deltaZ,Gi,Gj,A = set_matrices(N)
deltaZ = deltaZ.reshape((N, 1))
A = A.todense()

X = np.arange(0, W, dx)
Y = np.arange(0, L, dy)
X, Y = np.meshgrid(X, Y)
fig = plt.figure()
ax = fig.gca(projection='3d')
#wav = ax.plot_surface(X,Y,surf,rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)


surf = surf_solver(surf, deltaZ,Gi,Gj,A,g,dt,dx,dy)
print(surf)
U = solve_u(surf,deltaZ,Gi,A,g,dt,dx)
V = solve_v(surf,deltaZ,Gj,A,g,dt,dy)

#wav = ax.plot_surface(X,Y,surf,rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)


def advance(surf):
    surf = surf_solver(surf, deltaZ,Gi,Gj,A,g,dt,dx,dy)
    print(surf)
    U = solve_u(surf,deltaZ,Gi,A,g,dt,dx)
    V = solve_v(surf,deltaZ,Gj,A,g,dt,dy)

total_time = 900

def animate(i):
    advance(surf)
    ax.plot_surface(X,Y,surf,rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)

ani = animation.FuncAnimation(fig, animate, frames = total_time, interval = dt)
plt.show()
    

 