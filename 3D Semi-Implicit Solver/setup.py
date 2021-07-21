import numpy as np
import math
from scipy.sparse import diags

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
def F_op(u, v, f, nu, dt, dx, dy, N):
    Fu = np.zeros((N, N, N))
    Fv = np.zeros((N, N, N))
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
                Fv[i,j,k] = v_interp + nu*(dt/(dy**2))*(v_interpP - 2*v_interp + v_interpM) + (1/(dy**2))*(v_interp_jP - 2*v_interp + v_interp_jM) + f*dt*u_interp

    return Fu, Fv

# Returns built tridiagonal matrix A
def tridiag_maker(a,b,c,N):
    data = [a.tolist(), b.tolist(), c.tolist()]
    pos = [0,1,-1]
    tridiag = (diags(data, pos, (N,N))).toarray()
    return tridiag


# Setting up all 3D arrays/matrices
def set_matrices(U,V,Cz,Fu,Fv,Gu,Gv,deltaZ,N,dt,dz,nu_v):
    print('setting up matrices and 3D arrays...')
    str_x = 5
    str_y = 5
    bot_str = ((9.81*dt)/Cz**2)*(math.sqrt(U[0,0,-1]**2 + V[0,0,-1]**2))
    # Setting up G arrays
    for i in range(N):
        for j in range(N):
            Gu[i,j] = deltaZ*Fu[i,j]
            Gv[i,j] = deltaZ*Fv[i,j]
    Gu[:,:,0] = Gu[:,:,0] + dt*str_x
    Gv[:,:,0] = Gv[:,:,0] + dt*str_y 
    # Setting up A matrix
    a = (dz + nu_v*(dt/dz))*np.ones(N)
    b = (-nu_v*(dt/dz))*np.ones(N-1)
    c = b
    for i in range(1,N-1):
        a[i] = a[i] + nu_v*(dt/dz)
    a[-1] = a[-1] + bot_str
    A = tridiag_maker(a,b,c,N)
    return U,V,deltaZ,Gu,Gv,A