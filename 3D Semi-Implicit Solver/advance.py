import numpy as np
from numpy.linalg import inv, lstsq
from scipy.sparse import diags
import pentapy as pp
import math
from setup import F_op

def advance(surf,N,U,V,deltaZ,Gu,Gv,A,dt,dx,dy):
    surf = solve_surf(surf,N,deltaZ,Gu,Gv,A,dt,dx,dy)
    U = solve_U(surf,U,N,deltaZ,Gu,A,dt,dx)
    V = solve_V(surf,V,N,deltaZ,Gv,A,dt,dy)
    return surf, U, V

def fivediag_maker(a,b,c,d,e,N):
    data = [a, b, d, c, e]
    pos = [0,1,-1,2,-2]
    five_diag = (diags(data, pos, (N*N,N*N))).toarray()
    return five_diag

def solve_surf(surf,N,deltaZ,Gu,Gv,A,dt,dx,dy):
    d = np.zeros((N,N))
    q = np.zeros((N,N))
    q_vect = []
    d_vect = []
    surf_vect = []
    si_vect = []
    sj_vect = []
    #print('surface calcs...')
    #print('computing compact terms...')
    termZ = np.transpose(deltaZ)@inv(A)@deltaZ
    si = (9.81*((dt**2)/(dx**2))*termZ)*np.ones((N,N))
    sj = (9.81*((dt**2)/(dy**2))*termZ)*np.ones((N,N))
    for i in range(N-1):
        for j in range(N-1):
            termGu = np.transpose(deltaZ)@inv(A)@Gu[i+1,j]
            termGv = np.transpose(deltaZ)@inv(A)@Gv[i,j]
            termGu_mn = np.transpose(deltaZ)@inv(A)@Gu[i,j+1]
            termGv_mn = np.transpose(deltaZ)@inv(A)@Gv[i,j]
            d[i,j] = 1 + si[i+1,j] + si[i,j] + sj[i,j+1] + sj[i,j]
            q[i,j] = surf[i,j] - (dt/dx)*(termGu - termGu_mn) - (dt/dy)*(termGv - termGv_mn)
    for j in range(N):
        for i in range(N):
            #surf_vect.append(surf[i,j])
            q_vect.append(q[i,j])
            d_vect.append(d[i,j])
            si_vect.append(si[i,j])
            sj_vect.append(sj[i,j])
    print('Vectorizing and solving surface...')
    Mat = fivediag_maker(d_vect, si_vect[1:], sj_vect[1:], si_vect[0:-1], sj_vect[0:-1],  N)
    res = q_vect
    surf_vect = pp.solve(Mat, res)
    if np.isnan(surf_vect).any():
        print("WARNING: NaN in surface vector.")
    for j in range(N):
        surf[:,j] = surf_vect[j*N:(j+1)*N]
    surf[:,0] = surf[:,1]
    surf[:,-1] = surf[:,-2] 
    surf[0,:] = surf[1,:] 
    surf[-1,:] = surf[-2,:]
    return surf

def solve_U(surf,U,N,deltaZ,Gu,A,dt,dx):
    for i in range(N-1):
        for j in range(N):
            b = Gu[i,j] - (((9.81*dt)/dx)*(surf[i+1,j]-surf[i,j])*deltaZ)#.reshape(N)
            U[i,j] = lstsq(A, b, rcond=None)[0]
    return U

def solve_V(surf,V,N,deltaZ,Gv,A,dt,dy):
    for i in range(N):
        for j in range(N-1):
            b = Gv[i,j] - (((9.81*dt)/dy)*(surf[i,j+1]-surf[i,j])*deltaZ)#.reshape(N)
            V[i,j] = lstsq(A, b, rcond=None)[0]
    return V

def update_matrices(U,V,f,nu_h,Gu,Gv,deltaZ,N,dt,dx,dy,dz,nu_v):
    print('updating Gu, Gv, 3D arrays...')
    Fu, Fv = F_op(U, V, f, nu_h, dt, dx, dy, dz, N)
    # Setting up G arrays
    for i in range(N):  
        for j in range(N):
            Gu[i,j] = deltaZ*Fu[i,j]
            Gv[i,j] = deltaZ*Fv[i,j]
            str_x = nu_v*(U[i,j,0] - U[i,j,1])/dz
            str_y = nu_v*(V[i,j,0] - V[i,j,1])/dz
            Gu[i,j,0] = Gu[i,j,0] + dt*str_x
            Gv[i,j,0] = Gv[i,j,0] + dt*str_y 
    return Gu,Gv
    
#def solve_W(surf,deltaZ,Gu,Gv,A,dt,dz):
 #   return W
