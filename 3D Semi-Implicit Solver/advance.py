import numpy as np
from numpy.linalg import inv, lstsq
from scipy.sparse import diags
import pentapy as pp


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
    #ai = np.zeros((N,N))
    #aj = np.zeros((N,N))
    q = np.zeros((N,N))
    #b = np.zeros((N,N))
    #e = np.zeros((N,N))
    #Mp = np.zeros((N,N))
    #r = np.zeros((N,N))
    print('surface calcs...')
    print('computing compact terms...')
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
            #ai[i,j] = si[i,j]/(math.sqrt(d[i,j]*d[i,j]))
            #ai[i+1,j] = si[i+1,j]/(math.sqrt(d[i,j]*d[i+1,j]))
            #aj[i,j] = sj[i,j]/(math.sqrt(d[i,j]*d[i,j]))
            #aj[i,j+1] = sj[i,j+1]/(math.sqrt(d[i,j]*d[i,j+1]))
            #b[i,j] = q[i,j]/(math.sqrt(d[i,j]))
    q_vect = []
    d_vect = []
    surf_vect = []
    si_vect = []
    sj_vect = []
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
    for j in range(N):
        surf[:,j] = surf_vect[j*N:(j+1)*N]
    
    #p = np.zeros((N,N))
    #print('setting up gradient method...')
    #for i in range(N-1):
    #    for j in range(N-1):
    #        e[i,j] = surf[i,j]*math.sqrt(d[i,j])
    #        r[i,j] = e[i,j] - ai[i+1,j]*e[i+1,j] - ai[i,j]*e[i-1,j] - aj[i,j+1]*e[i,j+1] - aj[i,j]*e[i,j-1] - b[i,j]
    #        p[i,j] = r[i,j]
    #        Mp[i,j] = p[i,j] - ai[i+1,j]*p[i+1,j] - ai[i,j]*p[i-1,j] - aj[i,j+1]*p[i,j+1] - aj[i,j]*p[i,j-1]
    #print('running conjugate method algorithm...')
    #enew = np.zeros((N,N))
    #var = np.dot(r,r)/np.dot(p,Mp)
    #while np.dot(r.all(),r.all()) > epsilon:
    #    print("in?")
    #    alpha = np.dot(r,r)/np.dot(p,Mp)
    #    enew = e - alpha*p
    #    rnew = r - alpha*Mp
    #    beta = np.dot(rnew,rnew)/np.dot(r,r)
    #    p = r + beta*p
    #surf = enew*np.sqrt(d)
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


#def solve_W(surf,deltaZ,Gu,Gv,A,dt,dz):
 #   return W
