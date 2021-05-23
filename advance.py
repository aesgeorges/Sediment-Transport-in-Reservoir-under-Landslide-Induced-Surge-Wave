import numpy as np
from numpy.core.fromnumeric import shape
from parameters import N, A, B, C, E, H
import math
# Time advancement functions
# Each function performs a single time-step for a corresponding component
# adv function calls all function to perform a time-step for all components
# Components are:
# particles, u, v, temp, rho, q2, q21, l, kz, nu_t, kq

aU = np.zeros(N)
bU = np.zeros(N)
cU = np.zeros(N)
dU = np.zeros(N)
u = np.zeros(N)

aQ2 = np.zeros(N)
bQ2 = np.zeros(N)
cQ2 = np.zeros(N)
dQ2 = np.zeros(N)

aQ2l = np.zeros(N)
bQ2l = np.zeros(N)
cQ2l = np.zeros(N)
dQ2l = np.zeros(N)

q = np.zeros(N)
q2 = np.zeros(N)
q2l = np.zeros(N)
L = np.zeros(N)
kz = np.zeros(N)
kq = np.zeros(N)
nu_t = np.zeros(N)
sh = np.zeros(N)
sm = np.zeros(N)

q2p = []
q2lp = []
lp = []
kqp = []
kzp = []
nu_tp = []
up = []
vp = []
n_bvp = []
#u_star = 0

# Function to advance all components of column
def advance(c, N, C_D, kappa, beta, dt, px0, t_px, t, SMALL, zb, Sq, nu):
    [up, vp, nu_tp, q2p, q2lp, lp, kqp, kzp, n_bvp] = previous(c)
    u_star = abs(c.u[0]*math.sqrt(C_D))
    Px = pressure(N, t_px, px0, t)
    [c.sm, c.sh, c.gh] = update_param(c, SMALL, A, B, C)
    c.u = velocity(c, N, C_D, kappa, beta, dt, Px, up, nu_tp)  
    c.q2 = turb_q2(c, N, beta, dt, B, u_star, up, nu_tp, kqp, kzp, lp, n_bvp, q2p, SMALL)
    c.q2L = turb_q2l(c, N, beta, dt, B, E, u_star, kappa, H, up, nu_tp, zb, kqp, kzp, q2lp, lp, n_bvp, q2p, SMALL)
    [c.kq, c.nu_t, c.kz, c.q, c.q2L, c.L] = turbmix(c, SMALL, zb, Sq, nu)

# Save previous step
def previous(c):
    up = c.u
    vp = c.v
    nu_tp = c.nu_t
    q2p = c.q2
    q2lp = c.q2L
    lp = c.L
    kqp = c.kq
    kzp = c.kz
    n_bvp = c.n_bv
    return up, vp, nu_tp, q2p, q2lp, lp, kqp, kzp, n_bvp

# Pressure update
def pressure(N, t_px, px0, t):
    Px = np.ones(N)
    if t_px == 0:
        Px = Px * px0 # Steady and constant forcing
    else:
        Px = Px * px0*math.cos((2*math.pi*t)/(3600*t_px))
    return Px

# Update model paramaeters, Sm and Sh
def update_param(c, SMALL, A, B, C):
    for i in range(N):
        gh = -(c.n_bv[i]*c.L[i]/(c.q[i]+SMALL))
        gh = min(gh, 0.0233)
        gh = max(gh, -0.28)
        num = B[0]**(-1/3) - A[0]*A[1]*gh*((B[1]-3*A[1])*(1-6*A[0]/B[0])-3*C[0]*(B[1]+6*A[0]))
        dem = (1-3*A[1]*gh*(B[1]+6*A[0]))*(1-9*A[0]*A[1]*gh)
        sm[i] = num/dem
        sh[i] = A[1]*(1-6*A[0]/B[0])/(1-3*A[1]*gh*(B[1]+6*A[0]))
    return sm, sh, gh

def TDMA(a, b, c, d):
    n = len(d)
    w= np.zeros(n-1,float)
    g= np.zeros(n, float)
    p = np.zeros(n,float)
    
    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]

    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p

# Velocity (U,V)
def velocity(c, N, C_D, kappa, beta, dt, Px, up, nu_tp): 
    for i in range(1, N-1):
        aU[i] = -0.5*beta*(nu_tp[i] + nu_tp[i-1])
        bU[i] = 1+0.5*beta*(nu_tp[i+1] + 2*nu_tp[i] + nu_tp[i-1])
        cU[i] = -0.5*beta*(nu_tp[i] + nu_tp[i+1])
        dU[i] = up[i] - dt*Px[i]
    # Bottom-Boundary Condition: log-law
    bU[0] = 1+0.5*beta*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
    cU[0] = -0.5*beta*(nu_tp[1] + nu_tp[0])
    dU[0] = up[0] - dt*Px[0]
    # Top Boundary Condition: no-stress
    aU[-1] = -0.5*beta*(nu_tp[-1] + nu_tp[N-2])
    bU[-1] = 1 + 0.5*beta*(nu_tp[-1] + nu_tp[N-2])
    dU[-1] = up[-1] - dt*Px[-1]
    # TDMA to solve for u
    #u = TDMA(aU, bU, cU, dU)
    for i in range(1, N):
        bU[i] = bU[i] - aU[i]/bU[i-1]*cU[i-1]
        dU[i] = dU[i] - aU[i]/bU[i-1]*dU[i-1]
    u[-1] = dU[-1]/bU[-1]
    for i in range(N-2, -1, -1):
        u[i] = (1/bU[i])*(dU[i] - cU[i]*u[i+1])
    return u

# q2 - Turbulence Paramater
def turb_q2(c, N, beta, dt, B, u_star, up, nu_tp, kqp, kzp, lp, n_bvp, q2p, SMALL):
    for i in range(1, N-1):
        diss = (2 * dt *(q2p[i]**0.5))/(B[0]*lp[i]) # Coefficient for linearized term
        aQ2[i] = -0.5*beta*(kqp[i] + kqp[i-1])
        bQ2[i] = 1+0.5*beta*(kqp[i+1] + 2*kqp[i] + kqp[i-1]) + diss
        cQ2[i] = -0.5*beta*(kqp[i] + kqp[i+1])
        dQ2[i] = q2p[i] + 0.25*beta*nu_tp[i]*(up[i+1]-up[i-1])**2 -dt*kzp[i]*(n_bvp[i]**2)
    # Bottom-Boundary Condition 
    q2bot = B[0]**(2/3) * u_star**2
    bdryterm = 0.5*beta*kqp[0]*q2bot
    diss =  2 * dt *((q2p[0]**0.5)/(B[0]*lp[0]))
    bQ2[0] = 1+0.5*beta*(kqp[1] + kqp[0]) + diss
    cQ2[0] = -0.5*beta*(kqp[1] + kqp[0])
    dQ2[0] = q2p[0] + dt*((u_star**4)/nu_tp[0]) - dt*kzp[0]*(n_bvp[0]**2) + bdryterm
    # Top-Boundary Condition
    diss =  2 * dt *((q2p[-1]**0.5)/(B[0]*lp[-1]))
    aQ2[-1] = -0.5*beta*(kqp[-1] + kqp[N-2])
    bQ2[-1] = 1+0.5*beta*(kqp[-1] + 2*kqp[-1] + c.kq[N-2]) + diss
    dQ2[-1] = q2p[-1] + 0.25*beta*nu_tp[-1]*((up[-1] - up[N-2])**2) -4*dt*kzp[-1]*(n_bvp[-1]**2)
    # TDMA to solve for q2
    #q2 = TDMA(aQ2, bQ2, cQ2, dQ2)
    for i in range(1, N):
        bQ2[i] = bQ2[i] - aQ2[i]/bQ2[i-1]*cQ2[i-1]
        dQ2[i] = dQ2[i] - aQ2[i]/bQ2[i-1]*dQ2[i-1]
    q2[-1] = dQ2[-1]/bQ2[-1]
    for i in range(N-2, -1, -1):
        q2[i] = (1/bQ2[i])*(dQ2[i] - cQ2[i]*q2[i+1])
    for i in range(N):
        if q2[i] < SMALL:
            q2[i] = SMALL
    return q2

def turb_q2l(c, N, beta, dt, B, E, u_star, kappa, H, up, nu_tp, zb, kqp, kzp, q2lp, lp, n_bvp, q2p, SMALL):
    for i in range(1, N-1):
        diss = 2*dt*((q2p[i]**0.5) / (B[0]*lp[i]))*(1+E[1]*(lp[i]/(kappa*abs(-H-c.z[i])))**2 + E[2]*(lp[i]/(kappa*abs(c.z[i])))**2)
        aQ2l[i] = -0.5*beta*(kqp[i] + kqp[i-1])
        bQ2l[i] = 1+0.5*beta*(kqp[i+1] + 2*kqp[i] + kqp[i-1]) + diss
        cQ2l[i] = -0.5*beta*(kqp[i] + kqp[i+1])
        dQ2l[i] = q2lp[i] + 0.25*beta*nu_tp[i]*E[0]*lp[i]*(up[i+1]-up[i-1])**2 - 2*dt*lp[i]*E[0]*kzp[i]*(n_bvp[i]**2)
    # Bottom-Boundary Condition
    q2lbot = B[0]**(2/3) * u_star**2 * kappa * zb
    bdryterm = 0.5*beta*kqp[0]*q2lbot
    diss =  2 * dt *(q2p[0]**0.5)/(B[0]*lp[0])*(1+E[1]*(lp[0]/(kappa*abs(-H-c.z[0])))**2 + E[2]*(lp[0]/(kappa*abs(c.z[0])))**2)
    bQ2l[0] = 1+0.5*beta*(kqp[1] + kqp[0]) + diss
    cQ2l[0] = -0.5*beta*(kqp[1] + kqp[0])
    dQ2l[0] = q2lp[0] + dt*((u_star**4)/nu_tp[0])*E[0]*lp[0] - dt*lp[0]*E[0]*kzp[0]*n_bvp[0]**2 + bdryterm
    # Top-Boundary Condition
    diss =  2 * dt *(q2p[-1]**0.5)/(B[0]*lp[-1])*(1+E[1]*(lp[-1]/(kappa*abs(-H-c.z[-1])))**2 + E[2]*(lp[-1]/(kappa*abs(c.z[-1])))**2)
    aQ2l[-1] = -0.5*beta*(kqp[-1] + kqp[N-2])
    bQ2l[-1] = 1+0.5*beta*(kqp[-1] + 2*kqp[-1] + kqp[N-2]) + diss # Are we using kq or kqp here?
    dQ2l[-1] = q2lp[-1] + 0.25*beta*nu_tp[-1]*E[0]*lp[-1]*(up[-1]-up[N-2])**2 - 2*dt*lp[-1]*E[0]*kzp[-1]*(n_bvp[-1]**2)
    # TDMA to solve for q2
    #q2l = TDMA(aQ2l, bQ2l, cQ2l, dQ2l)
    for i in range(1, N):
        bQ2l[i] = bQ2l[i] - aQ2l[i]/bQ2l[i-1]*cQ2l[i-1]
        dQ2l[i] = dQ2l[i] - aQ2l[i]/bQ2l[i-1]*dQ2l[i-1]
    q2l[-1] = dQ2l[-1]/bQ2l[-1]
    for i in range(N-2, -1, -1):
        q2l[i] = (1/bQ2l[i])*(dQ2l[i] - cQ2l[i]*q2l[i+1])
    for i in range(N):
        if q2l[i] < SMALL:
            q2l[i] = SMALL
    return q2l


# Turbulent Lengthscale (L) and mixing coefficients (kz, nu_t, kq)
def turbmix(c, SMALL, zb, Sq, nu):
    for i in range(N):
        q[i] = math.sqrt(c.q2[i])
        L[i] = c.q2L[i]/(c.q2[i] + SMALL)
        # Limit due to stable stratification
        if ((L[i]**2)*c.n_bv[i])**2 > 0.281*c.q2[i]:
            # Adjust Q2L as well as L
            q2l[i] = c.q2[i]*math.sqrt(0.281*c.q2[i]/(c.n_bv[i]**2 + SMALL))
            L[i] = q2l[i] / c.q2[i]

        # Keep L from becoming zero -- zb=bottom roughness parameter
        if abs(L[i]) < zb:
            L[i] = zb
        # Update diffusivities 
        kq[i] = Sq*q[i]*L[i] + nu
        nu_t[i] = c.sm[i]*q[i]*L[i] + nu
        kz[i] = c.sh[i]*q[i]*L[i] + nu
    return kq, nu_t, kz, q, q2l, L
    