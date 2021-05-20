import numpy as np
#from parameters import *
import math
# Time advancement functions
# Each function performs a single time-step for a corresponding component
# adv function calls all function to perform a time-step for all components
# Components are:
# particles, u, v, temp, rho, q2, q21, l, kz, nu_t, kq

# Function to advance all components of column
def advance(c, N, C_D, kappa, beta, dt, px0, t_px, t):
    Px = pressure(N, t_px, px0, t)
    c.u = velocity(c, N, C_D, kappa, beta, dt, Px)

# Pressure update
def pressure(N, t_px, px0, t):
    Px = np.ones(N)
    if t_px == 0:
        Px = Px * px0 # Steady and constant forcing
    else:
        Px = Px * px0*math.cos((2*math.pi*t)/(3600*t_px))
    return Px

# Velocity (U,V)
def velocity(c, N, C_D, kappa, beta, dt, Px): 
    up = c.u
    vp = c.v
    nu_tp = c.nu_t
    aU = np.zeros(N)
    bU = np.zeros(N)
    cU = np.zeros(N)
    dU = np.zeros(N)
    u = np.zeros(N)
    for i in range(1, N-1):
        aU[i] = -0.5*beta*(nu_tp[i] - nu_tp[i-1])
        bU[i] = 1+0.5*beta*(nu_tp[i+1] + 2*nu_tp[i] + nu_tp[i-1])
        cU[i] = -0.5*beta*(nu_tp[i] + nu_tp[i+1])
        dU[i] = up[i] - dt*Px[i]
    # Bottom-Boundary Condition: log-law
    bU[0] = 1+0.5*beta*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
    cU[0] = -0.5*beta*(nu_tp[1] + nu_tp[0])
    dU[0] = up[0] - dt*Px[0]
    # Top Boundary Condition: no-stress
    aU[N-1] = -0.5*beta*(nu_tp[N-1] + nu_tp[N-2])
    bU[N-1] = 1 + 0.5*beta*(nu_tp[N-1] + nu_tp[N-2])
    cU[N-1] = up[N-1] - dt*Px[N-1]
    # TDMA to solve for u
    for i in range(1, N-1):
        bU[i] = bU[i] - aU[i]/bU[i-1]*cU[i-1]
        dU[i] = dU[i] - aU[i]/bU[i-1]*dU[i-1]
    u[N-1] = dU[N-1]/bU[N-1]
    for i in range(N-2, 0, -1):
        u[i] = 1/bU[i]*(dU[i] - cU[i]*c.u[i+1])
    return u

# q2 - Turbulence paramater

def turb_q2(c, N, C_D, kappa, beta, dt, Px):
    up = c.q2
    aQ2 = np.zeros(N)
    bQ2 = np.zeros(N)
    cQ2 = np.zeros(N)
    dQ2 = np.zeros(N)
    q2 = np.zeros(N)
    for i in range(1, N-1):
        aU[i] = -0.5*beta*(nu_tp[i] - nu_tp[i-1])
        bU[i] = 1+0.5*beta*(nu_tp[i+1] + 2*nu_tp[i] + nu_tp[i-1])
        cU[i] = -0.5*beta*(nu_tp[i] + nu_tp[i+1])
        dU[i] = up[i] - dt*Px[i]