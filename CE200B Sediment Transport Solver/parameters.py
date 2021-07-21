# Physical Parameters
z0 = 0.01
zb = 10*z0
g = 9.81
C_D = 0.0025
kappa = 0.4
nu = 1e-6
u_crit = 0.007
SMALL = 1e-6
Ar = 1.25 # Erosion coefficient
rho0 = 1000 # kg/m3 - water density
alpha = 0 # Thermal expansivity
px0 = 0.0001 # Magnitude of pressure gradient forcing
t_px = 0 # Periodic pressure forcing

# Model Parameters
H = 10 # Column depth (m)
Length = 10 # Plane Length (m)
N = 80  # No. of grid points
dt = 60  # Size of time step (s)
M = 900  # No. of time steps
pop = 100  # No. of particles in column
dz = H/N
dx = Length/N
beta = dt/(dz**2)

# 3 Columns
H1 = 25
H2 = 17.5
H3 = 10
dz1 = H1/N
dz2 = H2/N
dz3 = H3/N
beta1 = dt/(dz1**2)
beta2 = dt/(dz2**2)
beta3 = dt/(dz3**2)

# Turbulence Closure Parameters
A = [0.92, 0.74]
B = [16.6, 10.1]
C = [0.08]
E = [1.8, 1.33, 0.25]
Sq = 0.2