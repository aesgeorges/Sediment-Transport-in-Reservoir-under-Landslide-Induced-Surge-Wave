# Physical Parameters
z0 = 0.01
zb = 10*z0
g = 9.81
C_D = 0.0025
small = 1e-6
kappa = 0.4
nu = 1e-6
u_crit = 0.007
SMALL = 1e-6
Ar = 1.25 # Erosion coefficient
rho0 = 1000 # kg/m3 - water density
alpha = 0
# Model Parameters
H = 5  # Column depth (m)
N = 80  # No. of grid points
dt = 60  # Size of time step (s)
M = 900  # No. of time steps
pop = 100  # No. of particles in column
dz = H/N
beta = (dt/dz)**2
# Turbulence Closure Parameters
A = [0.92, 0.74]
B = [16.6, 10.1]
C = [0.08]
E = [1.8, 1.33, 0.25]
Sq = 0.2