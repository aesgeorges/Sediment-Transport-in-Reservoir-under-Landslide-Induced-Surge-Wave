from column import Column
from plot import plot
from parameters import *

c = Column(N, H, small, pop)
c.setup(N, H, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
plot(c.u, c.z)