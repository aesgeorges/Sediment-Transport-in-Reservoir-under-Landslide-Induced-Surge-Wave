# Contains parameters of water column. Included are:
# U, V, Q2, Q2L, Q2L, L, temp, Gh, Sm, Sh, Kq, nu_t, Kz, particles
import numpy as np
import math

from particle import Particle


class Column:
    def __init__(self, N, H, SMALL, pop):
        self.z = np.zeros(N)
        self.u = np.zeros(N)
        self.v = np.zeros(N)
        self.q2 = np.ones(N)*SMALL # Turbulent field is seeded with small values then let evolve
        self.q2L = np.ones(N)*SMALL 
        self.q = np.zeros(N)
        self.L = np.zeros(N)
        self.temp = np.zeros(N)
        self.gh = 0
        self.sm = np.zeros(N)
        self.sh = np.zeros(N)
        self.nu_t = np.zeros(N) # Turbulent viscosity
        self.kq = np.zeros(N) # Turbulent diffusivity for q2
        self.kz = np.zeros(N) # Turbulent scalar diffusivity
        self.n_bv = np.zeros(N)
        self.particles = []

    def setup(self, N, H, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop):
        dz = H/N

        # Initializing Grid
        for i in range(N):
            self.z[i] = -H + dz*((i+1) - 0.5)
        # Temperature Distribution
        for t in self.temp:
            t = 15
        
        # Density distribution from Temperature
        rho = rho0*(1-alpha*(self.temp - 15))

        # Brunt-Vaisala Frequency from density profile
        #brunt = math.sqrt((-g/rho0)*(rho(2)-rho(1))/(dz))
        #n_bv = np.zeros(N)
        self.n_bv[0] = math.sqrt((-g/rho0)*(rho[1]-rho[0])/(dz))
        for i in range(1,N-1):
            self.n_bv[i] = math.sqrt((-g/rho0)*(rho[i+1]-rho[i])/(dz))
        self.n_bv[N-1] = math.sqrt((-g/rho0)*(rho[N-1]-rho[N-2])/(dz))

        # Initial Conditions
        for i in range(N):
            self.q[i] = math.sqrt(self.q2[i])

            self.L[i] = -kappa*H*(self.z[i]/H)*(1-(self.z[i]/H))
        
            self.gh = ((self.n_bv[i]*self.L[i])/(self.q[i] + SMALL))**2
            self.gh = min(self.gh, 0.0233)
            self.gh = max(self.gh, -0.28)
        
            num= B[0]**(-1/3) - A[0]*A[1]*self.gh*((B[1]-3*A[1])*(1-6*A[0]/B[0])-3*C[0]*(B[1]+6*A[0]))
            dem= (1-3*A[1]*self.gh*(B[1]+6*A[0]))*(1-9*A[0]*A[1]*self.gh)
            self.sm[i] = num/dem
            self.sh[i] = A[1]*(1-6*A[0]/B[0])/(1-3*A[1]*self.gh*(B[1]+6*A[0]))
            self.nu_t[i] = self.sm[i] * self.q[i] * self.L[i] + nu
            self.kq[i] = Sq*self.q[i]*self.L[i] + nu
            self.kz[i] = self.sh[i] * self.q[i] * self.L[i] + nu

        # Particles generation
        # Will implement gaussian distribution later
        for i in range(N):
            rand_rhor = np.random.normal(loc=1.25, scale=0.75)
            rand_D = np.random.normal(loc=0.005025, scale=0.004975)
            rand_z = np.random.rand() * (-H/2)
            p = Particle(rand_z, rand_rhor, rand_D)
            self.particles.append(p)



        

        

    


        