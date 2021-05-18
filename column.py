# Contains parameters of water column. Included are:
# U, V, Q2, Q2L, Q2L, L, temp, Gh, Sm, Sh, Kq, nu_t, Kz, particles
import numpy as np
import math

from particle import Particle


class Column:
    def __init__(self, N, H, SMALL, pop):
        self.z = np.linspace(-H, 0, N)
        self.u = np.zeros(N)
        self.v = np.zeros(N)
        self.q2 = np.ones(N)*SMALL # Turbulent field is seeded with small values then let evolve
        self.q2L = np.ones(N)*SMALL 
        self.q = np.zeros(N)
        self.L = np.zeros(N)
        self.temp = np.zeros(N)
        self.gh = np.zeros(N)
        self.sm = np.zeros(N)
        self.sh = np.zeros(N)
        self.nu_t = np.zeros(N) # Turbulent viscosity
        self.kq = np.zeros(N) # Turbulent diffusivity for q2
        self.kz = np.zeros(N) # Turbulent scalar diffusivity
        self.particles = []

    def setup(self, N, H, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop):
        dz = H/N

        # Temperature Distribution
        for t in self.temp:
            t = 15
        
        # Density distribution from Temperature
        rho = rho0*(1-alpha*(self.temp - 15))

        # Brunt-Vaisala Frequency from density profile
        #brunt = math.sqrt((-g/rho0)*(rho(2)-rho(1))/(dz))
        n_bv = np.zeros(N)
        for i in range(N-2):
            n_bv[i] = math.sqrt((-g/rho0)*(rho[i+1]-rho[i])/(dz))
        n_bv[N-1] = math.sqrt((-g/rho0)*(rho[N-1]-rho[N-2])/(dz))

        # Initial Conditions
        for i in range(N-1):
            self.q[i] = math.sqrt(self.q2[i])

        self.L = -kappa*H*(self.z/H)*(1-(self.z/H))
        
        self.gh = ((n_bv*self.L)/(self.q + SMALL))**2
        for g in self.gh:
            g = min(g, 0.0233)
            g = max(g, -0.28)
        
        num= B[0]**(-1/3) - A[0]*A[1]*self.gh*((B[1]-3*A[1])*(1-6*A[0]/B[0])-3*C[0]*(B[1]+6*A[0]))
        dem= (1-3*A[1]*self.gh*(B[1]+6*A[0]))*(1-9*A[0]*A[1]*self.gh)
        self.sm = num/dem
        self.sh = A[1]*(1-6*A[0]/B[0])/(1-3*A[1]*self.gh*(B[1]+6*A[0]))
        self.nu_t = self.sm * self.q * self.L + nu
        self.kq = Sq*self.q*self.L + nu
        self.kz = self.sh * self.q * self.L + nu

        # Particles generation
        # Will implement gaussian distribution later
        for i in range(N-1):
            rand_rhor = np.random.normal(loc=1.25, scale=0.75)
            rand_D = np.random.normal(loc=0.005025, scale=0.004975)
            rand_z = np.random.rand() * (-H/2)
            p = Particle(rand_z, rand_rhor, rand_D)
            self.particles.append(p)
            print(p)



        

        

    


        