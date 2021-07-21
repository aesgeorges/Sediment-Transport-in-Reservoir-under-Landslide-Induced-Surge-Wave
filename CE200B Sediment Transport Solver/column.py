# Contains parameters of water column. Included are:
# U, V, Q2, Q2L, Q2L, L, temp, Gh, Sm, Sh, Kq, nu_t, Kz, particles
import numpy as np
import math
from particle import Particle

class Column:
    def __init__(self, N, H, SMALL):
        self.z = np.zeros(N) ##### In the future separate N in Nx and Nz
        self.x = np.zeros(N) ##### In the future separate N in Nx and Nz
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
        self.rho = np.zeros(N)
        #self.Laplace = np.zeros((N*N, N*N))
        self.particles = []

    def setup(self, N, H, Length, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop):
        dz = H/N
        dx = Length/N
        # Initializing Grid
        for i in range(N): ##### In the future separate N in Nx and Nz
            self.z[i] = -H + dz*((i+1) - 0.5)
        for i in range(N): ##### In the future separate N in Nx and Nz
            self.x[i] = Length - dx*((i+1) - 0.5)
        # Temperature Distribution
        de1C = 5 # Change in temperature at initial thermocline
        zde1C = -5 # Position of initial thermocline
        dzde1C = -2 # Width of initial thermocline 
        for i in range(N):
            self.temp[i] = 15
            if self.z[i] <= zde1C - 0.5*dzde1C:
                self.temp[i] = 15
            elif self.z[i] >= zde1C + 0.5*dzde1C:
                self.temp[i] = 15 + de1C
            else:
                self.temp[i] = 15 + de1C*(self.z[i] - zde1C + 0.5*dzde1C) / dzde1C
            self.rho[i] = rho0*(1-alpha*(self.temp[i] - 15))
        
        # Density distribution from Temperature
        self.rho = rho0*(1-alpha*(self.temp - 15))

        # Brunt-Vaisala Frequency from density profile
        self.n_bv[0] = math.sqrt((-g/rho0)*(self.rho[1]-self.rho[0])/(dz))
        for i in range(1,N-1):
            self.n_bv[i] = math.sqrt((-g/rho0)*(self.rho[i+1]-self.rho[i])/(dz))
        self.n_bv[N-1] = math.sqrt((-g/rho0)*(self.rho[N-1]-self.rho[N-2])/(dz))

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
        for i in range(pop):
            rand_rhor = np.random.normal(loc=2.5, scale=1.25)
            rand_D = np.random.normal(loc=0.0001, scale=0.00004975)
            rand_z = np.random.rand() * (-H/2)
            p = Particle(rand_z, rand_rhor, rand_D)
            self.particles.append(p)

    def laplacian(self, Nx, Nz, dxi, dzi):
        for j in range(Nz): 
            for i in range(Nx): 
                self.Laplace[i+(j-1)*Nx, i+(j-1)*Nx] = 2*(dxi**2) + 2*(dzi**2)
                for ii in range(i-1, i+1, 2):
                    if ii > 0 and ii <= Nx: # Interior Point
                        self.Laplace[i+(j-1)*Nx, ii + (j-1)*Nx] = -(dxi**2)
                    else:
                        self.Laplace[i+(j-1)*Nx, i+(j-1)*Nx] = self.Laplace[i+(j-1)*Nx, i+(j-1)*Nx] - (dxi**2)
                for jj in range(j-1, j+1, 2):
                    if jj > 0 and ii <= Nz: # Interior Point
                        self.Laplace[i+(j-1)*Nx, i+(jj-1)*Nx] = -(dzi**2)
                    else:
                        self.Laplace[i+(j-1)*Nx, i+(j-1)*Nx] = self.Laplace[i+(j-1)*Nx, i+(j-1)*Nx] - (dzi**2)
        self.Laplace[0,:] = 0
        self.Laplace[0, 0] = 1





        

        

    


        