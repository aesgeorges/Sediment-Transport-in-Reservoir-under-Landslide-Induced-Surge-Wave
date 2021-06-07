import math, numpy as np
class Particle:
    def __init__(self, z, rhor, D):
        self.x = 0
        self.z = z
        self.rhor = rhor
        self.D = D

    def tracking(self, g, nu, u, dz, dt, N, u_star, u_crit, H, kzp):
        ws = ((self.rhor - 1)*g* self.D**2)/(18*nu) # Settling velocity
        index = min(round(0.5+abs(self.z/dz)), N) - 1 # Get position index
        
        # Compute new position
        x = self.x + u[index]*dt
        z = self.z - ws*dt + math.sqrt(2*abs(kzp[index])*dt)*np.random.randn()

        # Keeping particle within domain
        if z > 0:
            z = 0
        # Resuspension
        psus = 0.04*(u_star**2 - u_crit)
        if z < -H:
            z = -H
            if u_star**2 < u_crit**2:
                psus = 0
            elif u_star**2 > (u_crit**2 + (1/0.04)):
                psus = 1
            if np.random.rand() < psus:
                z = z + dz
            else:
                z = -H
        return x, z


        
