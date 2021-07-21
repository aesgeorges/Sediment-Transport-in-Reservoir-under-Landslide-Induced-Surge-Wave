from column import Column
from plot import *
from parameters import *
from advance import *
from numba import jit
import matplotlib as mpl
import os

#dxi = 1/dx
#dzi = 1/dz
#Nx = N ####### Temporary
#Nz = N ####### Temporary
#c = Column(N, H, SMALL)
#c.setup(N, H, Length, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
#c.laplacian(Nx, Nz, dxi, dzi)
#print(c.Laplace)
#plot(c.u, c.z)
#plot_parts(c)

###### 3 Columns
c1 = Column(N, H1, SMALL)
c2 = Column(N, H2, SMALL)
c3 = Column(N, H3, SMALL)
c1.setup(N, H1, Length, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
c2.setup(N, H2, Length, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
c3.setup(N, H3, Length, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
t = []

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
norm = mpl.colors.Normalize(vmin=1.25, vmax=3.75)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='seismic_r'),
             ax=ax,  label='Specific Density of Particle')
ax.autoscale
#ax.set_facecolor("#000000") 

line, = ax.plot([], [], lw = 3) 


def init(): 
    line.set_data([], [])
    return line,

def animate(i):
    #for m in range(1,M):
    t.append(dt*(i))
    advance(c1, N, C_D, kappa, beta1, dt, px0, t_px, t[i], SMALL, zb, Sq, nu)
    #print('advance... t = '+str(t[i])+' s.')
    #print(c.kz)
    ax.set_title('Velocity U Profile')
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-H, 0)
    line.set_data(c.u, c.z)
    return line,

@jit
def animate_parts(i):
    ax.clear()
    t.append(dt*(i))
    cols = [2,1,0]
    c = [c1, c2, c3]
    beta = [beta1, beta2, beta3]
    H = [H1, H2, H3]
    for ytick in cols:
        advance(c[ytick], N, H[ytick], C_D, kappa, beta[ytick], dt, px0, t_px, t[i], SMALL, zb, Sq, nu)
        #print('advance... t = '+str(t[i])+' s.')
        x = []
        z = []
        col = []
        for particle in c[ytick].particles:
            x.append(particle.x)
            z.append(particle.z)
            col.append(particle.rhor)
        ax.scatter(x, ytick, z, marker='+', c=col, cmap='seismic_r')
        ax.set(xlabel='Position - x', ylabel='Depth - z')
        ax.set_title('Particle Tracking - dt=' + str(dt) + 's, N=' + str(N) + ', ' + str(pop) + ' particles.')
        ax.set_xlim(-50000, 0)
        ax.set_ylim(0, 2)
        ax.set_zlim(-H1, 0)
        ax.w_xaxis.set_pane_color((0,0,0,1))
        ax.w_yaxis.set_pane_color((0,0,0,1))
        ax.w_zaxis.set_pane_color((0,0,0,1))
        #line.set_data(c.u, c.z)
        #return line,
    

ani = FuncAnimation(fig, animate_parts, frames = M, interval=dt, repeat=False)

#ani.save('Visualization/particles_3D.gif', writer='imagemagick', fps=24)

plt.show()

