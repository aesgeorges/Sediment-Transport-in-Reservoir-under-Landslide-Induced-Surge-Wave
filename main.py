from column import Column
from plot import *
from parameters import *
from advance import *
import matplotlib as mpl
import os

c = Column(N, H, SMALL)
c.setup(N, H, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
#plot(c.u, c.z)
#plot_parts(c)
t = []

fig, ax = plt.subplots()
ax = plt.axes(ylim =(-H, 0), xlim =(-50000, 50000)) 
norm = mpl.colors.Normalize(vmin=1.25, vmax=3.75)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='seismic_r'),
             ax=ax,  label='Specific Density of Particle')
ax.set(xlabel='Position - x', ylabel='Depth - z')
ax.set_title('Particle Tracking')
ax.autoscale
ax.set_facecolor("#000000") 

line, = ax.plot([], [], lw = 3) 


def init(): 
    line.set_data([], [])
    return line,

def animate(i):
    #for m in range(1,M):
    t.append(dt*(i))
    advance(c, N, C_D, kappa, beta, dt, px0, t_px, t[i], SMALL, zb, Sq, nu)
    print('advance... t = '+str(t[i])+' s.')
    #print(c.kz)
    line.set_data(c.u, c.z)
    return line,

def animate_parts(i):
    ax.clear()
    t.append(dt*(i))
    advance(c, N, C_D, kappa, beta, dt, px0, t_px, t[i], SMALL, zb, Sq, nu)
    #print('advance... t = '+str(t[i])+' s.')
    x = []
    z = []
    col = []
    for particle in c.particles:
        x.append(particle.x)
        z.append(particle.z)
        col.append(particle.rhor)
    ax.scatter(x, z, marker='+', c=col, cmap='seismic_r')
    ax.set(xlabel='Position - x', ylabel='Depth - z')
    ax.set_title('Particle Tracking')
    ax.set_xlim(-50000, 50000)
    ax.set_ylim(-H, 0)
    line.set_data(c.u, c.z)
    return line,
    

ani = FuncAnimation(fig, animate_parts, frames = M, interval=dt, repeat=False)

#ani.save('Visualization/particles_mix.gif', writer='imagemagick', fps=30)

plt.show()

