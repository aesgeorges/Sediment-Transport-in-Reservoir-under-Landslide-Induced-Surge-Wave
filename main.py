from column import Column
from plot import *
from parameters import *
from advance import *

c = Column(N, H, SMALL)
c.setup(N, H, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
#plot(c.u, c.z)
#plot_parts(c)
t = []

fig, ax = plt.subplots()
ax = plt.axes(ylim =(-H, 0), xlim =(-15000, 15000)) 
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

def animate_scatter(i):
    ax.clear()
    t.append(dt*(i))
    advance(c, N, C_D, kappa, beta, dt, px0, t_px, t[i], SMALL, zb, Sq, nu)
    #print('advance... t = '+str(t[i])+' s.')
    x = []
    z = []
    for particle in c.particles:
        x.append(particle.x)
        z.append(particle.z)
    ax.scatter(x, z, marker='+')
    ax.set(xlabel='Position - x', ylabel='Depth - z')
    ax.set_title('Particle Tracking')
    ax.set_xlim(-20000, 20000)
    ax.set_ylim(-H, 0)
    line.set_data(c.u, c.z)
    return line,
    

ani = FuncAnimation(fig, animate_scatter, frames = M, interval=dt, repeat=True)


plt.show()

