from matplotlib.pyplot import xlabel, ylabel
from column import Column
from plot import *
from parameters import *
from advance import *

c = Column(N, H, SMALL, pop)
c.setup(N, H, A, B, C, E, Sq, kappa, SMALL, nu, g, z0, zb, u_crit, Ar, rho0, alpha, pop)
#plot(c.u, c.z)
#plot_parts(c)
t = [0]

fig, ax = plt.subplots()
ax = plt.axes(xlim =(-0.5, 0.5), ylim =(-H, 0)) 
ax.set(xlabel='Velocity U', ylabel='Depth - z')
ax.set_title('Animation of velocity with depth')
ax.grid()
line, = ax.plot([], [], lw = 3) 

def init(): 
    line.set_data([], [])
    return line,

def animate(i):
    #for m in range(1,M):
    t.append(dt*(i-1))
    advance(c, N, C_D, kappa, beta, dt, px0, t_px, t[i])
    print('advance... t = '+str(t[i])+' s.')
        #print(c.u)
    line.set_data(c.u, c.z)
    return line,

ani = FuncAnimation(fig, animate, frames = M, interval=dt)
plt.style.use('ggplot')

plt.show()

