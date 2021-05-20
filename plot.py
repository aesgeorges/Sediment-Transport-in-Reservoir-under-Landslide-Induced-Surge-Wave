import matplotlib.pyplot as plt
from matplotlib.pyplot import subplot
from matplotlib.animation import FuncAnimation

def plot(x,y):
    #fig, ax = plt.subplots()
    plt.plot(x,y)
    plt.show()

def plot_parts(c):
    partx = []
    partz = []
    fig, ax = plt.subplots()
    for particle in c.particles:
        partx.append(particle.x)
        partz.append(particle.z)
    ax.scatter(partx, partz)
    plt.show()
