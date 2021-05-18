import numpy as np

class Particle:
    def __init__(self, y, rhor):
        self.x = 0
        self.y = np.random.rand()*y
        self.rhor = rhor
