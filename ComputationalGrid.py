import numpy as np

class ComputationalGrid(object):
    def __init__(self):
        self.nx = 500
        self.x = np.arange(0, self.nx) * 1.0 / self.nx
