


import matplotlib

matplotlib.use('macosx')

import matplotlib.pyplot as plt

class DataPlotter(object):
    def __init__(self):
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, 1)
        plt.show(False)
        plt.draw()
        pass