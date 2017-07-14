from sys import platform

import matplotlib

matplotlib.use('TkAgg')
if platform == "linux" or platform == "linux2":
    matplotlib.use('TkAgg')
    # linux
elif platform == "darwin":
    matplotlib.use('macosx')
    # OS X
elif platform == "win32":
    matplotlib.use('TkAgg')

import matplotlib.pyplot as plt


class DataPlotter(object):
    def __init__(self):
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, sharex=True)
        plt.show(False)
        plt.draw()

    def plot_all(self, unk_n, xaxis, step):
        rho = unk_n[0, :]
        vel = unk_n[1, :] / rho
        p = (unk_n[2, :] - .5 * rho * vel ** 2)
        self.ax1.cla()
        self.ax1.plot(xaxis, rho, '-b')
        self.ax1.set_ylabel('rho')
        self.ax2.cla()
        self.ax2.plot(xaxis, vel, '-r')
        self.ax2.set_ylabel('v')
        self.ax3.cla()
        self.ax3.plot(xaxis, p, '-g')
        self.ax3.set_ylabel('p')
        plt.pause(0.05)
        self.fig.canvas.draw()
        plt.savefig('output' + str(step).zfill(5) + '.png')
