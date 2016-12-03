#!/usr/bin/env python3
from __future__ import print_function
import numpy as np
import matplotlib

matplotlib.use('macosx')

import matplotlib.pyplot as plt
import time


class mysim(object):
    def __init__(self):
        self.unk = np.ndarray(shape=(3, 10), dtype=float)
        print("test")
        self.fig, self.ax = plt.subplots(1, 1)
        plt.show(False)
        plt.draw()
        self.main_loop()
        time.sleep(10)
        pass

    def set_initial_conditions(self):
        rho=10
        vx=10
        p=10
        self.unk[0,:]=rho
        self.unk[1,:]=vx
        self.unk[2,:]=p



    def getflux(self):
        rho = self.unk[0, :]
        v1 = self.unk[1, :]
        p = self.unk[2, :]
        self.mass_flux = rho * v1
        self.mom_flux = rho * v1 * v1 + p
        self.en_flux = rho * v1 ^ 3 + v1 * p

    def time_advance(self):
        self.unk_n = self.unk

    def main_loop(self):
        self.set_initial_conditions()
        for step in range(0,10):
            self.time_advance()
            self.plot_all()
        pass

    def plot_all(self):
        rho = self.unk[0, :]
        print("plotting")
        pts = self.ax.plot(rho)
        plt.pause(0.05)
        self.fig.canvas.draw()


if __name__ == '__main__':
    a = mysim()
