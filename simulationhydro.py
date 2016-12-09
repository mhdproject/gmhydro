#!/usr/bin/env python3
from __future__ import print_function
import numpy as np
import matplotlib

matplotlib.use('macosx')

import matplotlib.pyplot as plt
import time


class Mesh(object):
    def __init__(self):
        self.nx = 50


class SimulationHydro(object):
    def __init__(self):
        self.nx = 50
        self.steps=2
        self.dtdx = 1.0
        self.nvar = 3
        self.cfl = 0.5
        self.gamma=5./3.
        self.unk = np.ndarray(shape=(self.nvar, self.nx), dtype=float)
        self.unk_n = np.empty_like(self.unk)
        self.flux = np.empty_like(self.unk)
        self.sound_speed = np.empty_like(self.unk)
        print("test")
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, 1)
        plt.show(False)
        plt.draw()
        self.main_loop()
        time.sleep(10)
        pass

    def set_initial_conditions(self):
        nx = self.nx
        rho = 1.
        vx = 0
        p = 10.
        self.unk[0, 0:nx / 2] = rho
        self.unk[1, 0:nx / 2] = vx
        self.unk[2, 0:nx / 2] = p
        self.unk[0, nx / 2:] = rho
        self.unk[1, nx / 2:] = vx
        self.unk[2, nx / 2:] = p / 10.

    def get_sound_speed(self):
        rho = self.unk[0, :]
        p = self.unk[2, :]
        self.sound_speed = np.sqrt(self.gamma * p / rho)

    def get_max_speed(self):
        self.get_sound_speed()
        max_speed= self.unk[1,:]+self.sound_speed
        self.dtdx = 1. / max_speed.max()
        pass

    def get_flux(self):
        rho = self.unk[0, :]
        v1 = self.unk[1, :]/rho
        ke = 0.5*rho*v1*v1
        e_int = self.unk[2, :]-ke
        p =e_int *(self.gamma-1)
        self.flux[0, :] = rho * v1
        self.flux[1, :] = rho * v1 * v1 + p
        self.flux[2, :] = 0.5 * rho * v1 * v1 * v1 + v1 * p

    def boundary_conditions(self):
        self.unk_n[:,0]=self.unk[:,1]
        self.unk_n[:,-1]=self.unk[:,-2]

    def time_advance(self):
        self.get_flux()
        for i in range(1, self.nx - 1):
            self.unk_n[:, i] = 0.5 * (self.unk[:, i - 1] + self.unk[:, i + 1]) - self.cfl * self.dtdx * self.flux[:,
                                                                                                        i + 1]

        self.boundary_conditions()
        self.unk = self.unk_n

    def main_loop(self):
        self.set_initial_conditions()
        for step in range(0, self.steps):
            print("Timestep: ", step)
            self.time_advance()
            self.plot_all()
        pass

    def plot_all(self):
        rho = self.unk[0, :]
        vel = self.unk[1, :]
        p = self.unk[2, :]
        pts = self.ax1.plot(rho)
        pts = self.ax2.plot(vel)
        pts = self.ax3.plot(p)
        plt.pause(0.05)
        self.fig.canvas.draw()


if __name__ == '__main__':
    a = SimulationHydro()
