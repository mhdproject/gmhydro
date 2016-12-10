#!/usr/bin/env python3
from __future__ import print_function
import numpy as np
import matplotlib

matplotlib.use('macosx')

import matplotlib.pyplot as plt
import time


class ComputationalGrid(object):
    def __init__(self):
        self.nx = 50


class DataPlotter(object):
    def __init__(self):
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, 1)
        plt.show(False)
        plt.draw()
        pass


class SimulationHydro(object):
    def __init__(self):
        self.grid = ComputationalGrid()
        # self.grid.nx = 50
        self.dx = 1. / self.grid.nx
        self.steps = 10
        self.step = 0
        self.dtdx = 1.0
        self.nvar = 3
        self.cfl = 0.5
        self.gamma = 7. / 5.
        self.unk = np.ndarray(shape=(self.nvar, self.grid.nx), dtype=float)
        self.unk_n = np.empty_like(self.unk)
        self.flux = np.empty_like(self.unk)
        self.sound_speed = np.empty_like(self.unk)
        print("test")
        self.plotter = DataPlotter()

        self.main_loop()
        time.sleep(10)
        pass

    def set_initial_conditions(self):
        nx = self.grid.nx
        rho = 1.
        vx = 0
        p = 1.
        pr = p / 10.
        ke = 0.5 * rho * vx * vx
        self.unk[0, 0:nx / 2] = rho
        self.unk[1, 0:nx / 2] = rho * vx
        self.unk[2, 0:nx / 2] = ke + p / (self.gamma - 1)
        self.unk[0, nx / 2:] = 0.125
        self.unk[1, nx / 2:] = rho * vx
        self.unk[2, nx / 2:] = ke + pr / (self.gamma - 1)

    def get_sound_speed(self):
        rho = self.unk[0, :]
        vx = self.unk[1, :]
        ke = 0.5 * rho * vx * vx
        e_tot = self.unk[2, :]
        e_int = e_tot - ke
        p = e_int * (self.gamma - 1)
        self.sound_speed = np.sqrt(self.gamma * p / rho)

    def get_prim(self, unk):
        rho_r = unk[0]
        v_r = unk[1]
        etot_r = unk[2]
        esp_r = etot_r / rho_r
        p_r = (self.gamma - 1) * (etot_r - (rho_r * v_r ^ 2) / 2)
        h_r = (e_tot_r + p_r) / rho_r
        return v_r, p_r, h_r

    def riemann_solver(self):
        # the culbert b. laney way!
        # step 1: compute primitives
        rho_l = self.unk[0, i]
        rho_r = self.unk[0, i + 1]
        v_r, p_r, h_r = get_prim(self.unk[i])
        v_l, p_l, h_l = get_prim(self.unk[i + 1])

        # step :2 roe averages
        srl = np.sqrt(rho_l)
        srr = np.sqrt(rho_r)
        rho_roe = srl*srr
        v_roe = (srr * v_r + srl * v_l) / (srl + srr)
        h_roe = (srr * h_r + srl * h_l) / (srl + srr)
        a_roe = np.sqrt((self.gamma - 1) * (h_roe - 0.5 * v_roe * v_roe))
        a_roe2 = a_roe * a_roe
        # step:3 compute eigenvalues
        self.ev[0] = v_roe
        self.ev[1] = v_roe + a_roe
        self.ev[2] = v_roe - a_roe
        # step:4 compute wave strengths
        drho = rho_r - rho_l
        dp = p_r - p_l
        dv = v_r - v_l
        dv[0] = drho - dp / a_roe2
        dv[1] = dv + dp / (rho_roe * a_roe)
        dv[2] = dv - dp / (rho_roe * a_roe)
        # step:5 construct right characteristic eigenvectors
        rev = np.empty([3, 3])
        rev[0, :] = np.array([1, v_roe, 0.5 * v_roe * v_roe])
        rev[1, :] = np.array([1, v_roe + a_roe, h_roe + a_roe * v_roe])
        rev[2, :] = np.array([1, v_roe - a_roe, h_roe - a_roe * v_roe])
        # Step 7: compute flux
        fl = left
        for i in range(0, 3):
            fl += rev[i, :] * np.max(ev[i], 0) * dv[i]
            # add solution to flux

    def get_max_speed(self):
        self.get_sound_speed()
        max_speed = np.abs(self.unk[1, :]) + self.sound_speed
        self.dt = max_speed.max() * self.dx
        self.dtdx = 1. / max_speed.max()
        pass

    def get_primitive_variables(self):
        rho = self.unk[0, :]
        v1 = self.unk[1, :] / rho
        ke = 0.5 * rho * v1 * v1
        e_int = self.unk[2, :] - ke
        p = e_int * (self.gamma - 1)

    def get_flux(self):
        rho = self.unk[0, :]
        v1 = self.unk[1, :] / rho
        ke = 0.5 * rho * v1 * v1
        e_int = self.unk[2, :] - ke
        p = e_int * (self.gamma - 1)
        e_tot = ke + e_int
        self.flux[0, :] = rho * v1
        self.flux[1, :] = rho * v1 * v1 + p
        self.flux[2, :] = v1 * (p + e_tot)

    def boundary_conditions(self):
        self.unk_n[:, 0] = self.unk[:, 1]
        self.unk_n[:, -1] = self.unk[:, -2]

    def time_advance(self):
        self.get_flux()
        for i in range(1, self.grid.nx - 1):
            self.unk_n[:, i] = 0.5 * (self.unk[:, i - 1] + self.unk[:, i + 1]) - self.cfl * self.dtdx * (
                self.flux[:, i + 1] - self.flux[:, i - 1])

        self.boundary_conditions()
        self.unk = self.unk_n

    def main_loop(self):
        info_str = "a"
        self.set_initial_conditions()
        for self.step in range(0, self.steps):
            self.get_max_speed()
            info_str = "Timestep =" + str(self.step)
            info_str += ", cfl=" + str(self.dtdx)
            info_str += ", dt=%4.1e " % (self.dt)
            print(info_str)
            self.time_advance()
            self.plot_all()
        pass

    def plot_all(self):
        rho = self.unk[0, :]
        vel = self.unk[1, :] / rho
        p = self.unk[2, :]
        self.plotter.ax1.cla()
        pts = self.plotter.ax1.plot(rho)
        self.plotter.ax2.cla()
        pts = self.plotter.ax2.plot(vel)
        self.plotter.ax3.cla()
        pts = self.plotter.ax3.plot(p)
        plt.pause(0.05)
        self.plotter.fig.canvas.draw()
        plt.savefig('output' + str(self.step) + '.png')


if __name__ == '__main__':
    a = SimulationHydro()
