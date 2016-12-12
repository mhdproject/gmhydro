#!/usr/bin/env python3
from __future__ import print_function

import sys

import matplotlib
import numpy as np

matplotlib.use('macosx')

import matplotlib.pyplot as plt
import time

from ComputationalGrid import ComputationalGrid

from DataPlotter import DataPlotter


class SimulationHydro(object):
    def __init__(self):
        self.grid = ComputationalGrid()
        # self.grid.nx = 50
        self.dx = 1. / self.grid.nx
        self.steps = 4
        self.step = 0
        self.dtdx = 1.0
        self.nvar = 3
        self.cfl = 0.5
        self.gamma = 7. / 5.
        self.unk = np.ndarray(shape=(self.nvar, self.grid.nx), dtype=float)
        self.unk_n = np.empty_like(self.unk)
        self.flux = np.empty_like(self.unk)
        self.ev = np.ndarray(shape=(3), dtype=float)
        self.dv = np.empty_like(self.ev)
        self.fl = np.empty_like(self.ev)
        self.sound_speed = np.empty_like(self.unk)
        print("test")
        self.plotter = DataPlotter()

        self.main_loop()
        time.sleep(10)
        pass

    def set_initial_conditions(self):
        nx = self.grid.nx
        rho_l = 1.
        p_l = 1.
        v_l = 0
        rho_r = 0.125
        p_r = 0.1
        v_r = 0
        self.unk[0, 0:nx / 2] = rho_l
        self.unk[1, 0:nx / 2] = rho_l * v_l
        self.unk[2, 0:nx / 2] = 0.5 * rho_l * v_l ** 2 + p_l / (self.gamma - 1)
        self.unk[0, nx / 2:] = rho_r
        self.unk[1, nx / 2:] = rho_r * v_r
        self.unk[2, nx / 2:] = 0.5 * rho_r * v_r ** 2 + p_r / (self.gamma - 1)

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
        p_r = (self.gamma - 1) * (etot_r - (rho_r * v_r * v_r) / 2)
        h_r = (etot_r + p_r) / rho_r
        return v_r, p_r, h_r

    def riemann_solver(self, left_state, right_state):
        # the culbert b. laney way!
        # step 1: compute primitives
        rho_l = left_state[0]
        rho_r = right_state[0]
        v_l, p_l, h_l = self.get_prim(left_state)
        v_r, p_r, h_r = self.get_prim(right_state)

        # step :2 roe averages
        srl = np.sqrt(rho_l)
        srr = np.sqrt(rho_r)
        rho_roe = srl * srr
        v_roe = (srr * v_r + srl * v_l) / (srl + srr)
        h_roe = (srr * h_r + srl * h_l) / (srl + srr)
        a_roe = np.sqrt((self.gamma - 1) * (h_roe - 0.5 * v_roe ** 2))
        a_roe2 = a_roe * a_roe
        # step:3 compute eigenvalues
        self.ev[0] = v_roe
        self.ev[1] = v_roe + a_roe
        self.ev[2] = v_roe - a_roe
        # step:4 compute wave strengths
        drho = rho_r - rho_l
        dp = p_r - p_l
        dvel = v_r - v_l
        self.dv[0] = drho - dp / a_roe2
        self.dv[1] = dvel + dp / (rho_roe * a_roe)
        self.dv[2] = dvel - dp / (rho_roe * a_roe)
        # step:5 construct right characteristic eigenvectors
        rev = np.empty([3, 3])
        rev[0, :] = np.array([1, v_roe, 0.5 * v_roe ** 2])
        rev[1, :] = (rho_roe / 2 / a_roe) * np.array([1, v_roe + a_roe, h_roe + a_roe * v_roe])
        rev[2, :] = (-rho_roe / 2 / a_roe) * np.array([1, v_roe - a_roe, h_roe - a_roe * v_roe])
        # Step 7: compute flux
        self.fl = self.get_flux(left_state)
        #
        if (drho != 0):
            print('rho', rho_l, rho_r)
            print('p', p_l, p_r)
            print('v', v_l, v_r)
            print('a_roe', a_roe)
            print('rho_roe', rho_roe)
            print('v_roe', v_roe)
            print('h_roe', h_roe)
            print('ev', self.ev)
            print('eigenstrength', self.dv)
            print('drho', drho, 'dp', dp, 'a_roe2', a_roe2)
            print('rev', rev)
            print('self.fl', self.fl)
        for i in range(0, 3):
            self.fl += rev[i, :] * np.minimum(self.ev[i], 0) * self.dv[i]
            if (drho != 0):
                print(self.fl[0], rev[i, 0], self.ev[i], self.dv[i])
                # add solution to flux
        return self.fl

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

    def get_flux(self, u):
        rho = u[0]
        v1 = u[1] / rho
        ke = 0.5 * rho * v1 * v1
        e_int = u[2] - ke
        p = e_int * (self.gamma - 1)
        e_tot = ke + e_int
        fl = np.empty_like(self.ev)
        fl[0] = rho * v1
        fl[1] = rho * v1 * v1 + p
        fl[2] = v1 * (p + e_tot)
        return fl

    def boundary_conditions(self):
        self.unk_n[:, 0] = self.unk[:, 1]
        self.unk_n[:, -1] = self.unk[:, -2]

    def time_advance(self):
        # self.get_flux()
        for i in range(1, self.grid.nx - 1):
            leftstate = self.unk[:, i - 1]
            rightstate = self.unk[:, i]
            fl1 = self.riemann_solver(leftstate, rightstate)
            leftstate = self.unk[:, i]
            rightstate = self.unk[:, i + 1]
            fl2 = self.riemann_solver(leftstate, rightstate)
            # self.unk_n[:, i] = 0.5 * (self.unk[:, i - 1] + self.unk[:, i + 1]) - self.cfl * self.dtdx * (
            #     self.flux[:, i + 1] - self.flux[:, i - 1])

            self.unk_n[:, i] = self.unk[:, i] - self.cfl * self.dtdx * (
                fl2 - fl1)

            rho = self.unk_n[0, i]
            if (rho < 0):
                print("Negative density")
                sys.exit()

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
        rho = self.unk_n[0, :]
        vel = self.unk_n[1, :] / rho
        p = self.unk_n[2, :]
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
