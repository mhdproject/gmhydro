#!/usr/bin/env python3
from __future__ import print_function

import sys

import matplotlib
import numpy as np

matplotlib.use('macosx')

import time

from ComputationalGrid import ComputationalGrid

from DataPlotter import DataPlotter


class SimulationHydro(object):
    def __init__(self):
        self.grid = ComputationalGrid()
        # self.grid.nx = 50
        self.dx = 1. / self.grid.nx
        self.steps = 4000
        self.step = 0
        self.dtdx = 1.0
        self.t = 0
        self.tend = 0.2
        self.dt = 1.0
        self.nvar = 3
        self.cfl = 0.49
        self.gamma = 7. / 5.
        self.unk = np.ndarray(shape=(self.nvar, self.grid.nx), dtype=float)
        self.unk_n = np.empty_like(self.unk)
        self.flux = np.empty_like(self.unk)
        self.ev = np.ndarray(shape=3, dtype=float)
        self.dv = np.empty_like(self.ev)
        self.fl = np.empty_like(self.ev)
        self.debug = 0
        self.sound_speed = np.empty_like(self.unk)
        print("test")
        self.plotter = DataPlotter()

        # self.main_loop()
        time.sleep(10)
        pass

    def set_initial_conditions(self):
        nx = self.grid.nx
        rho_l = 1.
        p_l = 1.
        u_l = 0
        rho_r = 0.125
        p_r = 0.1
        u_r = 0
        self.unk[0, 0:nx / 2] = rho_l
        self.unk[1, 0:nx / 2] = rho_l * u_l
        self.unk[2, 0:nx / 2] = 0.5 * rho_l * u_l ** 2 + p_l / (self.gamma - 1)
        self.unk[0, nx / 2:] = rho_r
        self.unk[1, nx / 2:] = rho_r * u_r
        self.unk[2, nx / 2:] = 0.5 * rho_r * u_r ** 2 + p_r / (self.gamma - 1)

    def get_sound_speed(self):
        rho = self.unk[0, :]
        vx = self.unk[1, :] / rho
        ke = 0.5 * rho * vx * vx
        e_tot = self.unk[2, :]
        e_int = e_tot - ke
        p = e_int * (self.gamma - 1)
        self.sound_speed = np.sqrt(self.gamma * p / rho)

    def get_prim(self, con):
        rho = con[0]
        mom = con[1]
        etot = con[2]
        v = mom / rho
        ke = 0.5 * rho * v ** 2
        p = (self.gamma - 1) * (etot - ke)
        h = (etot + p) / rho
        return v, p, h

    def get_flux(self, u):
        v, p, h = self.get_prim(u)
        rho = u[0]
        e_tot = u[2]

        fl = np.empty_like(self.ev)
        fl[0] = rho * v
        fl[1] = rho * v ** 2 + p
        fl[2] = v * (p + e_tot)
        return fl

    def riemann_solver(self, left_state, right_state):
        # the culbert b. laney way!
        # step 1: compute primitives
        rho_l = left_state[0]
        rho_r = right_state[0]
        u_l, p_l, h_l = self.get_prim(left_state)
        u_r, p_r, h_r = self.get_prim(right_state)

        # step :2 roe averages
        srl = np.sqrt(rho_l)
        srr = np.sqrt(rho_r)
        rho_RL = srl * srr
        u_RL = (srr * u_r + srl * u_l) / (srl + srr)
        h_RL = (srr * h_r + srl * h_l) / (srl + srr)
        a_RL2 = (self.gamma - 1) * (h_RL - 0.5 * u_RL ** 2)
        a_RL = np.sqrt(a_RL2)
        # step:3 compute eigenvalues
        self.ev[0] = u_RL
        self.ev[1] = u_RL + a_RL
        self.ev[2] = u_RL - a_RL
        # step:4 compute wave strengths
        drho = rho_r - rho_l
        dp = p_r - p_l
        dvel = u_r - u_l
        self.dv[0] = drho - dp / a_RL2
        self.dv[1] = dvel + dp / (rho_RL * a_RL)
        self.dv[2] = dvel - dp / (rho_RL * a_RL)
        # step:5 construct right characteristic eigenvectors
        rev = np.empty([3, 3])
        rev[0, :] = np.array([1, u_RL, 0.5 * u_RL ** 2])
        rev[1, :] = (rho_RL / 2 / a_RL) * np.array([1, u_RL + a_RL, h_RL + a_RL * u_RL])
        rev[2, :] = (-rho_RL / 2 / a_RL) * np.array([1, u_RL - a_RL, h_RL - a_RL * u_RL])
        # Step 7: compute flux
        self.fl = 0.5 * (self.get_flux(left_state) + self.get_flux(right_state))
        #
        if (self.debug == 1):
            print('rho', rho_l, rho_r)
            print('p', p_l, p_r)
            print('v', u_l, u_r)
            print('a_RL', a_RL)
            print('rho_RL', rho_RL)
            print('u_RL', u_RL)
            print('h_RL', h_RL)
            print('ev', self.ev)
            print('eigenstrength', self.dv)
            print('drho', drho, 'dp', dp, 'a_RL2', a_RL2)
            print('rev', rev)
            print('self.fl', self.fl)
        for i in range(0, 3):
            self.fl -= 0.5 * rev[i, :] * np.abs(self.ev[i]) * self.dv[i]
            if (self.debug == 1):
                print(self.fl[0], rev[i, 0], self.ev[i], self.dv[i])
                # add solution to flux
        return self.fl

    def get_max_speed(self):
        self.get_sound_speed()
        max_speed = np.abs(self.unk[1, :]) + self.sound_speed
        self.dt = max_speed.max() * self.dx
        self.dtdx = 1. / max_speed.max()
        pass

    def boundary_conditions(self):
        self.unk_n[:, 0] = self.unk[:, 1]
        self.unk_n[:, -1] = self.unk[:, -2]

    def time_advance(self):
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
            info_str += ", cfl=%5.3f" % self.dtdx
            info_str += ", dt=%4.1e " % self.dt
            info_str += ", t=%4.1e " % self.t
            print(info_str)
            self.t += self.dt
            if (self.t > self.tend):
                sys.exit("End of simulation")
            self.time_advance()
            self.plotter.plot_all(self.unk_n, self.grid.x, self.step)
        pass


if __name__ == '__main__':
    a = SimulationHydro()
    a.main_loop()
