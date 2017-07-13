from __future__ import print_function

import unittest

import numpy.testing as npt

from .context import SimulationHydro


class TestBasic(unittest.TestCase):
    @staticmethod
    def test_cons_to_primitive():
        a = SimulationHydro.SimulationHydro()
        test_primitive = [0, 0.4, 1.4]
        con = [1, 0, 1]
        prim = a.get_prim(con)
        print(prim)
        npt.assert_allclose(prim, test_primitive, rtol=1e-5)

    @staticmethod
    def test_riemann():
        a = SimulationHydro.SimulationHydro()
        left = [1, 0.0, 1.4]
        right = left
        test_flux = [0, 0.56, 0]
        flux = a.riemann_solver(left, right)

        npt.assert_allclose(flux, test_flux, rtol=1e-5)


if __name__ == '__main__':
    unittest.main()
