import unittest

from .context import SimulationHydro
import numpy.testing as npt


class TestBasic(unittest.TestCase):
    def test_cons_to_primitive(self):
        a = SimulationHydro.SimulationHydro()
        testprim = [0, 0.4, 1.4]
        con = [1, 0, 1]
        prim = a.get_prim(con)
        print (prim)
        npt.assert_allclose(prim, testprim, rtol=1e-5)

    def test_riemann(self):
        a = SimulationHydro.SimulationHydro()
        left = [1, 0.0, 1.4]
        right = left
        testflux = [0, 0.56, 0]
        flux = a.riemann_solver(left, right)

        npt.assert_allclose(flux, testflux, rtol=1e-5)


if __name__ == '__main__':
    unittest.main()
