import unittest
from Vassal.Cell import Cell
import numpy.testing as np_tst
import numpy as np

class testCell(unittest.TestCase):
    cell = None

    # Create one instance per test.
    def setUp(self):
        self.cell = Cell()

    # Destroy instance as soon as test is over.
    def tearDown(self):
        self.cell = None

    def test_set_basis(self):
        # Test using angles and lattice parameters as input.
        self.cell.set_basis(lengths=[5.643, 6.621,4.885], angles=[91.83,
                            93.58, 107.69])
        self.assertAlmostEquals(173.30, self.cell.volume(), delta=1e-2)
        np_tst.assert_array_almost_equal([5.643, 6.621,4.885],
                                         self.cell.get_lattice_parameters())
        np_tst.assert_array_almost_equal([91.83, 93.58, 107.69],
                    self.cell.get_lattice_angles_radians(radians=False))

        # Simple test with a primitive cell.
        basis = np.zeros((3, 3))
        basis[0] = np.array([0, 2.986, 2.986])
        basis[1] = np.array([2.986, 0, 2.986])
        basis[2] = np.array([2.986, 2.986, 0])

        self.cell.set_basis(basis=basis)
        self.assertAlmostEquals(13.312*4, self.cell.volume(), delta=1e-3)
        np_tst.assert_array_almost_equal([4.223, 4.223, 4.223],
                                         self.cell.get_lattice_parameters(),
                                         decimal=3)
        np_tst.assert_array_almost_equal([60, 60, 60],
                    self.cell.get_lattice_angles_radians(radians=False))