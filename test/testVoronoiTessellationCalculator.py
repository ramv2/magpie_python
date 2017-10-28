import unittest
import numpy.testing as np_tst
from Vassal.Atom import Atom
from Vassal.Cell import Cell
from Vassal.VoronoiTessellationCalculator import VoronoiTessellationCalculator


class testVoronoiTessellationCalculator(unittest.TestCase):
    def test_simple_cubic(self):
        # Create the simulation cell.
        structure = Cell()
        atom = Atom([0, 0, 0], 0)
        structure.add_atom(atom)

        # Run tessellation.
        result = VoronoiTessellationCalculator.compute(structure, radical=False)

        # Test results.
        self.assertEquals(structure.n_atoms(), len(result))
        self.assertEquals(6, len(result[0].get_faces()))
        self.assertAlmostEquals(structure.volume(), result[0].get_volume(),
                                delta=1e-6)
        poly_index = result[0].get_polyhedron_shape()
        self.assertEquals(6, poly_index[4])
        poly_index = result[0].get_coordination_shell_shape(result)
        self.assertEquals(6, poly_index[0])

        # Test out the nearest neighbor shells.
        # 1st NN shell.
        nns = result[0].get_neighbor_shell(result, 1)
        self.assertEquals(6, len(nns))

        # 2nd NN shell.
        nns = result[0].get_neighbor_shell(result, 2)
        self.assertEquals(18, len(nns))

        # 3rd - 5th NN shell.
        for s in range(3, 6):
            nns = result[0].get_neighbor_shell(result, s)
            for image in nns:
                cell = image.get_supercell()
                self.assertAlmostEquals(s, sum([abs(cell[i]) for i in range(
                    3)]), delta=1e-6)

        # Test path NNs.
        # 0th NN.
        paths = result[0].get_neighbors_by_walks(result, 0)
        self.assertEquals(1, len(paths))
        total_weight = sum(paths.values())
        self.assertAlmostEquals(1.0, total_weight, delta=1e-6)

        # 1st NN.
        paths = result[0].get_neighbors_by_walks(result, 1)
        self.assertEquals(6, len(paths))
        total_weight = sum(paths.values())
        np_tst.assert_array_almost_equal([1/6.0]*6, paths.values())
        self.assertAlmostEquals(1.0, total_weight, delta=1e-6)