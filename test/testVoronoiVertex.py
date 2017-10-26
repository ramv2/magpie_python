import unittest
from Vassal.Atom import Atom
from Vassal.AtomImage import AtomImage
from Vassal.Cell import Cell

class testVoronoiVertex(unittest.TestCase):
    def test_creation(self):
        # Make a simple crystal.
        cell = Cell()
        cell.add_atom(Atom([0, 0, 0], 0))

        # Initialize faces.
        image = AtomImage(cell.get_atom(0), [0, 0, 1])