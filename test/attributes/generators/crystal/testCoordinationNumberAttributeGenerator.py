import unittest
from attributes.generators.crystal.CoordinationNumberAttributeGenerator import CoordinationNumberAttributeGenerator
from data.materials.AtomicStructureEntry import AtomicStructureEntry
from vassal.data.Atom import Atom
from vassal.data.Cell import Cell

class testCoordinationNumberAttributeGenerator(unittest.TestCase):
    def get_generator(self):
        return CoordinationNumberAttributeGenerator()
    
    def expected_count(self):
        return 4
    
    def test_results(self):
        # Create Dataset.
        dataset = []

        # Create primitive cell for B2-AlNi.
        structure1 = Cell()
        structure1.set_basis(lengths=[2.88, 2.88, 2.88], angles=[90, 90, 90])
        structure1.add_atom(Atom([0, 0, 0], 0))
        structure1.add_atom(Atom([0.5, 0.5, 0.5], 1))
        structure1.set_type_name(0, "Al")
        structure1.set_type_name(1, "Ni")
        entry1 = AtomicStructureEntry(structure1, name="Primitive", radii=None)
        dataset.append(entry1)

        # Create Scaled Cell.
        structure2 = Cell()
        structure2.set_basis(lengths=[3.0, 3.0, 3.0], angles=[90, 90, 90])
        structure2.add_atom(Atom([0, 0, 0], 0))
        structure2.add_atom(Atom([0.5, 0.5, 0.5], 1))
        structure2.set_type_name(0, "Al")
        structure2.set_type_name(1, "Ni")
        entry2 = AtomicStructureEntry(structure2, name="Scaled", radii=None)
        dataset.append(entry2)

        # Create a cell where A & B are swapped.
        structure3 = Cell()
        structure3.set_basis(lengths=[3.0, 3.0, 3.0], angles=[90, 90, 90])
        structure3.add_atom(Atom([0, 0, 0], 0))
        structure3.add_atom(Atom([0.5, 0.5, 0.5], 1))
        structure3.set_type_name(0, "Al")
        structure3.set_type_name(1, "Ni")
        entry3 = AtomicStructureEntry(structure3, name="Scaled", radii=None)
        dataset.append(entry3)

        # Create a 2x1x1 superce.
        structure4 = Cell()
        structure4.set_basis(lengths=[6.0, 3.0, 3.0], angles=[90, 90, 90])
        structure4.add_atom(Atom([0, 0, 0], 0))
        structure4.add_atom(Atom([0.5, 0, 0], 0))
        structure4.add_atom(Atom([0.25, 0.5, 0.5], 1))
        structure4.add_atom(Atom([0.75, 0.5, 0.5], 1))
        structure4.set_type_name(0, "Ni")
        structure4.set_type_name(1, "Al")
        entry4 = AtomicStructureEntry(structure4, name="Primitive", radii=None)
        dataset.append(entry4)

        # Generate features.
        gen = self.get_generator()
        features = gen.generate_features(dataset)

        # Make sure the correct number were generated.
        self.assertAlmostEquals(4, len(dataset))

        # Make sure scaling doesn't effect it.
        for i in range(len(dataset)):
            self.assertAlmostEquals(features.values[0][i], features.values[1][i], delta=1e-6)

        # Make sure its permutationally-invariant.
        for i in range(len(dataset)):
            self.assertAlmostEquals(features.values[0][i], features.values[2][i], delta=1e-6)

        # Make sure it passes supercell.
        for i in range(len(dataset)):
            self.assertAlmostEquals(features.values[0][i], features.values[3][i], delta=1e-6)
