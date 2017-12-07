from ..data.Atom import Atom
from ..data.Cell import Cell

class MPStructure:
    """
    Class to import a structure from a Materials Project structure's string
    representation. According to pymatgen's documentation, the coordinates of
    atoms used in MPStructure instance is of type fractional.
    """

    @classmethod
    def import_structure(self, lines):
        """
        Function to import a structure from a Materials Project structure's
        string representation.
        :param lines: String representation of a Materials Project structure
        object.
        :return: Structure as Cell.
        """

        # Initialize the cell.
        structure = Cell()

        lengths = map(float, lines[2].strip().split()[2:])
        angles = map(float, lines[3].strip().split()[1:])

        # Check if the lengths of lists are 3.
        if len(lengths) != 3 or len(angles) != 3:
            raise Exception("Expected a 3-dimensional vector.")

        # Set the basis.
        structure.set_basis(lengths=lengths, angles=angles)

        # Parse atoms and add them to structure.
        t = 0
        prev_name = ""
        for i in range(7, len(lines)):
            words = lines[i].strip().split()
            this_name = words[1]
            this_pos = map(float, words[2:5])
            if prev_name != "" and prev_name != this_name:
                structure.set_type_name(t, prev_name)
                t += 1
            atom = Atom(this_pos, t)
            structure.add_atom(atom)
            prev_name = this_name
        structure.set_type_name(t, prev_name)

        return structure