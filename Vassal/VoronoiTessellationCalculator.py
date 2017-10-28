from Vassal.PairDistanceAnalysis import PairDistanceAnalysis
from Vassal.VoronoiCell import VoronoiCell

class VoronoiTessellationCalculator:
    @classmethod
    def compute(self, cell, radical):
        # Initialize Voronoi cells.
        output = [VoronoiCell(atom, radical) for atom in cell.get_atoms()]

        # Create tool to find closest images.
        image_finder = PairDistanceAnalysis()
        atom_length_scale = (cell.volume()/ cell.n_atoms()) ** (1/3.0)

        image_finder.analyze_structure(cell)

        # Generate cells.
        for c in output:
            c.compute_cell(image_finder, cutoff=atom_length_scale * 6)

        # Read the.
        return output