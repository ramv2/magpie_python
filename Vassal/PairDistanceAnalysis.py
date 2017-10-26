import numpy as np
from numpy.linalg import norm
from sympy import Plane
from Vassal.AtomImage import AtomImage
from Vassal.util.VectorCombinationComputer import VectorCombinationComputer
import math

class PairDistanceAnalysis:
    def __init__(self):
        # Cutoff distance. No pairs farther than this value are considered.
        self.cutoff_distance = 6.5

        # Periodic images to consider when searching for neighbors.
        self.supercells = []

        # Lattice vectors corresponding to each image.
        self.lattice_vectors = []

        # Link to structure being evaluated.
        self.structure = None

    def precompute(self):
        lat_vectors = self.structure.get_lattice_vectors()
        p0 = Plane(lat_vectors[0] * 0.5, normal_vector=lat_vectors[0])
        p1 = Plane(lat_vectors[1] * 0.5, normal_vector=lat_vectors[1])
        p2 = Plane(lat_vectors[2] * 0.5, normal_vector=lat_vectors[2])
        max_image_dist = norm(np.array(p0.intersection(p1)[0].intersection(
            p2)[0].evalf(), dtype=float))

        computer = VectorCombinationComputer(lat_vectors, max_image_dist +
                                             self.cutoff_distance)
        self.supercells = computer.get_supercell_coordinates()
        self.lattice_vectors = computer.get_vectors()

    def get_cutoff_distance(self):
        return self.cutoff_distance

    def set_cutoff_distance(self, d):
        self.cutoff_distance = d
        if self.structure is not None:
            self.precompute()

    def find_all_images(self, center_atom, neighbor_atom):
        # Get the two atoms.
        center_pos = self.structure.get_atom(
            center_atom).get_position_cartesian()

        # Find all images.
        output = []
        closest_image = self.structure.get_minimum_distance(
            center=center_atom, neighbor=neighbor_atom)

        B_sub_A = closest_image.get_position()
        B_sub_A -= center_pos

        closest_supercell = closest_image.get_supercell()

        # Loop over each periodic image to find those within range.
        cutoff_distance_sq = self.cutoff_distance ** 2
        for img in range(len(self.supercells)):
            supercell_vec = self.lattice_vectors[img]
            new_pos = B_sub_A.copy() + supercell_vec
            dist = new_pos[0] ** 2 + new_pos[1] ** 2 + new_pos[2] ** 2

            if dist < cutoff_distance_sq and dist > 1e-8:
                ss = self.supercells[img].copy() + closest_supercell
                output.append((AtomImage(closest_image.get_atom(), ss),
                               math.sqrt(dist)))

        return output

    def get_all_neighbors_of_atom(self, index):
        output = []
        for atom in self.structure.get_atoms():
            output += self.find_all_images(index, atom.get_id())
        return output

    def compute_PRDF(self, n_bin):
        if n_bin <= 0:
            raise ValueError("Number of bins must be 0.")

        n_t = self.structure.n_types()
        n_a = self.structure.n_atoms()

        # Initialize arrays.
        output = np.zeros((n_t, n_t, n_bin))

        # Find all pairs within the cutoff distance.
        n_type = np.zeros(n_t)

        for i in range(n_a):
            for j in range(i, n_a):
                # Find images.
                images = self.find_all_images(i, j)
                i_type = self.structure.get_atom(i).get_type()
                n_type[i_type] += 1
                j_type = self.structure.get_atom(j).get_type()

                # For each image, assign it to bin.
                for img in images:
                    bin_ = int(math.floor(img[1] * n_bin / self.cutoff_distance))
                    if bin_ >= n_bin:
                        # Happens if dist equals cutoff.
                        continue
                    output[i_type][j_type][bin_] += 1
                    output[j_type][i_type][bin_] += 1

        # Normalizing data.
        bin_spacing = self.cutoff_distance / n_bin
        for b in range(n_bin):
            vol = 4.0 / 3.0 * ((b + 1) ** 3 - b ** 3) * bin_spacing ** 3 * \
                  math.pi
            for i in range(n_t):
                for j in range(n_t):
                    output[i][j][b] /= vol * n_type[i]

        return output

    def analyze_structure(self, s):
        self.structure = s
        self.precompute()

    def recompute(self):
        self.precompute()