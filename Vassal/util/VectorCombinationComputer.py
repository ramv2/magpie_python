import numpy as np
from scipy.linalg import lu_factor, lu_solve
from numpy.linalg import det, norm
from math import sqrt, ceil, floor

class VectorCombinationComputer:
    def __init__(self, in_vectors, cutoff_distance, include_zero=True):
        if len(in_vectors) != 3:
            raise ValueError("Expecting exactly three vectors.")
        self.input_vectors = list(in_vectors)
        self.cutoff_distance_sq = cutoff_distance ** 2
        self.include_zero = include_zero
        self.super_cells = []
        self.vectors = []
        self.get_all_vectors()

    def compute_vector(self, x):
        return np.matmul(x, self.input_vectors)

    def get_all_vectors(self):
        # Create a matrix of basis vectors.
        basis = np.transpose(np.array(self.input_vectors))

        # Create ability to invert it.
        det_basis = det(basis)
        if det_basis == 0 or det_basis < 1e-14:
            raise RuntimeError("Vectors are not linearly independent.")

        fac = lu_factor(basis)

        # Compute range of each variable.
        cutoff_distance = sqrt(self.cutoff_distance_sq)
        step_range = []

        for i in range(3):
            max_disp = 0.0
            for j in range(3):
                max_disp += np.dot(self.input_vectors[i], self.input_vectors[
                    j]) / norm(self.input_vectors[i])
            step_range.append(int(ceil(max_disp / cutoff_distance)) + 1)

        # Ensure that we have sufficient range to get the cutoff distance
        # away from the origin by checking that we have large enough range to
        #  access a point cutoff distance away along the direction of xy,
        # xz and yz cross products.
        for dir in range(3):
            point = np.cross(self.input_vectors[dir], self.input_vectors[(dir
                                    + 1) % 3])
            point *= cutoff_distance / norm(point)
            sln = lu_solve(fac, point)
            step_range = [max(step_range[i], int(ceil(abs(sln[i])))) for i in
                          range(3)]

        # Create the initial vector.
        for x in range(-step_range[0], 1 + step_range[0]):
            for y in range(-step_range[1], 1 + step_range[1]):
                for z in range(-step_range[2], 1 + step_range[2]):
                    a = np.array([x, y, z])
                    l = self.compute_vector(a)
                    dist_sq = l[0] ** 2 + l[1] ** 2 +  l[2] ** 2
                    if dist_sq <= self.cutoff_distance_sq:
                        if not self.include_zero and x == 0 and y == 0 and z \
                                == 0:
                            continue
                        self.super_cells.append(a)
                        self.vectors.append(l)

    def get_vectors(self):
        return list(self.vectors)

    def get_supercell_coordinates(self):
        return list(self.super_cells)