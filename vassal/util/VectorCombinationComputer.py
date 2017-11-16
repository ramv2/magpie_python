from gmpy2 import mpfr

import gmpy2
import numpy as np
from numpy.linalg import norm
from scipy.linalg import lu_factor, lu_solve
from vassal.data.Cell import Cell

class VectorCombinationComputer:
    """
    Class to find all combinations of 3 vectors that are shorter than a
    certain length.
    """
    def __init__(self, in_vectors, cutoff_distance, include_zero=True):
        """
        Function to create the tool to compute all combinations of input
        vectors shorter than cutoff distance.
        :param in_vectors: Vectors to be combined. Must be exactly 3
        linearly-independent vectors.
        :param cutoff_distance: Desired cutoff distance.
        :param include_zero: Whether to include the zero vector in the list.
        """

        if len(in_vectors) != 3:
            raise ValueError("Expecting exactly three vectors.")

        # Vectors to be added.
        self.input_vectors = list(in_vectors)

        # Square of cutoff distance.
        self.cutoff_distance_sq = mpfr(cutoff_distance ** 2)

        # Whether to include the zero vector in the list.
        self.include_zero = include_zero

        # x, y, z coordinates of each vector shorter than cutoff.
        self.super_cells = []

        # All vectors shorter then cutoff.
        self.vectors = []

        self.get_all_vectors()

    def compute_vector(self, x):
        i_v = np.array([map(mpfr, iv) for iv in self.input_vectors])
        return np.array([x.dot(y) for y in i_v.T])


    def get_all_vectors(self):
        """
        Function to compute all vectors shorter than cutoff distance.
        :return:
        """

        # Create a matrix of basis vectors.
        basis = np.array([map(mpfr, x) for x in self.input_vectors],
                         dtype=object).T

        # Create ability to invert it.
        try:
            det_basis = np.linalg.det(basis)
        except TypeError:
            det_basis = Cell.get_determinant(basis)

        if det_basis == 0 or det_basis < 1e-14:
            raise RuntimeError("Vectors are not linearly independent.")

        m1, m2 = lu_factor(basis)
        m1_mpfr = np.array([map(mpfr, i) for i in m1], dtype=object)

        # Compute range of each variable.
        cutoff_distance = gmpy2.sqrt(self.cutoff_distance_sq)
        step_range = []

        for i in range(3):
            max_disp = 0.0
            for j in range(3):
                max_disp += np.dot(self.input_vectors[i], self.input_vectors[
                    j]) / Cell.get_mpfr_norm(self.input_vectors[i])
            step_range.append(int(gmpy2.ceil(max_disp / cutoff_distance)) + 1)

        # Ensure that we have sufficient range to get the cutoff distance
        # away from the origin by checking that we have large enough range to
        #  access a point cutoff distance away along the direction of xy,
        # xz and yz cross products.
        for dir in range(3):
            point = np.cross(self.input_vectors[dir], self.input_vectors[(dir
                                    + 1) % 3])
            point = point * cutoff_distance / Cell.get_mpfr_norm(point)
            sln = map(mpfr, lu_solve((m1_mpfr, m2), point))
            step_range = [max(step_range[i], int(gmpy2.ceil(abs(sln[i]))))
                          for i in range(3)]

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
        """
        Function to get the list of all vectors shorter than cutoff.
        :return: List of vectors.
        """
        return list(self.vectors)

    def get_supercell_coordinates(self):
        """
        Function to get the list of all image coordinates of vectors.
        :return: List of supercell coordinates.
        """
        return list(self.super_cells)