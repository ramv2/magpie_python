import types
import math
import numpy as np
import pandas as pd
from CompositionDistanceFilter import CompositionDistanceFilter
from EqualSumCombinations import EqualSumCombinations
import sys

class APEAttributeGenerator:
    def __init__(self, lp):
        self.lp = lp
        self.packing_threshold = 0.01
        self.n_nearest_to_eval = [1, 3, 5]
        self.radius_property = "MiracleRadius"
        self.max_n_types = 6

    def set_packing_threshold(self, threshold):
        if threshold < 0:
            raise ValueError("Threshold must be positive.")
        self.packing_threshold = threshold

    def set_n_nearest_to_eval(self, values):
        self.n_nearest_to_eval = values

    def set_radius_property(self, prop):
        self.radius_property = prop

    def get_ideal_radius_ratio(self, n_neighbors):
        if n_neighbors <= 3:
            return 0.154701
        elif n_neighbors == 4:
            return 0.224745
        elif n_neighbors == 5:
            return 0.361654
        elif n_neighbors == 6:
            return 0.414213
        elif n_neighbors == 7:
            return 0.518145
        elif n_neighbors == 8:
            return 0.616517
        elif n_neighbors == 9:
            return 0.709914
        elif n_neighbors == 10:
            return 0.798907
        elif n_neighbors == 11:
            return 0.884003
        elif n_neighbors == 12:
            return 0.902113
        elif n_neighbors == 13:
            return 0.976006
        elif n_neighbors == 14:
            return 1.04733
        elif n_neighbors == 15:
            return 1.11632
        elif n_neighbors == 16:
            return 1.18318
        elif n_neighbors == 17:
            return 1.24810
        elif n_neighbors == 18:
            return 1.31123
        elif n_neighbors == 19:
            return 1.37271
        elif n_neighbors == 20:
            return 1.43267
        elif n_neighbors == 21:
            return 1.49119
        elif n_neighbors == 22:
            return 1.54840
        elif n_neighbors == 23:
            return 1.60436
        else:
            return 1.65915


    def compute_APE(self, n_neighbors=None, center_radius=None,
                    neigh_eff_radius=None,
                    radii=None, center_type=None, shell_types=None):

        n_n = n_neighbors
        c_r = center_radius
        n_e_r = neigh_eff_radius

        if n_n is None and c_r is None and n_e_r is None:
            # Get the radius of the central atom.
            c_r = radii[center_type]

            # Get the mean radius of the 1st neighbor shell.
            n_e_r = np.dot(shell_types, radii)
            n_n = sum(shell_types)

            # for i in xrange(len(shell_types)):
            #     n_neighbors += shell_types[i]
            #     neighbor_radius += shell_types[i] * radii[i]

            n_e_r /= n_n

        ideal_ratio = self.get_ideal_radius_ratio(n_n)
        actual_ratio = c_r / n_e_r
        return ideal_ratio / actual_ratio

    def get_cluster_range(self, radii, packing_threshold):

        # Get the maximum and minimum radius.
        biggest_radius = max(radii)
        smallest_radius = min(radii)

        # Compute the smallest possible cluster.
        cluster = np.zeros(len(radii))
        center_type = smallest_radius
        cluster[biggest_radius] = 3

        while self.compute_APE(radii=radii, center_type=center_type,
                               shell_types=cluster) < (1 - packing_threshold):
            cluster[biggest_radius] += 1
            if cluster[biggest_radius] > 24:
                raise RuntimeError("Smallest cluster > 24 atoms: "
                                   "packing_threshold must be too large")

        smallest_cluster = cluster[biggest_radius]

        # Compute the largest possible cluster.
        cluster[biggest_radius] = 0
        cluster[smallest_radius] = 24
        center_type = biggest_radius

        while self.compute_APE(radii=radii, center_type=center_type,
                               shell_types=cluster) < (1 + packing_threshold):
            cluster[smallest_radius] -= 1
            if cluster[smallest_radius] < 3:
                raise RuntimeError("Largest cluster < 3 atoms: "
                                   "packing_threshold must be too large")

        biggest_cluster = cluster[smallest_radius]

        return smallest_cluster, biggest_cluster

    def get_closest_compositions(self, target_composition,
                                 other_compositions, n_closest, p_norm):

        tmp_list = []
        for comp in other_compositions:
            cdf = CompositionDistanceFilter.compute_distance(comp,
                  target_composition, p_norm)
            tmp_list.append([cdf, comp])

        tmp_list.sort()
        values = []
        for i in xrange(n_closest+1):
            values.append(tmp_list[i][1])

        return values

    def find_efficiently_packed_clusters(self, radii, packing_threshold):
        output = []

        # Special case: Only one atom type.
        if len(radii) == 1:
            clusters =[]

            # Loop through all cluster sizes.
            for i in xrange(3, 24):
                ape = self.compute_APE(i, center_radius=1.0,
                                       neigh_eff_radius=1.0)
                if abs(ape - 1) < 0.05:
                    clusters.append([i])

            output.append(clusters)
            return output

        # Determine the minimum and maximum cluster sizes.
        min_cluster_size, max_cluster_size = self.get_cluster_range(radii,
                                            packing_threshold)

        # Loop through each atom as the central type.
        for central_type in xrange(len(radii)):
            clusters = []

            # Loop over possible ranges of cluster sizes (determined from radii)
            for cluster_size in xrange(min_cluster_size, max_cluster_size+1):
                esc = EqualSumCombinations(cluster_size, len(radii))

                # Loop through all combinations of atom types in the first
                # shell
                for shell in esc.get_combinations():
                    ape = self.compute_APE(radii=radii,
                                           center_type=central_type,
                                           shell_types=shell)

                    if abs(ape - 1) < packing_threshold:
                        clusters.append(list(shell))

            output.append(clusters)

        return output

    def compute_cluster_compositions(self, elements, clusters):
        output = []

        # Array to store the number of each atom type.
        fractions = np.zeros(len(elements))

        # Loop through clusters with each type of atom at the center.
        for ct in xrange(len(clusters)):
            for shell in clusters[ct]:
                for i in xrange(len(fractions)):
                    fractions[i] = shell[i]

                fractions[ct] += 1

            entry = {}
            for i in xrange(len(elements)):
                e_name = self.lp.element_names[elements[i]]
                entry[e_name] = fractions[i]

            output.append(entry)

        return output

    def determine_optimal_APE(self, central_atom_type, shell_composition,
                              radii):
        # Initialize output.
        output = float("inf")

        # Get radius of center, mean radius of outside.
        center_r = radii[central_atom_type]
        shell_r = np.average(radii, weights=shell_composition.values())

        # Loop through all atom sizes.
        for z in xrange(3, 24):
            ape = self.compute_APE(n_neighbors=z,
                                   center_radius=center_r,
                                   neigh_eff_radius=shell_r)
            if abs(ape - 1) < abs(output - 1):
                output = ape
            else:
                break

        return output

    def generate_features(self, entries, verbose=False):
        # Initialize lists of feature values and headers for pandas data frame.
        feat_values = []
        feat_headers = []

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Get the atomic radii.
        radii_lookup = self.lp.load_property(self.radius_property)

        # Insert header names here.
        for n in self.n_nearest_to_eval:
            feat_headers.append("APE_Nearest{}_Below{}".format(n,
                                        self.packing_threshold))

        feat_headers.append("APE_SystemAverage")
        feat_headers.append("APE_SystemAverageDeviation")

        # Find the largest number of clusters to be considered.
        largest_n = max(self.n_nearest_to_eval)

        # Get the entries, sort so that alloys for the same system are
        # grouped together.