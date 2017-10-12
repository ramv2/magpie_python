import types
from heapq import heappush, heappop
import numpy as np
import pandas as pd
from CompositionDistanceFilter import CompositionDistanceFilter
from CompositionEntry import CompositionEntry
from EqualSumCombinations import EqualSumCombinations
from LookUpData import LookUpData

class APEAttributeGenerator:
    """
    Class to compute features using Atomic Packing Efficiency (APE) of nearby
    clusters.The APE, as defined by Laws et al.
    (http://www.nature.com/doifinder/10.1038/ncomms9123), is determined based
    on the ideal and actual ratio between the central and shell atoms of an
    atomic cluster with a certain number of atoms. Often, the packing
    efficiency is described as a ratio between these two quantities:

    [Packing efficiency] = [Ideal radius ratio] / [Actual radius ratio]

    The ideal ratio is determined based on the ratio between the size of a
    central atom and the neighboring atoms such that the packing around the
    central atom is maximized. These optimal ratios for clusters for
    different numbers of atoms have been tabulated by Miracle et al.
    (http://dx.doi.org/10.2320/matertrans.47.1737)

    The actual ratio is computed by dividing the radius of the central atom
    by the average of the central atoms.

    We currently use this framework to create two types of features:

    Distance to nearest clusters with a packing efficiency better than a
    certain threshold. If there are fewer than a requested number of
    efficiently packed clusters in an alloy system, the average is taken to
    be the average distance to all of the clusters. These features are
    designed to measure the availability of efficiently-packed atomic
    configurations in the liquid.

    Mean packing efficiency of the system assuming that the composition of
    the first nearest neighbor shell is equal to the composition of the
    system. Each atom type is surrounded by the number of atoms that
    maximizes the packing efficiency. As shown in recent work by Laws et al.
    (http://www.nature.com/doifinder/10.1038/ncomms9123), bulk metallic
    glasses are known to form when the clusters around all types of atom
    have the same composition as the alloy and are efficiently packed. We
    compute the average APE for each atom in the system, under this
    assumption, and the average deviation from perfect packing.

    This algorithm currently evaluates all possible clusters provided a list
    of elements. As the number of clusters scales with N!, the runtime of
    this algorithm scales with N!.

    For now, we only search for clusters with up to 7 atoms (max_n_types) in
    order to avoid this combinatorial problem. In practice, the algorithm
    picks the top 7 alloys with the highest fractions. While not idea,
    this might work in practice since most alloys have fewer than 7 main
    components. Many alloys >10 elements in the specification, but many are
    impurities that may not be present in large enough amounts to really
    affect the determination of efficiently packed clusters.

    """

    # Threshold at which to define a cluster as efficiently packed.
    # Packing efficiency is defined by |APE - 1|. Default value for this
    # parameter is 0.01.
    packing_threshold = 0.01

    # Number of nearest clusters to assess.
    n_nearest_to_eval = [1, 3, 5]

    # Maximum number of types over which to search for clusters. If an
    # alloy has more than this number of elements, the code will only
    # search for clusters with the most prevalent elements.
    max_n_types = 6

    # Name of elemental property to use as atomic radius. By default,
    # uses the radii from doi: 10.1179/095066010X12646898728200.
    radius_property = "MiracleRadius"

    def set_packing_threshold(self, threshold):
        """
        Function to define the threshold at which a cluster is considered
        efficiently packed.
        :param threshold: Desired threshold. Default: 0.01
        :return:
        """
        if threshold < 0:
            raise ValueError("Threshold must be positive.")
        self.packing_threshold = threshold

    def set_n_nearest_to_eval(self, values):
        """
        Function to define the number of nearest neighbor clusters to
        evaluate when computing features.
        :param values: A list containing the numbers of neighbor clusters.
        Default: [1, 3, 5]
        :return:
        """
        self.n_nearest_to_eval = values

    def set_radius_property(self, prop):
        """
        Function to set the name of the elemental property used to define
        radii. By default uses the "MiracleRadius" property which is from an
        assessment by Miracle et al.
        (http://dx.doi.org/10.2320/matertrans.47.1737)
        :param prop: Name of property used to define radii.
        :return:
        """
        self.radius_property = prop

    @classmethod
    def compute_APE(self, n_neighbors=None, center_radius=None,
                    neigh_eff_radius=None,
                    radii=None, center_type=None, shell_types=None):
        """
        Function to compute the APE of a cluster, given the identities of the
        central and 1st neighbor atoms or just the number of neighbors and
        the radii. Here, we follow the formulation given by Laws et al.
        (http://www.nature.com/doifinder/10.1038/ncomms9123).

        APE = ideal radius ratio / (radius of central atom / effective radius
        of nearest neighbors)

        :param n_neighbors: Number of 1st nearest neighbors in the cluster.
        :param center_radius: Radius of the central atom.
        :param neigh_eff_radius: Effective radius of the 1st shell. Usually
        computed as the average radius of all atoms in the shell.
        :param radii: Radius of each atom type.
        :param center_type: Type of atom in the center.
        :param shell_types: Number of atoms of each type in the outer shell.
        Must be same length as radii.
        :return: APE, as defined above, as a float.
        """
        n_n = n_neighbors
        c_r = center_radius
        n_e_r = neigh_eff_radius

        if n_n is None and c_r is None and n_e_r is None:
            # Get the radius of the central atom.
            c_r = radii[center_type]

            # Get the mean radius of the 1st neighbor shell.
            # n_e_r = np.dot(shell_types, radii)
            # n_n = np.sum(shell_types)
            n_e_r = 0.0
            n_n = 0.0
            for i in xrange(len(shell_types)):
                n_e_r += shell_types[i] * radii[i]
                n_n += shell_types[i]


            # for i in xrange(len(shell_types)):
            #     n_neighbors += shell_types[i]
            #     neighbor_radius += shell_types[i] * radii[i]

            n_e_r /= n_n

        #  The ideal radius ratio is only known for clusters with 3 and 24 (
        # inclusive) neighbors. If you request outside of this range,
        # the value of 3 is set for anything less than 3 and the value of 24
        # is set for anything larger than 24.
        if n_n <= 3:
            ideal_ratio = 0.154701
        elif n_n == 4:
            ideal_ratio = 0.224745
        elif n_n == 5:
            ideal_ratio = 0.361654
        elif n_n == 6:
            ideal_ratio = 0.414213
        elif n_n == 7:
            ideal_ratio = 0.518145
        elif n_n == 8:
            ideal_ratio = 0.616517
        elif n_n == 9:
            ideal_ratio = 0.709914
        elif n_n == 10:
            ideal_ratio = 0.798907
        elif n_n == 11:
            ideal_ratio = 0.884003
        elif n_n == 12:
            ideal_ratio = 0.902113
        elif n_n == 13:
            ideal_ratio = 0.976006
        elif n_n == 14:
            ideal_ratio = 1.04733
        elif n_n == 15:
            ideal_ratio = 1.11632
        elif n_n == 16:
            ideal_ratio = 1.18318
        elif n_n == 17:
            ideal_ratio = 1.24810
        elif n_n == 18:
            ideal_ratio = 1.31123
        elif n_n == 19:
            ideal_ratio = 1.37271
        elif n_n == 20:
            ideal_ratio = 1.43267
        elif n_n == 21:
            ideal_ratio = 1.49119
        elif n_n == 22:
            ideal_ratio = 1.54840
        elif n_n == 23:
            ideal_ratio = 1.60436
        else:
            ideal_ratio = 1.65915

        actual_ratio = c_r / n_e_r
        return ideal_ratio / actual_ratio

    @classmethod
    def get_cluster_range(self, radii, packing_threshold):
        """
        Function compute the maximum and minimum possible cluster sizes,
        given a list of radii. The smallest possible cluster has the smallest
        atom in the center and the largest in the outside. The largest
        possible has the largest in the inside and the smallest in the outside.

        :param radii: List of radii of elements in the system.
        :param packing_threshold: APE defining maximum packing threshold.
        :return: min_cluster_size, max_cluster_size as defined by number of
        atoms in the shell.
        """
        l_r = len(radii)

        # Get the indices of the maximum and minimum radius.
        biggest_radius = np.argmax(radii)
        smallest_radius = np.nanargmin(radii)

        # Compute the smallest possible cluster.
        cluster = np.zeros(l_r)
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
                               shell_types=cluster) > (1 + packing_threshold):
            cluster[smallest_radius] -= 1
            if cluster[smallest_radius] < 3:
                raise RuntimeError("Largest cluster < 3 atoms: "
                                   "packing_threshold must be too large")

        biggest_cluster = cluster[smallest_radius]

        return int(smallest_cluster), int(biggest_cluster)

    def get_closest_compositions(self, target_composition,
                                 other_compositions, n_closest, p_norm):
        """
        Function to get closest compositions from a given target composition.

        :param target_composition: Composition from which to measure
        distance.
        :param other_compositions: A list of compositions whose distance from
        the target will be ranked.
        :param n_closest: Number of closest compounds to return.
        :param p_norm: P-norm to use when computing distance.
        :return: dist, comps containing a list of distances and compositions.
        """
        tmp_list = []
        for comp in other_compositions:
            cdf = -CompositionDistanceFilter.compute_distance(comp,
                  target_composition, p_norm)
            heappush(tmp_list, (cdf, comp))
            if len(tmp_list) > n_closest:
                heappop(tmp_list)

        comps = []
        dist = []
        while tmp_list:
            tup = heappop(tmp_list)
            dist.insert(0, -tup[0])
            comps.insert(0, tup[1])
        return dist, comps

    @classmethod
    def find_efficiently_packed_clusters(self, radii, packing_threshold):
        """
        Function to find all clusters with better APE than a certain
        threshold, given a list of atomic radii. The packing efficiency is
        defined as abs(1 - APE).

        :param radii: Radii of elements.
        :param packing_threshold: Desired packing limit threshold. A
        "default" choice would be 0.05.
        :return: A list of efficiently packed structures for each atom type
        as the central atom. Ex: x[0][1] is the 2nd efficiently packed
        cluster with atom type 0 as the central atom.
        """
        output = []
        l_r = len(radii)

        # Special case: Only one atom type.
        if l_r == 1:
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

        esc = EqualSumCombinations(max_cluster_size - 1, l_r)
        for cluster_size in xrange(min_cluster_size,
                                   max_cluster_size):
            if not esc.dp[cluster_size][l_r]:
                shells = esc.get_combinations(cluster_size, l_r)

        # Loop through each atom as the central type.
        for central_type in xrange(l_r):
            clusters = []

            # Loop over possible ranges of cluster sizes (determined from
            # radii).
            for cluster_size in xrange(min_cluster_size, max_cluster_size):

                # Loop through all combinations of atom types in the first
                # shell.
                for shell in esc.dp[cluster_size][l_r]:
                    ape = self.compute_APE(radii=radii,
                                           center_type=central_type,
                                           shell_types=shell)

                    if abs(ape - 1) < packing_threshold:
                        clusters.append(list(shell))

            output.append(clusters)

        return output

    @classmethod
    def compute_cluster_compositions(self, e_ids, clusters):
        """
        Function to compute the compositions of a list of atomic clusters.
        The composition includes both atoms in the first nearest neighbor
        shell and the atom in the center of the cluster.

        :param e_ids: Ids of the elements from which clusters are composed.
        :param clusters:  Clusters to convert. List of identity shell
        compositions for each type of central atom. Ex: clusters[1][2] is an
        array defining the number of atoms of each type for clusters with an
        atom of type 1 in the center.
        :return:
        """
        output = []
        l_c = len(clusters)

        # Loop through clusters with each type of atom at the center.
        for ct in xrange(l_c):
            for shell in clusters[ct]:
                fractions = list(shell)
                fractions[ct] += 1.0

                entry = CompositionEntry(element_ids=e_ids, fractions=fractions)
                output.append(entry)

        return output

    @classmethod
    def determine_optimal_APE(self, central_atom_type, shell_composition,
                              radii):
        """
        Function to compute the optimal APE for a cluster with a certain atom
        type in the center and composition in the cell. This algorithm finds
        the number of atoms in the shell such that the APE of the cluster is
        closest to 1. Note: This calculation assumes that sites in the first
        nearest-neighbor shell can be partially-occupied.

        :param central_atom_type: Element id (Z - 1) of th central atom.
        :param shell_composition: Composition of the nearest-neighbor shell.
        :param radii: Lookup table of elemental radii.
        :return: output: The optimal APE.
        """
        # Initialize output.
        output = float("inf")

        # Get radius of center, mean radius of outside.
        center_r = radii[central_atom_type]
        tmp_r = [radii[elem] for elem in shell_composition.get_element_ids()]
        shell_r = np.average(tmp_r,
                             weights=shell_composition.get_element_fractions())

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

    def generate_features(self, entries, lookup_path, verbose=False):
        """
        Function to generate features as mentioned in the class description.

        :param entries: A list of CompositionEntry's.
        :param verbose: Flag that is mainly used for debugging. Prints out a
        lot of information to the screen.
        :return features: Pandas data frame containing the names and values
        of the descriptors.
        """

        # Initialize lists of feature values and headers for pandas data frame.
        feat_values = []
        feat_headers = []

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of "
                             "CompositionEntry's")
        elif (entries and not isinstance(entries[0], CompositionEntry)):
            raise ValueError("Argument should be of type list of "
                             "CompositionEntry's")

        # Get the atomic radii.
        radii_lookup = LookUpData.load_property(self.radius_property,
                                                lookup_dir=lookup_path)

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
        entries.sort()

        # Compute features for each entry.
        last_elements = []
        clusters = []
        for entry in entries:
            # print entry
            tmp_list = []
            cur_elements = list(entry.get_element_ids())
            cur_fractions = list(entry.get_element_fractions())

            # If list of elements is greater than max_n_types, pick only the
            # most prevalent.
            if len(cur_elements) > self.max_n_types:
                cur_elements = [x for _,x in sorted(zip(cur_fractions,
                                        cur_elements), reverse=True)]
                cur_elements = cur_elements[:self.max_n_types]
                # cur_elements = sorted(cur_elements, key=entry.get)[
                #                :self.max_n_types]

            # Sort elements by atomic number.
            cur_elements.sort()

            # Get radii of those elements.
            radii = [radii_lookup[e_id] for e_id in cur_elements]

            # If current element list doesn't equal last_elements, recompute
            # list of nearest clusters.
            if cur_elements != last_elements:
                tmp = self.find_efficiently_packed_clusters(radii,
                                                        self.packing_threshold)
                clusters = self.compute_cluster_compositions(cur_elements, tmp)
                last_elements = cur_elements

            # Find the closest clusters and distance to our cluster.
            distances, closest_clusters = self.get_closest_compositions(entry,
                                            clusters, largest_n, 2)
            l_d = len(distances)
            if l_d == 0:
                for n in self.n_nearest_to_eval:
                    tmp_list.append(1000.0)
            else:
                for n in self.n_nearest_to_eval:
                    tmp_list.append(np.mean(distances[:min(n, l_d - 1)]))

            # Compute the packing efficiency of clusters around each atom
            # assuming that the composition of the first nearest-neighbor
            # shell is equal to the composition of the alloy.

            entry_elements = list(entry.get_element_ids())
            entry_fractions = list(entry.get_element_fractions())
            cluster_APEs = [self.determine_optimal_APE(id, entry,
                                    radii_lookup) for id in entry_elements]

            # Compute the composition-weighted average and average deviation
            # from 1.
            avg = np.average(cluster_APEs, weights=entry_fractions)
            avg_dev = np.average([abs(1.0 - c) for c in cluster_APEs],
                                 weights=entry_fractions)
            tmp_list.append(avg)
            tmp_list.append(avg_dev)

            feat_values.append(tmp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    # entries = [{'Fe': 0.4, 'Cn': 0.6}, {'H': 0.6666666666666666,
    #             'O': 0.3333333333333333}, {'Na': 0.5,'Cl': 0.5},
    #            {"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}]
    entries = [{"Sc": 0.25, "Ti": 0.25, "P": 0.125, "Si": 0.125, "C": 0.125,
                "N": 0.125}]
    lp = LookUpData()
    x = APEAttributeGenerator(lp)
    x.generate_features(entries, verbose=True)
    # cProfile.run('x.generate_features(entries, verbose=True)')
