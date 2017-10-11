import types
import numpy as np
from LookUpData import LookUpData

class CompositionDistanceFilter:
    """
    Class to filter compositions based on distance from a target composition.
    Filters any composition where the maximum change in any element is less
    than a certain value.
    """
    def __init__(self, lp):
        """
        Initialize field variables here.
        :param lp: Instance of the LookUpData class required by this class.
        """

        # Target composition.
        self.target_composition = {}

        # Threshold distance.
        self.threshold = 0.0

        self.lp = lp

    def set_target_composition(self, entry):
        """
        Function to define the target composition.
        :param entry: Desired target composition with element names and
        fractions as keys and values respectively.
        :return:
        """
        self.target_composition = entry

    def set_distance_threshold(self, distance):
        """
        Function to define the threshold composition distance. Here distance
        is defined as the maximum change in the fraction of any element.
        :param distance: Target threshold in %
        :return:
        """
        self.threshold = distance / 100.0

    @classmethod
    def compute_distance(self, entry_1, entry_2, p):
        """
        Function to compute the distance between two entries. Distance is
        defined as the L_p norm of the distance between element fractions.
        :param entry_1: Dictionary 1 with element names and fractions as keys
        and values respectively.
        :param entry_2: Dictionary 1 with element names and fractions as keys
        and values respectively.
        :param p: Desired norm.
        :return: dist: Distance between two entries.
        """

        # Get the list of common elements with fractions greater than zero to
        # consider.
        elements = []
        for e in entry_1:
            if entry_1[e] > 0.0:
                elements.append(e)
        for e in entry_2:
            if e not in elements and entry_2[e] > 0.0:
                elements.append(e)

        # Compute differences
        dist = 0.0
        for e in elements:
            f1 = entry_1[e] if e in entry_1 else 0
            f2 = entry_2[e] if e in entry_2 else 0
            diff = f1 - f2
            if p == 0:
                if diff != 0:
                    dist += 1
            elif p == -1:
                dist = max(dist, diff)
            else:
                dist += abs(diff)**p

        # Compute distance
        if p == 0 or p == 1:
            return dist

        return dist**(1.0 / p)

    def label(self, entries):
        """
        Function to compute labels of composition entries indicating whether
        or not they are within the threshold of the target composition.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        :return: label: A numpy array containing True if the entry is within
        bounds of the target composition and False otherwise.
        """

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        label = np.array([False]*len(entries))
        for i,entry in enumerate(entries):
            within_bounds = True
            for e in self.target_composition:
                v1 = entry[e] if e in entry else 0
                if abs(v1 - self.target_composition[e]) > self.threshold:
                    within_bounds = False
                    break
            label[i] = within_bounds

        return label
if __name__ == "__main__":
    y = LookUpData()
    x = CompositionDistanceFilter(y)
    e_1 = {"Sc":0.25, "Ti":0.25, "P":0.125, "Si":0.125, "C":0.125, "N":0.125}
    e_2 = {'C': 0.25, 'P':0.0, 'Sc': 0.75}

    print x.compute_distance(e_1, e_2, 2)