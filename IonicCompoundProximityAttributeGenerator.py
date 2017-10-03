import types
from CompositionDistanceFilter import CompositionDistanceFilter
from IonicCompoundFinder import IonicCompoundFinder
import pandas as pd
from LookUpData import LookUpData

class IonicCompoundProximityAttributeGenerator:
    """
    Class to generate attributes based on the distance of a composition from
    a compositions that can form charge-neutral ionic compounds. This
    generator only computes a single feature: the L_1 distance between the
    composition of an entry and the nearest ionic compound (determined using
    IonicCompoundFinder). For a compound where it is not possible to form an
    ionic compound (e.g., only metallic elements), the entry is assigned
    arbitrarily large distance (equal to the number of elements in the alloy).
    The one adjustable parameter in this calculation is the maximum number of
    atoms per formula unit used when looking for ionic compounds. For binary
    compounds, the maximum conceivable number of elements in a formula unit
    is for a compound with a 9+ and a 5- species, which has 14 atoms in the
    formula unit. Consequently, we recommend using 14 or larger for this
    parameter.
    """
    def __init__(self, lp):
        """
        Initialize field variables here.
        :param lp: Instance of the LookUpData class required by this class.
        """
        self.lp = lp

        # Maximum number of atoms per formula unit.
        self.max_formula_unit = 14

    def set_max_formula_unit(self, size):
        """
        Function to define the maximum number of atoms per formula unit.
        :param size: Desired size.
        :return:
        """
        self.max_formula_unit = size

    def generate_features(self, entries, verbose=False):
        """
        Function to generate features as mentioned in the class description.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
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
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Insert header names here.
        feat_headers.append("IonicCompoundDistance_MaxSize"+str(
            self.max_formula_unit))

        # Get ionic compound finder.
        finder = IonicCompoundFinder(self.lp)
        finder.set_max_formula_unit_size(self.max_formula_unit)

        cdf = CompositionDistanceFilter(self.lp)
        for entry in entries:

            # Set the maximum distance to be equal to the number of elements.
            #  The maximum possible L_1 distance for an N-element system is N.
            finder.set_maximum_distance(len(entry))

            # If the material has only 1 element, set feature to 1.0.
            if len(entry) == 1:
                feat_values.append(1.0)
            else:
                # Get the list of all ionic compounds in the system.
                finder.set_nominal_composition(entry)
                ionic_compounds = finder.find_all_compounds()

                # Find the distance to the closest one. If no other
                # compounds, set distance to be the maximum possible.
                if not ionic_compounds:
                    feat_values.append(len(entry))
                else:
                    # print ionic_compounds[0]
                    # print cdf.compute_distance(entry, ionic_compounds[0], 1)
                    feat_values.append(cdf.compute_distance(entry,
                                                            ionic_compounds[
                                                                0], 1))

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    y = LookUpData()
    x = IonicCompoundProximityAttributeGenerator(y)
    entries = [{'Fe': 0.4, 'Cn': 0.6}, {'H': 0.6666666666666666,
                                        'O': 0.3333333333333333}, {'Na': 0.5,
                                                                   'Cl':
                                                                       0.5},
               {"Sc": 0.25, "Ti": 0.25, "P": 0.125, "Si": 0.125, "C": 0.125,
             "N": 0.125}]
    x.set_max_formula_unit(6)
    x.generate_features(entries, verbose=True)