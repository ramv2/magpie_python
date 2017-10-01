import types
import pandas as pd
from LookUpData import LookUpData

class ElementalFractionAttributeGenerator:
    """
    Class to set the element fractions as the features of materials.
    """
    def __init__(self, lp):
        self.lp = lp

    def generate_features(self, entries, verbose=False):
        """
        Function to set element fractions as features of materials.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        :param verbose: Flag that is mainly used for debugging. Prints out a
        lot of information to the screen.
        :return features: Pandas data frame containing the names and values
        of the descriptors.
        """

        # Raise exception if input argument is not of type list of
        # dictionaries or if it is empty.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Initialize lists of feature values and headers for pandas data frame.
        feat_values = []
        feat_headers = []

        # Insert feature headers here.
        for elem in self.lp.element_ids:
            feat_headers.append("X_"+elem)

        # Insert feature values here.
        for entry in entries:
            tmp_list = []
            for elem in self.lp.element_ids:
                if elem in entry:
                    tmp_list.append(entry[elem])
                else:
                    tmp_list.append(0.0)
            feat_values.append(tmp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    entry = [{"Sc": 0.25, "Ti": 0.25, "P": 0.125, "Si": 0.125, "C": 0.125,
              "N": 0.125}]
    y = LookUpData()
    x = ElementalFractionAttributeGenerator(y)
    x.generate_features(entry, True)