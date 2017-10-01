import types
import math
import pandas as pd

from LookUpData import LookUpData
from OxidationStateGuesser import OxidationStateGuesser


class IonicityAttributeGenerator:
    """
    Class to generate the attributes based on the ionicity of a compound.
    Creates attributes based on whether it is possible to form a
    charge-neutral ionic compound, and two features based on a simple measure
    "bond ionicity" I(x,y)
    (see-http://www.wiley.com/WileyCDA/WileyTitle/productCd-EHEP002505.html
    W.D. Callister's text):

    I(x,y) = 1 - exp(-0.25* (chi_x - chi_y)^2)

    Maximum ionic character: Max I(x,y) between any two constituents.
    Mean ionic character: Sum x_i*x_j* I(i,j) where x_i is the fraction of
    element i and chi_x is the electronegativity of element x.
    """
    def __init__(self, lp):
        self.lp = lp

    def generate_features(self, entries, verbose=False):
        """
        Function to generate the three "ionicity" based features as described
        in the description of the class.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        :param verbose: Flag that is mainly used for debugging. Prints out a
        lot of information to the screen.
        :return features: Pandas data frame containing the names and values
        of the descriptors.
        """
        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Initialize lists of feature values and headers for pandas data frame.
        feat_values = []
        feat_headers = []

        # Insert header names here.
        feat_headers.append("CanFormIonic")
        feat_headers.append("MaxIonicChar")
        feat_headers.append("MeanIonicChar")

        # Instantiate and initialize oxidation state guesser with
        # electronegativity and oxidation state values.
        ox_guesser = OxidationStateGuesser()
        en = self.lp.load_property("Electronegativity")
        ox_guesser.set_electronegativity(en)
        ox_guesser.set_oxidationstates(self.lp.load_special_property(
            "OxidationStates"))

        for entry in entries:
            tmp_list = []
            # Can it form an ionic compound?
            tmp_list.append(0 if ox_guesser.get_possible_states(entry).size
                                 == 0 else 1)
            tmp_en = []

            # Compute mean ionic character.
            mean_ionic = 0.0
            for elem1 in entry:
                e1_id = self.lp.element_ids[elem1]
                e1_fraction = entry[elem1]
                tmp_en.append(en[e1_id])
                for elem2 in entry:
                    e2_id = self.lp.element_ids[elem2]
                    e2_fraction = entry[elem2]
                    m = 1 - math.exp(-0.25*(en[e1_id] - en[e2_id])**2)

                    mean_ionic += e1_fraction*e2_fraction*m

            # Compute max ionic character.
            max_ionic = 1 - math.exp(-0.25 * (max(tmp_en) - min(tmp_en)) ** 2)
            tmp_list.append(max_ionic)
            tmp_list.append(mean_ionic)
            feat_values.append(tmp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()

        return features

if __name__ == "__main__":
    entry = [{"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}]
    y = LookUpData()
    x = IonicityAttributeGenerator(y)
    x.generate_features(entry, True)