import pandas as pd
import types

# TODO: Implement more rigorous tests
class StoichiometricAttributeGenerator:
    """
    Class to set up and generate descriptors based on the stoichiometry of a
    given material. Includes features that are only based on fractions of
    elements, but not what those elements actually are.
    """
    def __init__(self):
        """
        List of p norms to compute.
        """
        self.p_norms = []

    def clear_p_norms(self):
        """
        Clear out the list of p norms to be computed.
        :return:
        """
        del self.p_norms[:]

    def add_p_norm(self, norm):
        """
        Add a p norm to be computed.
        :param norm: desired norm.
        :return:
        """
        if (norm == 0):
            return
        elif (norm == 1):
            raise ValueError("L1 norm is always 1. Useless as attribute.")
        self.p_norms.append(norm)

    def add_p_norms(self, norms):
        """
        Add a list of p norms to be computed.
        :param norms: list of desired norms.
        :return:
        """
        for norm in norms:
            self.add_p_norm(norm)

    def generate_features(self, entries, verbose=False):
        """
        Function to generate the stiochiometric features. Computes the norms
        based on elemental fractions.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        :param verbose: Flag that is mainly used for debugging. Prints out a
        lot of information to the screen.
        :return features: Pandas data frame containing the names and values
        of the descriptors.
        """
        feat_values = []
        feat_headers = []

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Issue warning if no p norms are added.
        if (not self.p_norms):
            print "Warning: only L0 norm is computed."

        # Add in feature names.
        feat_headers.append("NComp")
        for p in self.p_norms:
            feat_headers.append("Comp_L"+str(p)+"Norm")

        # Compute features.
        for entry in entries:
            tmp_list = []
            n_comp = 0
            for f in entry.values():
                if (f > 0):
                    n_comp += 1
            # Number of components.
            tmp_list.append(n_comp)

            # Lp norms.
            for p in self.p_norms:
                tmp = 0.0
                for f in entry.values():
                    tmp += f**p
                tmp_list.append(tmp**(1.0/p))
            feat_values.append(tmp_list)

        # features as a pandas data frame.
        features = pd.DataFrame(feat_values, columns=feat_headers)
        if (verbose):
            print features.head()
        return features

if __name__ == "__main__":
    x = StoichiometricAttributeGenerator()
    x.add_p_norms([2,3,5,7,10])
    entries = [{'Fe': 0.4, 'O': 0.6}, {'H': 0.6666666666666666,
                                       'O': 0.3333333333333333}, {'Na': 0.5,
                                                                  'Cl': 0.5}]
    x.generate_features(entries, True)