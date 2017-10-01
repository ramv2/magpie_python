import types
import numpy as np
import pandas as pd
from LookUpData import LookUpData

class ValenceShellAttributeGenerator:
    """
    Class that generates attributes based on fraction of electrons in valence
    shell of constituent elements. Creates 4 feature: [Composition-weighted
    mean # of electrons in the {s,p,d,f} shells]/[Mean # of Valence Electrons]

    Originally presented by:
    http://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.094104
    Meredig et al. Physical Review B (2015)
    """
    def __init__(self, lp):
        self.lp = lp

    def generate_features(self, entries, verbose=False):
        """
        Function that generates the attributes mentioned in the class
        description above.
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

        shell = ['s','p','d','f']
        n_valence = np.zeros((len(shell), len(self.lp.element_ids)))

        # Read in the data from tables and insert feature headers here.
        for i in xrange(len(shell)):
            s = shell[i]
            feat_headers.append("frac_"+s+"Valence")
            n_valence[i] = self.lp.load_property("N"+s+"Valence")

        for entry in entries:
            sum_e = 0.0
            total_e = []
            for i in xrange(len(shell)):
                tmp_valence = []
                for elem in entry:
                    elem_id = self.lp.element_ids[elem]
                    tmp_valence.append(n_valence[i][elem_id])

                # Fraction weighted average # of electrons in this shell.
                x = np.average(tmp_valence, weights=entry.values())
                sum_e += x
                total_e.append(x)

            for i in xrange(len(total_e)):
                total_e[i] /= sum_e

            feat_values.append(total_e)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    entry = [{"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}]
    y = LookUpData()
    x = ValenceShellAttributeGenerator(y)
    x.generate_features(entry, True)