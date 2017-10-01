import types
import pandas as pd
import numpy as np
from LookUpData import LookUpData

class MeredigAttributeGenerator:
    """
    Class to generate attributes as described by Meredig et al.
    http://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.094104
    This class is meant to be used in conjunction with
    ElementFractionAttributeGenerator and ValenceShellAttributeGenerator.

    To match the attributes from the Meredig et al. paper, use all three
    attribute generators.
    """
    def __init__(self, lp):
        self.lp = lp

    def generate_features(self, entries, verbose=False):
        """
        Function to generate attributes as described by Meredig et al.
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

        # Insert feature headers here.
        feat_headers.append("mean_AtomicWeight")
        feat_headers.append("mean_Column")
        feat_headers.append("mean_Row")
        feat_headers.append("maxdiff_AtomicNumber")
        feat_headers.append("mean_AtomicNumber")
        feat_headers.append("maxdiff_CovalentRadius")
        feat_headers.append("mean_CovalentRadius")
        feat_headers.append("maxdiff_Electronegativity")
        feat_headers.append("mean_Electronegativity")
        feat_headers.append("mean_NsValence")
        feat_headers.append("mean_NpValence")
        feat_headers.append("mean_NdValence")
        feat_headers.append("mean_NfValence")

        # Load all property tables.
        mass = self.lp.load_property("AtomicWeight")
        column = self.lp.load_property("Column")
        row = self.lp.load_property("Row")
        number = self.lp.load_property("Number")
        radius = self.lp.load_property("CovalentRadius")
        en = self.lp.load_property("Electronegativity")
        s = self.lp.load_property("NsValence")
        p = self.lp.load_property("NpValence")
        d = self.lp.load_property("NdValence")
        f = self.lp.load_property("NfValence")

        for entry in entries:
            tmp_list = []
            elem_fractions = entry.values()
            tmp_mass = []
            tmp_column = []
            tmp_row = []
            tmp_number = []
            tmp_radius = []
            tmp_en = []
            tmp_s = []
            tmp_p = []
            tmp_d = []
            tmp_f = []
            for elem in entry:
                elem_id = self.lp.element_ids[elem]
                tmp_mass.append(mass[elem_id])
                tmp_column.append(column[elem_id])
                tmp_row.append(row[elem_id])
                tmp_number.append(number[elem_id])
                tmp_radius.append(radius[elem_id])
                tmp_en.append(en[elem_id])
                tmp_s.append(s[elem_id])
                tmp_p.append(p[elem_id])
                tmp_d.append(d[elem_id])
                tmp_f.append(f[elem_id])

            tmp_list.append(np.average(tmp_mass, weights=elem_fractions))
            tmp_list.append(np.average(tmp_column, weights=elem_fractions))
            tmp_list.append(np.average(tmp_row, weights=elem_fractions))
            tmp_list.append(max(tmp_number) - min(tmp_number))
            tmp_list.append(np.average(tmp_number, weights=elem_fractions))
            tmp_list.append(max(tmp_radius) - min(tmp_radius))
            tmp_list.append(np.average(tmp_radius, weights=elem_fractions))
            tmp_list.append(max(tmp_en) - min(tmp_en))
            tmp_list.append(np.average(tmp_en, weights=elem_fractions))
            tmp_list.append(np.average(tmp_s, weights=elem_fractions))
            tmp_list.append(np.average(tmp_p, weights=elem_fractions))
            tmp_list.append(np.average(tmp_d, weights=elem_fractions))
            tmp_list.append(np.average(tmp_f, weights=elem_fractions))

            feat_values.append(tmp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    entry = [{"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}]
    y = LookUpData()
    x = MeredigAttributeGenerator(y)
    x.generate_features(entry, True)