import pandas as pd
import os
import sys
import types

class ElementalPropertyAttributeGenerator:
    def __init__(self, use_default_properties=True):
        self.elemental_properties = []
        if (use_default_properties):
            self.elemental_properties = ["Number", "MendeleevNumber",
                                         "AtomicWeight", "MeltingT",
                                         "Column", "Row", "CovalentRadius",
                                         "Electronegativity", "NsValence",
                                         "NpValence", "NdValence",
                                         "NfValence", "NValance",
                                         "NsUnfilled", "NpUnfilled",
                                         "NdUnfilled", "NfUnfilled",
                                         "NUnfilled", "GSvolume_pa",
                                         "GSbandgap", "GSmagmom",
                                         "SpaceGroupNumber"]
        self.lookup_location = "lookup-data/"
        self.lookup_data = {}
        self.element_names = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F",
                              "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                              "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
                              "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
                              "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb",
                              "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
                              "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                              "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
                              "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
                              "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
                              "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
                              "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
                              "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
                              "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"]
        self.load_lookup_data()

    def load_lookup_data(self):
        for file in os.listdir(self.lookup_location):
            if (file.endswith(".table")):
                prop = file.split(".table")[0]
                values = []
                try:
                    prop_file = open(file, 'r')
                except IOError:
                    print "File {} doesn't exist!!! Please make sure you " \
                          "specify the correct file name".format(
                        file)
                    sys.exit(1)
                else:
                    for line in prop_file.readlines():
                        if (not line.strip().startswith("Miss")):
                            values.append(float(line.strip()))
                        else:
                            values.append(float('nan'))
                    self.lookup_data[prop] = values

    def get_features(self, entries):
        feat_values = []
        feat_headers = []

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        n_statistics = 6
        for prop in self.elemental_properties:
            feat_headers.append("mean_"+prop)
            feat_headers.append("maxdiff_" + prop)
            feat_headers.append("dev_" + prop)
            feat_headers.append("max_" + prop)
            feat_headers.append("min_" + prop)
            feat_headers.append("most_" + prop)

        # Handle missing data here
        for entry in entries:
            tmp_list = []
            for prop in self.elemental_properties:

        features = pd.DataFrame(feat_values, columns=feat_headers)
        return features