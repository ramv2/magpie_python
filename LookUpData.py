import numpy as np
import sys

class LookUpData:
    """
    Class to look up properties of elements stored in files.
    """
    # Location of lookup tables.
    lookup_location = "lookup-data/"
    pair_location = "lookup-data/pair/"

    # Element indices of the periodic table.
    element_ids = {"H": 0, "He": 1, "Li": 2, "Be": 3, "B": 4, "C": 5,
                   "N": 6, "O": 7, "F": 8, "Ne": 9, "Na": 10, "Mg": 11,
                   "Al": 12, "Si": 13, "P": 14, "S": 15, "Cl": 16,
                   "Ar": 17, "K": 18, "Ca": 19, "Sc": 20, "Ti": 21,
                   "V": 22, "Cr": 23, "Mn": 24, "Fe": 25, "Co": 26,
                   "Ni": 27, "Cu": 28, "Zn": 29, "Ga": 30, "Ge": 31,
                   "As": 32, "Se": 33, "Br": 34, "Kr": 35, "Rb": 36,
                   "Sr": 37, "Y": 38, "Zr": 39, "Nb": 40, "Mo": 41,
                   "Tc": 42, "Ru": 43, "Rh": 44, "Pd": 45, "Ag": 46,
                   "Cd": 47, "In": 48, "Sn": 49, "Sb": 50, "Te": 51,
                   "I": 52, "Xe": 53, "Cs": 54, "Ba": 55, "La": 56,
                   "Ce": 57, "Pr": 58, "Nd": 59, "Pm": 60, "Sm": 61,
                   "Eu": 62, "Gd": 63, "Tb": 64, "Dy": 65, "Ho": 66,
                   "Er": 67, "Tm": 68, "Yb": 69, "Lu": 70, "Hf": 71,
                   "Ta": 72, "W": 73, "Re": 74, "Os": 75, "Ir": 76,
                   "Pt": 77, "Au": 78, "Hg": 79, "Tl": 80, "Pb": 81,
                   "Bi": 82, "Po": 83, "At": 84, "Rn": 85, "Fr": 86,
                   "Ra": 87, "Ac": 88, "Th": 89, "Pa": 90, "U": 91,
                   "Np": 92, "Pu": 93, "Am": 94, "Cm": 95, "Bk": 96,
                   "Cf": 97, "Es": 98, "Fm": 99, "Md": 100, "No": 101,
                   "Lr": 102, "Rf": 103, "Db": 104, "Sg": 105,
                   "Bh": 106, "Hs": 107, "Mt": 108, "Ds": 109,
                   "Rg": 110, "Cn": 111}

    def load_property(self, property, lookup_dir=lookup_location):
        """
        Function to load a specific property from the directory containing
        all the lookup tables.
        :param property: Property whose values need to be loaded.
        :param lookup_dir: Directory containing all the property value files.
        :return: values: A numpy array containing the property values for all
        the elements.
        """
        # IonizationEnergies and OxidationStates are 2-D arrays. So treat
        # them differently.
        if property == "IonizationEnergies" or property == "OxidationStates":
            print "Use special functions readIonizationEnergies or " \
                  "readOxidationStates to read these properties"
            sys.exit(1)

        # Initialize the numpy array.
        values = np.zeros(len(self.element_ids), dtype=np.float)
        values.fill(np.nan)

        # Property file name.
        file = lookup_dir+property+".table"
        try:
            prop_file = open(file, 'r')
        except IOError:
            print "File {} doesn't exist!!! Please make sure you " \
                  "specify the correct file name".format(file)
            sys.exit(1)
        else:
            for i in xrange(values.size):
                line = prop_file.readline().strip()
                if line[0].isdigit():
                    values[i] = float(line)
            prop_file.close()
        return values

    def load_pair_property(self, property, data_dir=pair_location):
        """
        Function to load property of a binary system.
        :param property: Property whose values need to be loaded.
        :param data_dir: Directory containing all the property value files.
        :return: values: A 2-D numpy array containing the property values for
        all the elements.
        """

        # Initialize the 2-D numpy array.
        values = np.zeros(len(self.element_ids), dtype=object)
        for i in xrange(len(self.element_ids)):
            values[i] = np.zeros(i, dtype=float)

        # Property file name.
        file = data_dir + property + ".table"

        try:
            prop_file = open(file, 'r')
        except IOError:
            print "File {} doesn't exist!!! Please make sure you " \
                  "specify the correct file name".format(file)
            sys.exit(1)
        else:
            for line in prop_file.readlines():
                words = line.strip().split()
                if len(words) < 3:
                    continue
                else:
                    if words[0] not in self.element_ids or words[1] not in \
                            self.element_ids:
                        continue
                    elemA = self.element_ids[words[0]]
                    elemB = self.element_ids[words[1]]
                    if words[2].endswith("\n"):
                        print line
                        sys.exit(1)
                    values[max(elemA,elemB)][min(elemA,elemB)] = float(words[2])
            prop_file.close()
        return values

    def load_pair_properties(self, properties, data_dir=pair_location):
        """
        Function to load multiple pair property values from the directory
        containing all the lookup tables.
        :param properties: A list of pair properties whose values need to be
        loaded.
        :param data_dir: Directory containing all the pair property value files.
        :return: values: A dictionary containing <Property,Info> as <key,
        value> where Info is a numpy array containing the pair property values.
        """

        # Initialize the dictionary.
        values = {}

        for prop in properties:
            values[prop] = self.load_pair_property(prop, data_dir)
        return values

    def load_properties(self, properties, lookup_dir=lookup_location):
        """
        Function to load multiple property values from the directory
        containing all the lookup tables.
        :param properties: A list of properties whose values need to be loaded.
        :param lookup_dir: Directory containing all the property value files.
        :return: values: A dictionary containing <Property,Info> as <key,
        value> where Info is a numpy array containing the property values.
        """

        # Initialize the dictionary.
        values = {}

        for prop in properties:
            values[prop] = self.load_property(prop, lookup_dir)
        return values

    def load_special_property(self, property, lookup_dir=lookup_location):
        """
        Function to load the special property files related to
        IonizationEnergies and OxidationStates.
        :param property: Property whose values need to be loaded.
        :param lookup_dir: Directory containing all the property value files.
        :return: values: A 2-D numpy array containing the property values
        whose dtype is object.
        """

        # Make sure property is either IonizationEnergies or OxidationStates.
        if not (property == "IonizationEnergies" or property == \
                "OxidationStates"):
            print "Use load_property function. This function is used " \
                  "exclusively for IonizationEnergies and OxidationStates"
            sys.exit(1)

        # Property file name.
        file = lookup_dir +property+".table"

        # Initialize the list.
        tmp_values = []

        try:
            prop_file = open(file, 'r')
        except IOError:
            print "File {} doesn't exist!!! Please make sure you " \
                  "specify the correct file name".format(file)
            sys.exit(1)
        else:
            for line in prop_file.readlines():
                words = line.strip().split()
                tmp_list = []
                for word in words:
                    tmp_list.append(float(word))
                tmp_values.append(tmp_list)
            prop_file.close()

        values = np.zeros(len(tmp_values), dtype=object)
        for i in xrange(len(tmp_values)):
            values[i] = np.asarray(tmp_values[i], dtype=float)
        return values

    def set_lookup_location(self, location):
        """
        Set a different lookup location.
        :param location: Path to new directory containing all the property
        value files.
        :return:
        """
        lookup_location = location

if __name__ == "__main__":
    x = LookUpData()
    # a = x.load_property("Electronegativity")
    # print a
    # b = x.load_special_property("IonizationEnergies")
    # for i in b:
    #     print i
    # c = x.load_pair_property("B2Volume")
    # print c[5]
    # d = ["Number", "MendeleevNumber"]
    # e = x.load_properties(d)
    # print e["MendeleevNumber"]