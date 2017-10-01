import pandas as pd
import sys
import types
import numpy as np

# TODO: Implement more rigorous tests
class ElementalPropertyAttributeGenerator:
    """
    Class to set up and generate descriptors based on elemental property
    statistics. Computes the mean, maximum, minimum, range, mode and mean
    absolute deviation of all elemental properties provided.
    """
    def __init__(self, lp, use_default_properties=True):
        """
        Class constructor.
        :param lp: Instance of the LookUpData class to be used by this class.
        :param use_default_properties: Flag to use default set of properties
        as defined below.
        """
        self.elemental_properties = []

        # Use default properties to reproduce Ward et al. descriptor values.
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
        self.lp = lp

        # Initialize dictionary that will contain all the property values.
        self.lookup_data = {}

    def load_lookup_data(self, lp):
        """
        Function to load the property values into self.lookup_data for the
        computation of features.
        :param lp: An instance of LookUpData required to load different
        property values into the dictionary.
        :return:
        """
        self.lookup_data = lp.load_properties(self.elemental_properties)

    def generate_features(self, entries, verbose=False):
        """
        Function to generate the elemental property based features. Computes
        6 statistics (mean, maximum, minimum, range, mode and mean absolute
        deviation) of all the elemental properties provided.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        :param verbose: Flag that is mainly used for debugging. Prints out a
        lot of information to the screen.
        :return features: Pandas data frame containing the names and values
        of the descriptors.
        """

        # Make sure that there is at least one elemental property provided.
        if not self.elemental_properties:
            print "No elemental property is set. Add at least one property " \
                  "to compute meaningful descriptors."
            sys.exit(1)

        # If the dictionary containing the property values is empty,
        # load values into it.
        if not self.lookup_data:
            self.load_lookup_data(self.lp)

        # Initialize lists of feature values and headers for pandas data frame.
        feat_values = []
        feat_headers = []

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Insert header names here.
        n_statistics = 6
        for prop in self.elemental_properties:
            feat_headers.append("mean_"+prop)
            feat_headers.append("maxdiff_" + prop)
            feat_headers.append("dev_" + prop)
            feat_headers.append("max_" + prop)
            feat_headers.append("min_" + prop)
            feat_headers.append("most_" + prop)

        missing_data = {}
        # Generate features for each entry.
        for entry in entries:
            elem_fractions = entry.values()
            max_f = max(elem_fractions)
            tmp_list = []
            # Look up values for each property.
            for prop in self.elemental_properties:
                tmp_prop = []
                for elem in entry:
                    elem_id = self.lp.element_ids[elem]
                    tmp_prop_value = self.lookup_data[prop][elem_id]
                    # If data is missing, make a note of it so that we can
                    # inform the user later.
                    if np.isnan(tmp_prop_value):
                        if not prop in missing_data:
                            missing_data[prop] = []
                        missing_data[prop].append(elem)
                    tmp_prop.append(tmp_prop_value)

                # If there is no missing data, compute statistics.
                if not np.isnan(tmp_prop).any():
                    mean_ = np.average(tmp_prop, weights=elem_fractions)
                    max_ = max(tmp_prop)
                    min_= min(tmp_prop)
                    max_diff_ = max_ - min_
                    avg_dev_ = np.average([abs(x - mean_) for x in tmp_prop],
                                          weights=elem_fractions)

                    indices = [i for i,f in enumerate(elem_fractions) if f >=
                               max_f]
                    most_ = sum(tmp_prop[i] for i in indices)/len(indices)

                    tmp_list.append(mean_)
                    tmp_list.append(max_diff_)
                    tmp_list.append(avg_dev_)
                    tmp_list.append(max_)
                    tmp_list.append(min_)
                    tmp_list.append(most_)
                else:
                    # Handle nan descriptors here from missing data.
                    for i in xrange(n_statistics):
                        tmp_list.append(np.nan)
            feat_values.append(tmp_list)

        # Issue warning to user about missing data here if it exists.
        if len(missing_data) > 0:
            sys.stderr.write("WARNING: There are " + str(len(missing_data)) +
                             " elemental properties with missing values: \n")
            for key in missing_data:
                sys.stderr.write("\t"+key+":")
                for elem in missing_data[key]:
                    sys.stderr.write(" "+elem)
                sys.stderr.write("\n")

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

    def add_elemental_property(self, property):
        """
        Function to provide an elemental property to be used to compute
        features.
        :param property: Property to be included.
        :return:
        """
        if property not in self.elemental_properties:
            self.elemental_properties.append(property)

    def add_elemental_properties(self, properties):
        """
        Function to provide a list of elemental properties to be used to
        compute
        features.
        :param properties: List of properties to be included.
        :return:
        """
        for prop in properties:
            self.add_elemental_property(prop)

    def remove_elemental_property(self, property):
        """
        Function to remove an elemental property from the list of elemental
        properties.
        :param property: Property to be removed.
        :return:
        """
        if property in self.elemental_properties:
            self.elemental_properties.remove(property)

    def remove_elemental_properties(self, properties):
        """
        Function to remove a list of elemental properties from the list of
        elemental properties.
        :param properties: List of properties to be removed.
        :return:
        """
        for prop in properties:
            self.remove_elemental_property(prop)

if __name__ == "__main__":
    entries = [{'Fe': 0.4, 'Cn': 0.6}, {'H': 0.6666666666666666,
                                       'O': 0.3333333333333333}, {'Na': 0.5,
                                                                   'Cl': 0.5}]
    x = ElementalPropertyAttributeGenerator(True)
    # x.add_elemental_property("AtomicWeight")
    x.generate_features(entries, verbose=True)