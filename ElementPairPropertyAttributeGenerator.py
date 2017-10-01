import types
import pandas as pd
import numpy as np
import sys

from LookUpData import LookUpData


class ElementPairPropertyAttributeGenerator:
    """
    Class to generate attributes based on the properties of constituent
    binary systems. Computes the minimum, maximum and range of all pairs in
    the material, and the fraction-weighted mean and variance of all pairs.
    Variance is defined as the mean absolute deviation from the mean over all
    pairs. If an entry has only one element, the value of NaN is used for all
    attributes.
    """
    def __init__(self, lp):
        self.lp = lp
        self.elemental_pair_properties = []
        self.pair_lookup_data = {}

    def load_pair_lookup_data(self):
        """
        Function to load the property values into self.lookup_data for the
        computation of features.
        :param lp: An instance of LookUpData required to load different
        property values into the dictionary.
        :return:
        """
        self.pair_lookup_data = self.lp.load_pair_properties(
            self.elemental_pair_properties)

    def add_elemental_pair_property(self, property):
        """
        Function to provide an elemental pair property to be used to compute
        features.
        :param property: Property to be included.
        :return:
        """
        if property not in self.elemental_pair_properties:
            self.elemental_pair_properties.append(property)

    def add_elemental_pair_properties(self, properties):
        """
        Function to provide a list of elemental pair properties to be used to
        compute features.
        :param properties: List of properties to be included.
        :return:
        """
        for prop in properties:
            self.add_elemental_pair_property(prop)

    def remove_elemental_pair_property(self, property):
        """
        Function to remove an elemental pair property from the list of elemental
        properties.
        :param property: Property to be removed.
        :return:
        """
        if property in self.elemental_pair_properties:
            self.elemental_pair_properties.remove(property)

    def remove_elemental_pair_properties(self, properties):
        """
        Function to remove a list of elemental pair properties from the list of
        elemental properties.
        :param properties: List of properties to be removed.
        :return:
        """
        for prop in properties:
            self.remove_elemental_pair_property(prop)

    def generate_features(self, entries, verbose=False):
        """
        Function to generate features of a binary material based on its
        constituent element properties.
        :param entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        :param verbose: Flag that is mainly used for debugging. Prints out a
        lot of information to the screen.
        :return features: Pandas data frame containing the names and values
        of the descriptors.
        """

        # Make sure that there is at least one elemental pair property provided.
        if not self.elemental_pair_properties:
            print "No elemental property is set. Add at least one property " \
                  "to compute meaningful descriptors."
            sys.exit(1)

        # If the dictionary containing the property values is empty,
        # load values into it.
        if not self.pair_lookup_data:
            self.load_pair_lookup_data()

        # Initialize lists of feature values and headers for pandas data frame.
        feat_values = []
        feat_headers = []

        # Raise exception if input argument is not of type list of dictionaries.
        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of dictionaries.")
        elif (entries and type(entries[0]) is not types.DictType):
            raise ValueError("Argument should be of type list of dictionaries.")

        # Insert header names here.
        n_statistics = 5
        for prop in self.elemental_pair_properties:
            feat_headers.append("binary_max_" + prop)
            feat_headers.append("binary_min_" + prop)
            feat_headers.append("binary_range_" + prop)
            feat_headers.append("binary_mean_" + prop)
            feat_headers.append("binary_variance_" + prop)

        for entry in entries:
            tmp_list = []
            elem_fractions = entry.values()
            elements = entry.keys()
            elem_ids = []
            for e in elements:
                elem_ids.append(self.lp.element_ids[e])
            if len(elem_fractions) == 1:
                for i in xrange(n_statistics):
                    tmp_list.append(np.nan)
                feat_values.append(tmp_list)
                continue

            pair_weights = []
            for i in xrange(len(elem_fractions)):
                for j in xrange(i):
                    pair_weights.append(elem_fractions[i]*elem_fractions[j])

            total_sum = sum(pair_weights)
            for i in xrange(len(pair_weights)):
                pair_weights[i] /= total_sum

            # Look up values for each pair property.
            for prop in self.elemental_pair_properties:
                tmp_prop = []

                for i in xrange(len(elem_fractions)):
                    e_i = elem_ids[i]
                    for j in xrange(i):
                        e_j = elem_ids[j]
                        idx_1 = max(e_i, e_j)
                        idx_2 = min(e_i, e_j)
                        tmp_prop.append(self.pair_lookup_data[prop][idx_1][
                                            idx_2])

                max_ = max(tmp_prop)
                min_= min(tmp_prop)
                range_ = max_ - min_
                mean_ = np.average(tmp_prop, weights=pair_weights)
                variance_ = np.average([abs(x - mean_) for x in tmp_prop],
                                      weights=pair_weights)
                tmp_list.append(max_)
                tmp_list.append(min_)
                tmp_list.append(range_)
                tmp_list.append(mean_)
                tmp_list.append(variance_)


            feat_values.append(tmp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    entry = [{"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}]
    y = LookUpData()
    x = ElementPairPropertyAttributeGenerator(y)
    x.add_elemental_pair_property("B2Volume")
    x.generate_features(entry, True)