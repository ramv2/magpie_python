from EqualSumCombinations import EqualSumCombinations
import numpy as np
from itertools import combinations as comb
from LookUpData import LookUpData

class PhaseDiagramCompositionEntryGenerator:
    """
    Class to generate composition entries at many points in many phase
    diagrams. Has two different ways of selecting compositions within phase
    diagrams.

    1. Even Spacing: Compositions are selected to be evenly-spaced throughout
    the phase diagram (e.g. A0.2B0.2C0.6, A0.4B0.2C0.4 etc.). This method is
    most appropriate for alloy systems.

    2. Simple Fractions: Compositions with the smallest denominator are
    selected (e.g. ABC, A2C, B2C, etc.). This method is most appropriate for
    phase diagrams that represent ordered crystalline compounds.
    """
    def __init__(self, lp):
        """
        Initialize field variables here.
        :param lp: Instance of the LookUpData class required fo getting the
        element names.
        """
        self.lp = lp

        # List of elements to use (id is Z-1).
        self.e_names = []

        # Minimum number of constituents.
        self.min_order = 1

        # Maximum number of constituents.
        self.max_order = 1

        # Whether to use even spacing or small integers.
        self.even_spacing = True

        # Either number of stops in each direction or max denominator.
        self.size = 3

    def set_elements_by_index(self, indices):
        """
        Function to define the list of elements to be included in the phase
        diagrams.
        :param indices: List of elements by index (Z-1).
        :return:
        """
        for i in indices:
            if i < len(self.lp.element_names):
                self.e_names.append(self.lp.element_names[indices])
            else:
                raise ValueError("Index out of range: "+str(i))

    def set_elements_by_name(self, names):
        """
        Function to define the list of elements to be included in the phase
        diagrams.
        :param names: List of element names.
        :return:
        """
        self.e_names = names

    def set_order(self, min_, max_):
        """
        Function to define the order of generated phase diagrams.
        :param min_: Minimum number of constituents.
        :param max_: Maximum number of constituents.
        :return:
        """
        if min_ < 1 or max_ < 1:
            raise ValueError("Orders must be greater than 1.")
        self.min_order = min_
        self.max_order = max_

    def set_even_spacing(self, es):
        """
        Function to define whether to use evenly-spaced compositions (or
        low-denominator).
        :param es: Boolean indicating the same.
        :return:
        """
        self.even_spacing = es

    def set_size(self, size):
        """
        Function to define the number of points per binary diagram (in using
        even spacing) or the maximum denominator (for low-denominator).
        :param size: Desired size parameter.
        :return:
        """
        if size < 2:
            raise ValueError("Size must be greater than 1.")
        self.size = size

    def generate_alloy_compositions(self):
        """
        Function to generate evenly-spaced compositions. Generates
        compositions for all diagrams up to the user-specified Minimum order.

        For example: If the user wants ternary diagrams with a minimum
        spacing of 0.25 this code will generate the following map:
        1 -> ([1.0])
        2 -> ([0.25, 0.75], [0.5, 0.5], [0.75, 0.25])
        3 -> ([0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5])
        :return: A dictionary containing <Order, Possible compositions> as
        <key,value> pairs. Here, Order is the number of elements and Possible
        compositions is a list of numpy arrays containing the fractions of
        elements.
        """

        output = {}

        # Add in diagrams of greater order.
        for order in xrange(self.min_order, self.max_order+1):
            if order == 1:
                tmp_list = []
                tmp_list.append(np.array([1.0]))
                output[order] = tmp_list
                continue

            tmp_list = []
            es = EqualSumCombinations(self.size-1, order)
            for compI in es.get_combinations():
                if 0 in compI:
                    # Don't add compositions from a lower-order diagram.
                    continue
                comp = np.zeros(order)
                for i in xrange(order):
                    comp[i] = compI[i]/ float(self.size - 1.0)
                tmp_list.append(comp)
            output[order] = tmp_list
        return output

    def generate_crystal_compositions(self):
        """
        Function to generate compositions with small denominators. Generates
        compositions for all diagrams up to the user-specified Minimum order.

        For example: If the user wants ternary diagrams with a maximum
        denominator of 3 this code will generate the following map:
        1 -> ([1])
        2 -> ([1/3, 2/3], [1/2, 1/2], [2/3, 1/3])
        3 -> ([1/3, 1/3, 1/3])
        :return: A dictionary containing <Order, Possible compositions> as
        <key,value> pairs. Here, Order is the number of elements and Possible
        compositions is a list of numpy arrays containing the fractions of
        elements.
        """

        output = {}

        # Add in diagrams of greater order.
        for order in xrange(self.min_order, self.max_order+1):
            if order == 1:
                tmp_list = []
                tmp_list.append(np.array([1.0]))
                output[order] = tmp_list
                continue

            tmp_list = []
            reduced_examples = []
            for d in xrange(order, self.size+1):
                es = EqualSumCombinations(d, order)
                for compI in es.get_combinations():
                    if 0 in compI:
                        # Don't add compositions from a lower-order diagram.
                        continue
                    comp = np.zeros(order)
                    red_comp = np.zeros(order)
                    for i in xrange(order):
                        comp[i] = float(compI[i])
                        red_comp[i] = comp[i] / d

                    # Check if this composition is already represented.
                    was_found = False
                    for ex in reduced_examples:
                        if np.array_equal(red_comp, ex):
                            was_found = True
                            break

                    if not was_found:
                        tmp_list.append(comp)
                        reduced_examples.append(red_comp)
            output[order] = tmp_list
        return output

    def generate_entries(self):
        """
        Function to generate the list of entries corresponding to the list of
        compositions, element names specified by the user and the mapping of
        number of elements to compositions.
        :return: entries: A list of dictionaries containing <Element name,
        fraction> as <key,value> pairs.
        """
        entries = []

        # Get the correct composition mapping.
        compositions = self.generate_alloy_compositions() if \
            self.even_spacing else self.generate_crystal_compositions()

        for order,list_of_fractions in compositions.iteritems():

            # Generate all possible combinations of elements of a given order.
            compounds = [list(i) for i in comb(self.e_names, order)]
            # print order, len(compounds)*len(list_of_fractions)
            for fractions in list_of_fractions:
                sum_ = sum(fractions)
                for compound in compounds:
                    entry = {}
                    for i in xrange(len(compound)):
                        entry[compound[i]] = float(fractions[i])/sum_
                    entries.append(entry)
        return entries

if __name__ == "__main__":
    y = LookUpData()
    x = PhaseDiagramCompositionEntryGenerator(y)
    x.set_elements_by_name(["C", "N", "P", "Si", "Sc", "Ti"])
    x.set_even_spacing(True)
    x.set_order(1, 2)
    x.set_size(5)
    # x.generate_entries()
    for i in x.generate_entries():
        print i