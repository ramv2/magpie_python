import numpy as np
import sys
import itertools
from LookUpData import LookUpData

class OxidationStateGuesser:
    """
    Class to predict the likely oxidation states of a material, given its
    input composition.
    """
    def __init__(self):
        """
        Initialize variables here.
        """
        self.electronegativity = np.zeros(0)
        self.oxidationstates = np.zeros(0, dtype=object)


    def set_electronegativity(self, values):
        """
        Function to set the electronegativity values.
        :param values: Numpy array containing electronegativity values for
        all th elements.
        :return:
        """
        self.electronegativity = values

    def set_oxidationstates(self, values):
        """
        Function to set the oxidation states values.
        :param values: 2-D numpy array containing oxidation states values for
        all th elements.
        :return:
        """
        self.oxidationstates = values

    def get_possible_states(self, entry):
        """
        Function to compute all the possible oxidation states of a material,
        given its input composition. The function works by finding all
        combinations of non-zero oxidation states for each element, computing
        which are the most reasonable, and finding which of those have minimum
        value of
        sum_{i,j} (chi_i - chi_j)*(c_i - c_j) for i < j
        where chi_i is the electronegativity and c_i is the oxidation. This
        biases the selection towards the more electronegative elements being
        more negatively charged.
        :param entry: Dictionary containing the element names and fractions
        as keys and values respectively.
        :return: output: A numpy array containing the list of possible
        oxidation states arranged in the order mentioned above.
        """

        # Make sure entry is not empty.
        if not entry:
            print "Input argument cannot be empty. Please pass a valid " \
                  "argument."
            sys.exit(1)

        # Make sure electronegativity and oxidation states are not empty.
        if not self.electronegativity.size or not self.oxidationstates.size:
            print "Electronegativity or OxidationStates values are not " \
                  "initialized. Set them and try again."
            sys.exit(1)

        # Make sure element fractions add up to 1.0.
        if sum(entry.values()) != 1.0:
            print "Entry should be a dictionary containing element names and " \
                  "fractions as keys and values respectively. Also, " \
                  "fractions should add up to 1.0."
            sys.exit(1)

        # Initialize list of possible states.
        possible_states = []

        # Get element ids.
        elem_ids = []
        for elem in entry:
            elem_ids.append(LookUpData.element_ids[elem])

        # List of all states.
        states = []
        for id in elem_ids:
            states.append(self.oxidationstates[id])

        for state in itertools.product(*states):
            charge = np.dot(np.asarray(state), np.asarray(entry.values()))
            # If charge is balanced, add state to the list of possible states.
            if abs(charge) < 1E-6:
                possible_states.append(list(state))

        if len(possible_states) < 2:
            return np.asarray(possible_states)

        # Compute the summation mentioned in the function description.
        rankVal = np.zeros(len(possible_states))
        for s in xrange(len(possible_states)):
            state = possible_states[s]
            tmp_val = 0.0
            for i in xrange(len(state)):
                for j in xrange(i+1,len(state)):
                    tmp_val += (self.electronegativity[elem_ids[i]] -
                                self.electronegativity[elem_ids[j]]) * (
                        state[i] - state[j])
            rankVal[s] = tmp_val

        # Order them based on electronegatvity rank.
        ranks = np.argsort(rankVal)
        output = []
        for i in xrange(len(possible_states)):
            output.append(possible_states[ranks[i]])

        return np.asarray(output)

if __name__ == "__main__":
    entry = {"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}
    x = OxidationStateGuesser()
    y = x.get_possible_states(entry)