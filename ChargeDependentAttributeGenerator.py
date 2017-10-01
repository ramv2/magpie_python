import types
import numpy as np
import pandas as pd
import sys
from LookUpData import LookUpData
from OxidationStateGuesser import OxidationStateGuesser

class ChargeDependentAttributeGenerator:
    """
    Class to generate attributes derived from the oxidation states of
    elements in a material. Based on work by Deml et al.
    http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.085142
    Deml et al. PRB. 93 (2016), 085142
    These features are based on the formal charges of materials determined
    using the OxidationStateGuesser. Currently implemented features:
    Statistics of formal charges (min, max, range, mean, variance)
    Cumulative ionization energies/ electron affinities
    Difference in electronegativities between cation and anion.

    For materials that the algorithm fails to find charge states, NaN is set
    for all features.
    """
    def __init__(self, lp):
        self.lp = lp

    def generate_features(self, entries, verbose=False):
        """
        Function to generate the charge dependent features as mentioned in
        the class description.
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

        # Insert feature headers here.
        n_features = 8
        feat_headers.append("min_Charge")
        feat_headers.append("max_Charge")
        feat_headers.append("maxdiff_Charge")
        feat_headers.append("mean_Charge")
        feat_headers.append("var_Charge")
        feat_headers.append("CumulativeIonizationEnergy")
        feat_headers.append("CumulativeElectronAffinity")
        feat_headers.append("AnionCationElectronegativtyDiff")

        # Load properties here.
        en = self.lp.load_property("Electronegativity")
        ea = self.lp.load_property("ElectronAffinity")
        ie = self.lp.load_special_property("IonizationEnergies")

        # Instantiate and initialize oxidation state guesser.
        ox_guesser = OxidationStateGuesser()
        ox_guesser.set_electronegativity(en)
        ox_guesser.set_oxidationstates(self.lp.load_special_property(
            "OxidationStates"))

        missing_data = {}
        for entry in entries:
            tmp_list = []

            # Get possible states with charges.
            possible_states = ox_guesser.get_possible_states(entry)
            elem_fractions = entry.values()

            # If there are no possible states, set all features to NaN.
            if len(possible_states) == 0:
                for i in xrange(n_features):
                    tmp_list.append(np.nan)

                feat_values.append(tmp_list)
                continue

            # Check that we have data for all ionization energies.
            any_missing = False
            tmp_charges = possible_states[0]
            for i,elem in enumerate(entry):
                elem_id = self.lp.element_ids[elem]
                if len(ie[elem_id]) < tmp_charges[i]:
                    if elem_id not in missing_data:
                        missing_data[elem_id] = []
                    missing_data[elem_id].append(possible_states[0][i])
                    any_missing = True
                    break

            # Compute statistics related to charges.
            min_ = min(tmp_charges)
            max_ = max(tmp_charges)
            max_diff_ = max_ -  min_
            mean_ = np.average([abs(x) for x in tmp_charges],
                               weights=elem_fractions)
            var_ = np.average([abs(abs(x) - mean_) for x in tmp_charges],
                              weights=elem_fractions)

            tmp_list.append(min_)
            tmp_list.append(max_)
            tmp_list.append(max_diff_)
            tmp_list.append(mean_)
            tmp_list.append(var_)

            if any_missing:
                tmp_list.append(np.nan)
                tmp_list.append(np.nan)
                tmp_list.append(np.nan)
                feat_values.append(tmp_list)
                continue

            # Compute features related to ionization/affinity.
            cation_fraction = anion_fraction = cation_ie_sum = anion_ea_sum = \
            mean_cation_en = mean_anion_en = 0.0
            for e,elem in enumerate(entry):
                elem_id = self.lp.element_ids[elem]
                if tmp_charges[e] < 0:
                    anion_fraction += elem_fractions[e]
                    mean_anion_en +=  en[elem_id] * elem_fractions[e]
                    anion_ea_sum -= tmp_charges[e] * ea[elem_id]* \
                                    elem_fractions[e]
                else:
                    cation_fraction += elem_fractions[e]
                    mean_cation_en += en[elem_id] * elem_fractions[e]
                    cation_ie_sum += sum(ie[elem_id][c] * elem_fractions[e] for
                                        c in xrange(int(tmp_charges[e])))

            mean_anion_en /= anion_fraction
            mean_cation_en /= cation_fraction
            anion_ea_sum /= anion_fraction
            cation_ie_sum /= cation_fraction

            tmp_list.append(cation_ie_sum)
            tmp_list.append(anion_ea_sum)
            tmp_list.append(mean_anion_en - mean_cation_en)

            feat_values.append(tmp_list)

        # Issue warning to user about missing data here if it exists.
        if len(missing_data) > 0:
            sys.stderr.write("WARNING: Missing ionization energy data for\n")
            for elem in missing_data:
                sys.stderr.write("\t" + elem + ":")
                for state in missing_data[elem]:
                    sys.stderr.write(" +" + state)
                sys.stderr.write("\n")

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    entry = [{"Sc":0.25,"Ti":0.25,"P":0.125,"Si":0.125,"C":0.125,"N":0.125}]
    y = LookUpData()
    x = ChargeDependentAttributeGenerator(y)
    x.generate_features(entry, True)