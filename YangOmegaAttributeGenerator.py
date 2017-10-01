import types
import numpy as np
import math
import pandas as pd
from LookUpData import LookUpData

class YangOmegaAttributeGenerator:
    """
    Class to compute the attributes Omega and delta developed by Yang and Zhang.
    http://dx.doi.org/10.1016/j.matchemphys.2011.11.021
    These parameters are based on the liquid formation enthalpy and atomic
    sizes of elements respectively and were originally developed to predict
    whether a metal alloy will form a solid solution of bulk metallic glass.

    Omega; is derived from the melting temperature, ideal mixing entropy,
    and regular solution solution interaction parameter (Omega_{i,j}) predicted
    by the Miedema model for binary liquids. Specifically, it is computed
    using the relationship:
    Omega = T_m Delta S_mix / |Delta H_mix|
    where T_m is the composition-weighted average of the melting temperature,
    Delta S_mix is the ideal solution entropy, and Delta H_mix is he mixing
    enthalpy. The mixing enthalpy is computed using the Miedema mixing
    enthalpies tabulated by Takeuchi and Inoue
    (https://www.jstage.jst.go.jp/article/matertrans/46/12/46_12_2817/_article)
    where:
    Delta H_mix = Sum Omega_{i,j} c_i c_j and Omega_{i,j} = 4 *
    Delta H_mix.

    delta is related to the polydispersity of atomic sizes, and is computed
    using the relationship:

    delta = (Sum c_i (1-r_i/r_average)^2)^0.5

    where r_i is the atomic size. Here, we use the atomic radii compiled by
    Miracle et al.
    http://openurl.ingenta.com/content/xref?genre=article&issn=0950-6608
    &volume=55&issue=4&spage=218
    rather than those compiled by Kittel, as in the original work.
    """

    def __init__(self, lp):
        self.lp = lp

    def generate_features(self, entries, verbose=False):
        """
        Function to generate the features as mentioned in the class description.
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

        # Insert header names here.
        feat_headers.append("Yang_Omega")
        feat_headers.append("Yang_delta")

        # Load property values here.
        radii = self.lp.load_property("MiracleRadius")
        meltingT = self.lp.load_property("MeltingT")
        miedema = self.lp.load_pair_property("MiedemaLiquidDeltaHf")

        for entry in entries:
            tmp_list = []
            tmp_radii = []
            tmp_meltingT = []
            tmp_miedema = []
            e_ids = []
            element_fractions = entry.values()
            for elem in entry:
                elem_id = self.lp.element_ids[elem]
                e_ids.append(elem_id)
                tmp_radii.append(radii[elem_id])
                tmp_meltingT.append(meltingT[elem_id])

            # Compute the average melting point.
            averageTm = np.average(tmp_meltingT, weights=element_fractions)

            # Compute the ideal entropy.
            entropy = 0.0
            for f in element_fractions:
                entropy += f*math.log(f) if f > 0 else 0.0
            entropy *= 8.314/1000

            # Compute the enthalpy
            enthalpy = 0.0
            for i in xrange(len(e_ids)):
                for j in xrange(i+1,len(e_ids)):
                    enthalpy += miedema[max(e_ids[i], e_ids[j])][min(e_ids[i],
                                                                    e_ids[
                                                                        j])] \
                                * element_fractions[i] * element_fractions[j]
            enthalpy *= 4

            # Compute omega
            tmp_list.append(abs(averageTm * entropy / enthalpy))

            # Compute delta
            delta_squared = 0.0
            average_r = np.average(tmp_radii, weights=element_fractions)
            for i in xrange(len(e_ids)):
                delta_squared += element_fractions[i] * (1 - tmp_radii[
                    i] / average_r)**2

            tmp_list.append(math.sqrt(delta_squared))

            feat_values.append(tmp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)
        if verbose:
            print features.head()
        return features

if __name__ == "__main__":
    entry = [{"Sc": 0.25, "Ti": 0.25, "P": 0.125, "Si": 0.125, "C": 0.125,
              "N": 0.125}]
    y = LookUpData()
    x = YangOmegaAttributeGenerator(y)
    x.generate_features(entry, True)