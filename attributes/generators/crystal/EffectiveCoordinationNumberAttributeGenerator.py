import types
import pandas as pd
import numpy as np
from data.materials.AtomicStructureEntry import AtomicStructureEntry
from vassal.analysis.VoronoiCellBasedAnalysis import VoronoiCellBasedAnalysis


class EffectiveCoordinationNumberAttributeGenerator:
    """
    Compute attributes based on the effective coordination number. The effective 
    coordination number can be thought of as a face-size-weighted coordination number.
    It is computed by the formula
    
    <center><i>N<sub>eff</sub></i> = 1 / sum[(<i>f<sub>i</sub></i>
    / <i>SA</i>)<sup>2</sup>]</center>
    
    where <i>f<sub>i</sub></i> is the area of face <i>i</i> and SA is the surface
    area of the entire cell. 
    
    <p>The effective coordination number has major benefit: stability against the
    additional of a very small face. Small perturbations in atomic positions
    can break symmetry in a crystal, and lead to the introduction of small faces.
    The conventional coordination number treats all faces equally, so the coordination
    number changes even when one of these small faces is added. 
    
    <p>One approach in the literature is to first apply a screen on small
    faces (e.g., remove any smaller than 1% of the total face area), which still
    runs into problems with discontinuity for larger displacements.
    
    <p>Our approach is differentiable with respect to the additional of a small face 
    (ask Logan if you want the math), and also captures another interesting effect
    small coordination numbers for Voronoi cells with a dispersity in face sizes.
    For example, BCC has 14 faces on its voronoi cell. 8 large faces, and 6 small ones.
    Our effective face size identifies a face size of closer to 8, the commonly-accepted
    value of the BCC coordination number, than 14 reported by the conventional measure.
    Additional, for systems with equal-sized faces (e.g., FCC), this measure
    agrees exactly with conventional reports.
    
    <usage><p><b>Usage</b>: *No options*</usage>
    
    @author Logan Ward
    """

    def mean_abs_dev(self, data):
        n = float(len(data))
        mean = sum(data) / n
        diff = [abs(x - mean) for x in data]
        return sum(diff) / n

    def generate_features(self, entries, verbose=False):
        # Raise exception if input argument is not of type list of
        # AtomicStructureEntry's.

        if (type(entries) is not types.ListType):
            raise ValueError("Argument should be of type list of "
                             "AtomicStructureEntry's")
        elif (entries and not isinstance(entries[0], AtomicStructureEntry)):
            raise ValueError("Argument should be of type list of "
                             "AtomicStructureEntry's")

        # Initialize lists of feature values and headers for pandas data frame.
        feat_headers = []
        feat_values = []

        feat_headers.append("mean_Coordination")
        feat_headers.append("var_Coordination")
        feat_headers.append("min_Coordination")
        feat_headers.append("max_Coordination")

        for entry in entries:
            temp_list = []
            output = entry.compute_voronoi_tessellation()
            N_eff = output.get_effective_coordination_numbers()

            mean = np.mean(N_eff)
            absdev = self.mean_abs_dev(data=N_eff)
            minimum = np.min(N_eff)
            maximum = np.max(N_eff)

            temp_list.append(mean)
            temp_list.append(absdev)
            temp_list.append(minimum)
            temp_list.append(maximum)

            feat_values.append(temp_list)

        features = pd.DataFrame(feat_values, columns=feat_headers)

        if verbose:
            print features.head()
            
        return features
