import types
from scipy.optimize import linprog
import numpy as np
from itertools import izip


class GCLPCalculator:
    def __init__(self, lp):
        self.lp = lp
        self.phases = {}
        # Add for each element.
        for elem in self.lp.element_names:
            self.phases[elem] = 0.0

    def set_mu(self, elem, mu):
        if elem not in self.lp.element_names:
            raise ValueError("Not an element: "+elem)
        self.phases[elem] = mu

    def add_phases(self, entries, energies):
        for entry,energy in izip(entries, energies):
            # if has measurement
            self.add_phase(entry, energy)

    def add_phase(self, entry, energy):
        if entry not in self.phases:
            self.phases[entry] = energy
        elif self.phases[entry] > energy:
            self.phases[entry] = energy

    def get_num_phases(self):
        return len(self.phases)

    def run_GCLP(self, composition):

        # List of dictionaries.
        components = []

        # List.
        energies = []

        # Get the current possible phases (i.e., those that contain
        # exclusively the elements in the current compound.
        for entry in self.phases:
            # Check whether this entry is in the target phase diagram.
            if set(entry.keys()) <= set(composition.keys()):
                components.append(entry)
                energies.append(self.phases[entry])

        # Set up constraints.
        # Type #1: Mass conservation.
        l_c = len(components)
        a_eq = np.ones((len(composition)+1, l_c))
        b_eq = np.ones(len(composition)+1)
        for i,elem in enumerate(composition):
            b_eq[i] = composition[elem]
            for j in xrange(l_c):
                a_eq[i][j] = components[j][elem]

        # Type #2: Normalization.
        # Taken care of when we initialized a_eq and b_eq to ones.
        c = np.array(energies)

        res = linprog(c=c, A_eq=a_eq, b_eq=b_eq)
        equilibrium =[[components[i], res.x[i]] for i in range(l_c) if res.x[
            i] > 1E-6]
        return res.fun, equilibrium