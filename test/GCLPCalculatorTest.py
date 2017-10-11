import unittest
from GCLPCalculator import GCLPCalculator
from LookUpData import LookUpData

class GCLPCalculatorTest(unittest.TestCase):
    def setUp(self):
        self.lp = LookUpData()
        self.calc = GCLPCalculator(self.lp)

    def tearDown(self):
        self.lp = None
        self.calc = None

    def test_initialization(self):
        n_elem = self.calc.get_num_phases()
        self.assertEqual(len(self.lp.element_names), n_elem, "Initial number "
                        "of phases should be equal to 112")

        # Add in NaCl
        NaCl = self.lp.get_sorted_and_normalized({"Na":1, "Cl":1})
        self.calc.add_phase(NaCl, -1)
