"""
Unit tests for quantum analysis. Takes in pre-made data with known results,
computes the results from the data and checks everything is correct.
"""
import unittest
import pickle
import numpy as np
import sys

sys.path.append("..")  # noqa
import pyEPR as epr

# Files location
save_file = "./data.npz"
correct_results = "./correct_results.pkl"


class TestQuantumAnalysis(unittest.TestCase):
    def setUp(self):
        self.epra = epr.QuantumAnalysis(save_file)
        with open(correct_results, "rb") as file:
            self.correct_res = pickle.load(file)

    def test_analyze_all_variations(self):
        """
        Check that the calculated results matches the known correct ones
        """
        results = self.epra.analyze_all_variations(
            cos_trunc=8, fock_trunc=15, print_result=False
        )[
            "0"
        ]  # Variation 0
        # TODO: Remove start/finish diagonalization messages (back_box_numeric L:153)

        for key, value in results.items():
            if key == "hfss_variables":  # All numeric-only datatypes
                return
            value = np.array(value)
            corr_value = np.array(self.correct_res[key])

            self.assertTrue(np.allclose(value, corr_value))
            epr.logger.info(key + " " + "-" * (13 - len(key)) + "-> OK!")

    def test_analyze_variation(self):
        pass

    def test_hamiltonian(self):
        pass  # TODO: Need to pass **kwargs to epr_num_diag for return_H option

    def test_properties(self):
        pass
