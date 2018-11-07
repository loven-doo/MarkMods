import os
import json
import unittest

import numpy as np

from markmods.constants import TEST_DIR
from markmods.models.hmm import base, sdhmm, mdhmm


PACKAGE_DIR = 'models_hmm'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestBase(unittest.TestCase):
    pass


class TestSDHMM(unittest.TestCase):

    data_array = np.array([
        {"s": "a"},
        {"s": "a"},
        {"s": "b"},
        {"s": "a"},
        {"s": "b"},
        {"s": "c"},
        {"s": "b"},
        {"s": "b"},
        {"s": "a"},
        {"s": "d"},
    ])
    with open(os.path.join(INPUT_DIR, "sdhmm_scheme.json")) as hmm_scheme_json:
        hmm_scheme = json.load(hmm_scheme_json)  # 1D HMM configuration
    # Let's add emission functions to configuration dict (it can be any function)
    hmm_scheme["states"]["state1"]["emis_func"] = lambda x: {"a": 0.99, "b": 0.01, "c": 0.003, "d": 0.02}[x["s"]]
    hmm_scheme["states"]["state1"]["emis_func"] = lambda x: {"a": 0.01, "b": 0.99, "c": 0.003, "d": 0.02}[x["s"]]
    hmm_scheme["states"]["state3"]["emis_func"] = lambda x: {"a": 0.01, "b": 0.1, "c": 0.3, "d": 0.01}[x["s"]]
    hmm_scheme["states"]["state4"]["emis_func"] = lambda x: {"a": 0.01, "b": 0.05, "c": 0.15, "d": 0.6}[x["s"]]
    single_dim_hmm = sdhmm.SDHMM(scheme=hmm_scheme)

    def test_init(self):
        states_trans_exp = {
            'state1': {
                'state1': 0.95,
                'state2': 0.2,
                'state3': 0.9,
                'state4': 0.1,
            },
            'state2': {
                'state1': 0.2,
                'state2': 0.9,
                'state3': 0.2,
                'state4': 0.05,
            },
            'state3': {
                'state1': 0.9,
                'state2': 0.2,
                'state3': 0.88,
                'state4': 0.1,
            },
            'state4': {
                'state1': 0.8,
                'state2': 0.05,
                'state3': 0.8,
                'state4': 0.25,
            },
        }
        states_trans = dict()
        for state_name in self.single_dim_hmm.states:
            states_trans[state_name] = self.single_dim_hmm.states[state_name].trans
        self.assertEqual(states_trans_exp, states_trans)


class TestMDHMM(unittest.TestCase):

    grid = np.array([  # The element of a matrix may be any type that your emission func can prepare
        [
            [{"s_set": ("a","b","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","a")}, {"s_set": ("a","a","a")}],
            [{"s_set": ("a","a","b")}, {"s_set": ("b","b","b")}, {"s_set": ("a","-","a")}, {"s_set": ("a","a","a")}],
            [{"s_set": ("a","b","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","b","a")}, {"s_set": ("a","a","b")}],
        ],
        [
            [{"s_set": ("a","a","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","b")}, {"s_set": ("a","b","a")}],
            [{"s_set": ("a","a","a")}, {"s_set": ("b","b","b")}, {"s_set": ("b","-","a")}, {"s_set": ("a","a","a")}],
            [{"s_set": ("a","b","-")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","a")}, {"s_set": ("a","a","b")}],
        ],
    ])
    with open(os.path.join(INPUT_DIR, "mdhmm_scheme.json")) as hmm_scheme_json:
        hmm_scheme = json.load(hmm_scheme_json)  # 3D HMM configuration
    # Let's add emission functions to configuration dict (it can be any function)
    hmm_scheme["states"]["state1"]["emis_func"] = lambda x: {
        frozenset(("a", "b",)): 0.01,
        frozenset(("a",)): 0.6,
        frozenset(("b",)): 0.7,
        frozenset(("-",)): 0.1,
        frozenset(("a", "b", "-")): 0.01,
        frozenset(("a", "-")): 0.99,
        frozenset(("b", "-")): 0.99,
    }[frozenset(x["s_set"])]
    hmm_scheme["states"]["state2"]["emis_func"] = lambda x: {
        frozenset(("a", "b",)): 0.8,
        frozenset(("a",)): 0.1,
        frozenset(("b",)): 0.1,
        frozenset(("-",)): 0.2,
        frozenset(("a", "b", "-")): 0.45,
        frozenset(("a", "-")): 0.1,
        frozenset(("b", "-")): 0.2,
    }[frozenset(x["s_set"])]
    hmm_scheme["states"]["state3"]["emis_func"] = lambda x: {
        frozenset(("a", "b",)): 0.45,
        frozenset(("a",)): 0.1,
        frozenset(("b",)): 0.1,
        frozenset(("-",)): 0.01,
        frozenset(("a", "b", "-")): 0.8,
        frozenset(("a", "-")): 0.1,
        frozenset(("b", "-")): 0.2,
    }[frozenset(x["s_set"])]
    mult_dim_hmm = mdhmm.MDHMM(scheme=hmm_scheme)

    def test_init(self):
        states_trans_exp = {
            'state1': {},
            'state2': {},
            'state3': {},
        }
        states_trans = dict()
        for state_name in self.mult_dim_hmm.states:
            states_trans[state_name] = self.mult_dim_hmm.states[state_name].trans
        self.assertEqual(states_trans_exp, states_trans)


if __name__ == "__main__":
    unittest.main()
