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
        {"s": "b"},
        {"s": "b"},
        {"s": "a"},
    ])
    with open(os.path.join(INPUT_DIR, "sdhmm_scheme.json")) as hmm_scheme_json:
        hmm_scheme = json.load(hmm_scheme_json)  # 1D HMM configuration
    # Let's add emission functions to configuration dict
    hmm_scheme["states"]["state1"]["emis_func"] = lambda x: {"a": 0.99, "b": 0.01}[x["s"]]  # it can be any function
    hmm_scheme["states"]["state1"]["emis_func"] = lambda x: {"a": 0.01, "b": 0.99}[x["s"]]
    hmm_scheme["states"]["state3"]["emis_func"] = lambda x: 0.1
    hmm_scheme["states"]["state4"]["emis_func"] = lambda x: 0.2
    single_dim_hmm = sdhmm.SDHMM(scheme=hmm_scheme)


class TestMDHMM(unittest.TestCase):

    grid = np.array([  # The element of a matrix may be any type that your emission func can prepare
        [
            [{"s_set": ("a","b","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","a")}, {"s_set": ("a","a","a")}],
            [{"s_set": ("a","a","b")}, {"s_set": ("b","b","b")}, {"s_set": ("a","a","a")}, {"s_set": ("a","a","a")}],
            [{"s_set": ("a","b","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","b","a")}, {"s_set": ("a","a","b")}],
        ],
        [
            [{"s_set": ("a","a","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","b")}, {"s_set": ("a","b","a")}],
            [{"s_set": ("a","a","a")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","a")}, {"s_set": ("a","a","a")}],
            [{"s_set": ("a","b","b")}, {"s_set": ("b","b","b")}, {"s_set": ("b","a","a")}, {"s_set": ("a","a","b")}],
        ],
    ])
    with open(os.path.join(INPUT_DIR, "mdhmm_scheme.json")) as hmm_scheme_json:
        hmm_scheme = json.load(hmm_scheme_json)  # 3D HMM configuration
    # Let's add emission functions to configuration dict
    hmm_scheme["states"]["state1"]["emis_func"] = lambda x: {
        frozenset(("a", "b")): 0.01,
        frozenset(("a",)): 0.99,
        frozenset(("b",)): 0.99,
    }[frozenset(x["s_set"])]  # it can be any function
    hmm_scheme["states"]["state2"]["emis_func"] = lambda x: 0.1
    hmm_scheme["states"]["state3"]["emis_func"] = lambda x: 0.2
    mult_dim_hmm = mdhmm.MDHMM(scheme=hmm_scheme)


if __name__ == "__main__":
    unittest.main()
