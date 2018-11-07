import os
import json
import unittest

from markmods.constants import TEST_DIR
from markmods import general


PACKAGE_DIR = 'markmods'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestGeneral(unittest.TestCase):

    def test_combine_parts(self):
        to_combine = [
            ["abc", "def", "gh"],
            ["ij", "klm"],
            ["edf", "abc", "vw"]
        ]
        combined_exp = [
            ["abc", "edf", "ij"],
            ["abc", "ij", "vw"],
            ["abc", "edf", "klm"],
            ["abc", "klm", "vw"],
            ["def", "edf", "ij"],
            ["abc", "def", "ij"],
            ["def", "ij", "vw"],
            ["def", "edf", "klm"],
            ["abc", "def", "klm"],
            ["def", "klm", "vw"],
            ["edf", "gh", "ij"],
            ["abc", "gh", "ij"],
            ["gh", "ij", "vw"],
            ["edf", "gh", "klm"],
            ["abc", "gh", "klm"],
            ["gh", "klm", "vw"],
        ]
        self.assertEqual(set(tuple(comb) for comb in combined_exp),
                         set(tuple(comb) for comb in general.combine_parts(to_combine, without_repeats=True, sort="+")))


if __name__ == "__main__":
    unittest.main()
