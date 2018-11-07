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
            ["ij", "KLM"],
            ["oprs", "tutu", "VW"]
        ]
        combined_exp = [
            "abc-ij-oprs",
            "abc-ij-tutu",
            "abc-ij-VW",
            "abc-KLM-oprs",
            "abc-KLM-tutu",
            "abc-KLM-VW",
            "def-ij-oprs",
            "def-ij-tutu",
            "def-ij-VW",
            "def-KLM-oprs",
            "def-KLM-tutu",
            "def-KLM-VW",
            "gh-ij-oprs",
            "gh-ij-tutu",
            "gh-ij-VW",
            "gh-KLM-oprs",
            "gh-KLM-tutu",
            "gh-KLM-VW",
        ]
        self.assertEqual(set(combined_exp), set(general.combine_str(to_combine, "-")))


if __name__ == "__main__":
    unittest.main()
