import os
import json
import unittest

from markmods.constants import TEST_DIR
from markmods.models import base


PACKAGE_DIR = 'models'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestBase(unittest.TestCase):

    exp_grouped_array = [[[1, 2, 3], [4, 5], [6]],
                         [[7, 8], [9, 10, 11, 12]]]
    exp_aggr_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    exp_group_lims = [[2, 4, 5, 7, 11], [2, 4]]

    def test_aggregate_array(self):
        aggr_array, group_lims = base.aggregate_array(self.exp_grouped_array, aggr_level=2)
        self.assertEqual(self.exp_aggr_array, aggr_array)
        self.assertEqual(self.exp_group_lims, group_lims)

    def test_group_array(self):
        grouped_array = base.group_array(self.exp_aggr_array, group_lims=self.exp_group_lims)
        self.assertEqual(self.exp_grouped_array, grouped_array)


if __name__ == "__main__":
    unittest.main()
