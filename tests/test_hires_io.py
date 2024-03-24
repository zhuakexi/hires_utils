import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")

import pandas as pd
from hires_utils.hires_io import parse_3dg, parse_pairs

class Test_hires_io(unittest.TestCase):
    def setUp(self):
        self.pairs_path = os.path.join(os.path.dirname(__file__), "data", "test.pairs.gz")
        self._3dg_path = os.path.join(os.path.dirname(__file__), "data", "test.3dg.gz")
    def test_parse_3dg(self):
        # test basic functionality of parse_3dg
        _3dg = parse_3dg(self._3dg_path)
        print(_3dg.head())
        self.assertIsInstance(_3dg, pd.DataFrame)
    def test_parse_3dg_s2m(self):
        # test basic functionality of parse_3dg
        _3dg = parse_3dg(self._3dg_path, s2m=True)
        print(_3dg.head())
        self.assertIsInstance(_3dg, pd.DataFrame)
    def test_parse_pairs(self):
        pairs = parse_pairs(
            "/shareb/ychi/repo/sperm_struct/ds_pipeline/smk/new_clean/c2m5_c2d3000000/pairs_c12/HuS08_HuSS3078.c12.pairs.gz"
        )
        print(pairs.head())
        self.assertIsInstance(pairs, pd.DataFrame)
if __name__ == "__main__":
    unittest.main()