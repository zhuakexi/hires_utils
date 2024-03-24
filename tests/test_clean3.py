import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
from io import StringIO

import pandas as pd
import numpy as np

from hires_utils.clean3 import clean3, get_legs, parse_3dg, parse_pairs
from hires_utils.hires_io import m2s_index

class TestClean3(unittest.TestCase):
    def setUp(self):
        self.clean_quantile = 0.26
        self.max_clean_distance = 10000
        
        self.mock_s_data = {
            "chr": ["chr1", "chr1", "chr1", "chr1"],
            "pos": [20000, 40000, 60000, 80000],
            "x": [10, 20, 30, 40],
            "y": [40, 50, 60, 70],
            "z": [70, 80, 90, 100]
        }
        self.mock_s = pd.DataFrame(self.mock_s_data)
        
        # will remove this when treat as middle
        bin1_pos = [29999, 29998]
        bin2_pos = [49999, 49998, 49997]
        # will remove this when treat as start
        bin3_pos = [50001, 50002, 50003, 50004]
        bin4_pos = [89999, 89998, 89997]
        chr1_pos = bin1_pos + bin2_pos + bin3_pos + bin4_pos
        chr2_pos = chr1_pos
        chr1 = ["chr1"] * len(chr1_pos)
        chr2 = ["chr1"] * len(chr2_pos)
        self.mock_pairs_data = {
            "readID" : ".",
            "chr1": chr1,
            "pos1": chr1_pos,
            "chr2": chr2,
            "pos2": chr2_pos
        }
        self.mock_pairs = pd.DataFrame(self.mock_pairs_data)
        
        self.s_name = StringIO()
        self.s_name.write("#unit: 0.1\n")
        self.s_name.write(self.mock_s.to_csv(sep="\t", index=False, header=False))
        self.s_name.seek(0)
        self.con_name = StringIO(
            self.mock_pairs.to_csv(sep="\t", index=False, header=False)
        )

    def test_clean3(self):
        # test basic functionality of clean3
        good_3dg = clean3(self.s_name, self.con_name, self.clean_quantile, self.max_clean_distance)

        self.assertIsInstance(good_3dg, pd.DataFrame)
        # old version, treat pos in 3dg as middle of bin
        # remain_pos = [40000, 60000, 80000]
        # new version, treat pos in 3dg as start of bin
        good_3dg = m2s_index(good_3dg)
        remain_pos = [20000, 40000, 80000]
        self.assertTrue(all(good_3dg.index.get_level_values("pos").isin(remain_pos)))
        # newer version, treat pos in 3dg as start of bin
    
    def test_clean3_real(self):
        _3dg_p = "/sharec/ychi/repo/sperm64_hg/3dg/HuS08_HuSS3095.20k.1.3dg"
        pairs_p = "/sharec/ychi/repo/sperm64_hg/pairs_c12/HuS08_HuSS3095.c12.pairs.gz"
        good_3dg = clean3(
            _3dg_p, pairs_p, self.clean_quantile, self.max_clean_distance
        
        )
        print(good_3dg)
        self.assertIsInstance(good_3dg, pd.DataFrame)

if __name__ == '__main__':
    unittest.main()