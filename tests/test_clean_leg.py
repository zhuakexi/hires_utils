import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
from io import StringIO

from hires_utils.clean_leg import clean_leg
from hires_utils.hires_io import parse_pairs, write_pairs

class TestCleanLeg(unittest.TestCase):
    def setUp(self):
        self.out_dir = os.path.join(os.path.dirname(__file__), "output")
        self.num_thread = 4
        self.max_distance = 1000
        self.max_count = 10
    def test_clean_leg(self):
        pairsp = "/sharec/ychi/repo/sperm64_hg/pairs_0/HuS08_HuSS3056.pairs.gz"
        outp = os.path.join(self.out_dir, "HuS08_HuSS3056.c1.pairs.gz")
        pairs = parse_pairs(pairsp)
        res = clean_leg(
            pairs, self.num_thread, self.max_distance, self.max_count
            )
        write_pairs(res, outp)
if __name__ == "__main__":
    unittest.main()