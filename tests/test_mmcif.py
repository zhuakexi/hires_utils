import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
from hires_utils.mmcif import threedg_to_cif

class TestClipBPymol(unittest.TestCase):
    def setUp(self):
        # Setup can include creating test data or configurations common to all methods
        self.input_file = "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg"
        self.output_dir = "/share/home/ychi/dev/hires_utils/tests/output"
        self.cpg_hom_file = "/share/home/ychi/software/dip-c/color/hg19.cpg.20k.hom.txt"
        self.cpg_file = "/share/home/ychi/software/dip-c/color/hg19.cpg.20k.txt"
        self.output_file_diploid = os.path.join(self.output_dir,"GMO1001.clean.20k.4.cpg.dip.cif")
        self.output_file_haploid = os.path.join(self.output_dir,"GMO1001.clean.20k.4.cpg.hap.png")

    def test_threedg_to_cif_chimera(self):
        threedg_to_cif(self.input_file, self.output_file_diploid, self.cpg_hom_file,
            maxGap=1000000, flavor="chimera", resolution=20000)
        self.assertTrue(os.path.exists(self.output_file_diploid))

    def tearDown(self):
        # Clean up code if necessary, e.g., deleting test files, etc.
        pass

if __name__ == '__main__':
    unittest.main()
