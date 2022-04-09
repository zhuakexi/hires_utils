from fileinput import filename
from itertools import count
import os
from hires_utils.gcount import count_pairs, count_fastq
def test_count_pairs(request):
    filename = os.path.join(request.fspath.dirname, "data", "test.pairs.gz")
    assert count_pairs(filename) == 1973
def test_count_fastq(request):
    filename = os.path.join(request.fspath.dirname, "data", "test_R1.fq.gz")
    # test.fq.gz has 50 reads, for default pair-end mode, library has 100 total reads
    assert count_fastq(filename) == 100