import gzip
import os
from hires_utils.countsam import count_tag
def test_count_tag(request):
    filename = os.path.join(request.fspath.dirname, "data", "FC.sam.gz")
    # get top 3 lines of FC.sam.gz
    with gzip.open(filename, 'rt') as f:
        sam = f.readlines()[:3]
    # test count_tag function
    assert count_tag(sam, "XS") == {
        ("Control_03","Unassigned_NoFeatures"):1,
        ("noEXO1_02", "Assigned"):1,
        ("noEXO1_05", "Unassigned_NoFeatures"):1
    }
def test_count_tag2(request):
    filename = os.path.join(request.fspath.dirname, "data", "FC.sam.gz")
    with gzip.open(filename, 'rt') as f:
        sam = f.readlines()
    # test count_tag function
    res = count_tag(sam, "AS")
    print(res)
    assert res[("Control_03","Unassigned_NoFeatures")] == 62