import os
from numpy.testing import assert_allclose
from hires_utils.seg_stat import seg_values
def test_seg_values(request):
    filename = os.path.join(request.fspath.dirname, "data", "test.seg.gz")
    hap1_phased, hap2_phased, biasedX_score, hap_score, yp, xp, cp_count = seg_values(filename)
    assert_allclose(
        [hap1_phased, hap2_phased, biasedX_score, hap_score, yp],
        [0.26067188684011117,0.2157110381409447,1.0,0.09437963944856839,0.0010103561505430665]
    )