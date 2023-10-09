from itertools import dropwhile
from toolz import reduceby
def count_tag(sam, tag):
    """
    Count the number of reads in a SAM file that have a given tag within each read group.
    Input:
        sam - a SAM file
        tag - a tag to count
    Output:
        a dictionary of read (group,value) : count
    """
    result = {}
    sam = dropwhile(lambda x: x.startswith('@'), map(str.strip, sam))
    def RG_tag(line):
        # return the read group and the value of the tag
        for field in line.split('\t')[11:]:
            if field.startswith('RG:'):
                rg = field.split(':')[2]
            if field.startswith(tag):
                value = field.split(':')[2]
        return rg, value
    def count(x, y):
        x += 1
        return x
    result = reduceby(RG_tag, count, sam, 0)
    return result
# import gzip
# filename = "/share/home/ychi/dev/hires_utils/tests/data/FC.sam.gz"
# with gzip.open(filename, 'rt') as f:
#         sam = f.readlines()
# print(count_tag(sam, "XS"))