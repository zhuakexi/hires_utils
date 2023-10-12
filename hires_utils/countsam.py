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
def dump_tsv(res, output):
    """
    Dump the result of count_tag to a TSV file.
    Input:
        res - dict, key: (read group, tag value), value: count
        output - the output file name
    """
    with open(output, 'w') as f:
        for key, value in res.items():
            f.write('\t'.join([key[0], key[1], str(value)]) + '\n')
# import gzip
# filename = "/share/home/ychi/dev/hires_utils/tests/data/FC.sam.gz"
# with gzip.open(filename, 'rt') as f:
#         sam = f.readlines()
# dump_tsv(count_tag(sam, "XS"), "out/count_tag_FC.tsv")
# # print(count_tag(sam, "XS"))