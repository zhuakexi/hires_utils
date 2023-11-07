from itertools import dropwhile
from toolz import reduceby
import sys
import time
def cli(args):
    input_file, input_fmt, tag, output_file = \
        args.input_file, args.input_fmt, args.tag, args.output_file
    if output_file == input_file:
        raise ValueError("output_file must be different from input_file")
    if input_fmt == 'STDIN':
        res = count_tag(sys.stdin, tag)
    elif input_fmt == 'SAM':
        with open(input_file, 'r') as f:
            res = count_tag(f, tag)
    elif input_fmt == "SAMgz":
        import gzip
        with gzip.open(input_file, 'rt') as f:
            res = count_tag(f, tag)
    elif input_fmt == 'BAM':
        import pysam
        with pysam.AlignmentFile(input_file, 'rb', threads=6, check_sq=False) as f:
            # transform AlignedSegment to SAM format
            sam = map(lambda x: x.to_string(), f)
            res = count_tag(sam, tag)
    else:
        raise ValueError("input_fmt must be one of 'STDIN', 'SAM', 'SAMgz', 'BAM'")
    dump_tsv(res, output_file)
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
        rg, value = None, None
        for field in line.split('\t')[11:]:
            if field.startswith('RG:'):
                rg = field.split(':')[2]
            if field.startswith(tag):
                value = field.split(':')[2]
        if rg is None or value is None:
            print(line)
            raise ValueError("RG or tag not found in the SAM file")
        return rg, value
    def count(x, y):
        x += 1
        return x
    time_start = time.time()
    result = reduceby(RG_tag, count, sam, 0)
    print("count_tag finished, running time:%.2f min" % ((time.time()-time_start)/60))
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
