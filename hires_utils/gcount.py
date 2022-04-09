# count fastq, pairs files; gzipped or not
import gzip
from itertools import takewhile, repeat
from subprocess import check_output, CalledProcessError
import os
from .hires_io import gen_record, divide_name

def count_comments(filename:str, comark:str):
    """
    Count starting comment lines.
    Input:
        filename: path of gzipped file (ends with .gz or .gzip) or txt file(other extend names)
        comark: comment line starts with ...
    Output:
        number of comment lines; int
    """
    if os.path.splitext(filename)[1] in [".gz",".gzip"]:
        f = gzip.open(filename, "rt")
    else:
        f = open(filename, "rt")
    comments = 0
    for line in takewhile(lambda x: x.startswith(comark), f):
            comments += 1
    f.close()
    return comments
def big_count_line(filename):
    """    
    Count line number of file using raw.read and byte count trick.
    Detect extend name.
    Input:
        filename: path of gzipped file (ends with .gz or .gzip) or txt file(other extend names)
    Return:
        line number; int 
    """
    if os.path.splitext(filename)[1] in [".gz",".gzip"]:
        with gzip.open(filename, "rb") as f:
            bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
            line_num = sum( buf.count(b'\n') for buf in bufgen if buf)
    else:
        with open(filename, "rb") as f:
            bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
            line_num = sum( buf.count(b'\n') for buf in bufgen if buf)
    return line_num
def count_pairs(filename:str)->int:
    # count number of contacts in 4DN pairs file
    return big_count_line(filename) - count_comments(filename, "#")
def count_fastq(filename:str, form:str="p")->int:
    # count fastq file
    # Input:
    #    form: s(ingle_end), p(air_end), m(erged)
    #    filename: fastq file. for pair_end,
    #             give R1 or R2
    # Output:
    #    number of sequencing reads
    if form == "s" or form == "m":
        # in fastq, @ID starts one read, but quality score line has ^@ too.
        # popular solution is line_num/4
        return big_count_line(filename) // 4
    elif form == "p":
        # include the paired fastq file 
        return big_count_line(filename) // 4 * 2
    else:
        return -1
def cli(args):
    filename, file_format, record_directory, sample_name, attr = \
        args.filename[0], args.fformat, args.record_directory, args.sample_name, args.attributes
    if file_format in ["pe_fastq"]:
        result = count_fastq(filename)
    elif file_format in ["se_fastq","merged_fastq"]:
        result = count_fastq(filename,"m")
    elif file_format in ["pairs"]:
        result = count_pairs(filename)
    else:
        print("lcount: format not supported yet.")
        return -1
    print(filename + ":" + str(result))

    if record_directory != None:
        infer_sample_name, _ = divide_name(filename)
        if sample_name != None:
            if attr != None:
                record = {sample_name:{attr:result}}
        else:
            record = {infer_sample_name:{"reads":result}}
        gen_record(record, record_directory)
        