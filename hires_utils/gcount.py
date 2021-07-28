import gzip
from subprocess import check_output
from .hires_io import gen_record


def zcount(filename:str,target:str)->int:
    # count zipped file line number
    output_bytes = check_output([
        "zgrep","-c",target,
        filename])
    return int(output_bytes.decode("utf-8").strip())
def count_pairs(filename:str)->int:
    # count number of contacts in 4DN pairs file
    comments = 0
    with gzip.open(filename,"rt") as f:
        for line in f:
            if line[0] == "#":
                comments += 1
    all_lines = zcount(filename, "$")
    return all_lines - comments
def count_fastq(filename:str, form:str="p")->int:
    # count fastq file
    # Input:
    #    form: s(ingle_end), p(air_end), m(erged)
    #    filename: fastq file. for pair_end,
    #             give R1 or R2
    # Output:
    #    number of sequencing reads
    if form == "s" or form == "m":
        # in fastq, @ID starts one read
        return zcount(filename,"@")
    elif form == "p":
        # include the paired fastq file 
        return zcount(filename,"@") * 2
    else:
        return -1
def cli(args):
    filename, file_format, record_directory, sample_name = \
        args.filename[0], args.fformat, args.record_directory, args.sample_name
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
        if sample_name != None:
            record = {sample_name:result}
        else:
            record = {filename:result}
        gen_record(record, record_directory)
        