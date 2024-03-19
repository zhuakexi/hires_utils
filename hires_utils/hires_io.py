import gzip
import json
import os
import sys
import time
import uuid
from functools import partial
from io import StringIO
from pathlib import Path
from pkgutil import get_data

import pandas as pd

from . import reference


'''
work with local pairs file generated from hickit
has 27 headline start with "#"
'''

# utils

## norm name
def converter_template(c_in:str,ref_dict:pd.DataFrame):
    # a reat_table converter function
    #print(ref_dict)
    if c_in in ref_dict:
        return ref_dict[c_in]
    else:
        return c_in
def fill_func_ref(template_func:callable, ref_file:str, index_col:str)->callable:
    # read in ref_file for template_fucn, generate new func
    # hope will boost new func's speed
    
    ## read in ref_file, get ref_dict in memory
    ref_df = pd.read_csv(ref_file, index_col=index_col)
    ref_dict = ref_df.iloc[:,0].to_dict()
    working_func = partial(template_func, ref_dict=ref_dict)
    return working_func
# converter: name -> standard name for pandas read_table 
norm_chr = fill_func_ref(
    converter_template,
    StringIO(
        ## get alias file in package
        ## reference is a "data module" with its own __init__.py
        get_data(
            reference.__name__,
            "chrom_alias.csv"
            ).decode()
        ),
    "alias"
    )
## get sample name
def divide_name(filename):
    if not (isinstance(filename, str) or isinstance(filename, Path)):
        return "", ""
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 

# parsers
def read_comments(file, next_line=False)->list:
    '''
    Read comments from file.
    Input:
        file: file path or StringIO
        next_line: whether to read next line after comments
    Output:
        list of comments
    '''
    comments = []
    if isinstance(file, StringIO):
        for line in file.getvalue().splitlines():
            if line[0] != "#":
                break
            comments.append(line.strip("#\n"))
    elif isinstance(file, str):
        try:
            with gzip.open(file,"rt") as f:
                for line in f:
                    if line[0] != "#":
                        break
                    comments.append(line.strip("#\n"))
        except OSError:
            with open(file,"rt") as f:
                for line in f:
                    if line[0] != "#":
                        break
                    comments.append(line.strip("#\n"))
    if next_line:
        return comments, line
    else:
        return comments
def parse_pairs(filename)->pd.DataFrame:
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    dtype_array = {"readID":"category",
            "chr1":"string",
            "pos1":"int",
            "chr2":"string",
            "pos2":"int",
            "strand1":"string",
            "strand2":"string",
            "phase0":"string",
            "phase1":"string",
            "phase_prob00":"float",
            "phase_prob01":"float",
            "phase_prob10":"float",
            "phase_prob11":"float"}
    #read comment line
    comments, line = read_comments(filename, next_line=True)
    #infer number of columns
    line_length = len(line.strip().split("\t"))
    #pick used eles from builtin arrays
    columns = name_array[0:line_length]
    dtypes = {key:value for key, value in dtype_array.items() if key in columns}
    #read table format data
    pairs = pd.read_table(
        filename, 
        header=None, 
        comment="#",
        dtype=dtypes,
        names=columns
        )
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # infer real sample name
    #assign column names
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def parse_gtf(filename:str) -> pd.DataFrame:
    # read gtf, get exons
    gencode = pd.read_table(filename, comment="#", header=None)
    gencode.columns="seqname source feature start end score strand frame group".split()
    return gencode
def parse_3dg(file:str, sorting=False)->pd.DataFrame:
    """
    Read in hickit 3dg file(or the .xyz file)
    Norm chr name alias
    Read into dataframe.attrs if has comments, treat last comment line as backbone_unit
    Input:
        filename: file path
        sorting: whether to sort chromosome and positions
    Output:
        3 col dataframe with 2-level multindex: (chrom, pos) x, y, z
    """

    ## read comments
    comments = read_comments(file, next_line=False)
    ## read real positions
    s = pd.read_table(file, 
                      comment="#",header=None,
                     index_col=[0,1],
                     converters={0:norm_chr})
    s.columns = "x y z".split()
    s.index.names = ["chr","pos"]
    ## assume last comment is backbone_unit
    if len(comments) > 0:
        s.attrs["comments"] = comments
        backbone_unit = float(comments[-1].split(":")[1].strip())
        s.attrs["backbone_unit"] = backbone_unit    
    if sorting:
        s.sort_index(inplace=True)
    return s
def parse_ref(filename:str, index=False, header=None, value_name="value") -> pd.DataFrame:
    """
    Parse reference file (csv file, with or without header or index, 
        only take first 3 non-index col), return a dataframe.
    Will convert homologous chromosomes to standard names.
    Input:
        filename: file path,
        index: whether the csv file has index col, will drop it
        header: give pd.read_table
        value_name: set column name for value col
    Output:
        pd.DataFrame with index: (chrom, pos) and column: value_name
    """
    ## read real positions
    ref = pd.read_table(
        filename,
        comment="#",
        converters={0:norm_chr},
        header=header
        )
    # treat index
    if index:
        ref = ref.iloc[:, 1:4]
    else:
        ref = ref.iloc[:, 0:3]
    ref.columns = ["chr", "pos", value_name]
    ref.set_index(["chr", "pos"], inplace=True)
    return ref
def parse_seg(filename:str) -> tuple:
    """
    read in .seg file in dip-c
    return comments and list of all segment element
    seg_format: read1/read1(. for read1, m for read2/mate), read_start, read_end, chrom, genome_start,
        genome_end, strand, haplotype(.for unknown, 0 for paternal, 1 for maternal)
    """
    comments = []
    segs = []
    with gzip.open(filename,"rt") as f:
        for line in f:
            if line[0] == "#":
                comments.append(line)
            else:
                for alignment in line.split()[1:]:
                    # remove read ID
                    segs.append(alignment)
    return comments, segs

# writers

def write_3dg(pairs:pd.DataFrame, outname:str):
    pairs.to_csv(outname, sep="\t",header=None)
    return 0
def write_pairs(pairs:pd.DataFrame, out_name:str):
    '''
    write dataframe to tab delimited zipped file
    reserve comment lines, no dataframe index
    headers store in last comment line
    need to sort with upper triangle label
    '''
    #sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        pairs.attrs["comments"].pop()
        pairs.attrs["comments"].append("#columns:" + "\t".join(pairs.columns) + "\n")
        f.write("".join(pairs.attrs["comments"]))
        pairs.to_csv(f, sep="\t", header=False, index=False, mode="a")

# [pipeline] json metrics recorder
def gen_record(record:dict, record_dir):
    uuid_str = uuid.uuid4().hex
    if not os.path.exists(record_dir):
        os.makedirs(record_dir) 
    with open(os.path.join(record_dir, uuid_str+".json"), "w") as f:
        json.dump(record, f)
def print_records(dat:dict,f=None):
    if f == None:
        for sample in dat:
            print(sample)
            for attr in dat[sample]:
                print(str(attr) + ":" + str(dat[sample][attr]))
    else:
        for sample in dat:
            print(sample,file=f)
            for attr in dat[sample]:
                print(str(attr) + ":" + str(dat[sample][attr]),file=f)
if __name__ == "__main__":
    from io import StringIO
    res = parse_3dg(
        "tests/data/test.3dg.gz",
        sorting=True
    )
    # suppose I wan't to switch x and y on the fly
    res = parse_3dg(
        StringIO(
            pd.read_table(
                "tests/data/test.3dg.gz",
                comment="#",
                header=None,
                names = "chr pos x y z".split()
            )[["chr","pos","y","x","z"]].to_csv(sep="\t",header=None,index=False)
            )
    )
    res = parse_3dg(
        "tests/data/test.3dg",
        sorting=True
    )