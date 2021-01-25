import time
import sys
import gzip
import os
import pandas as pd

from classes import Cell, Data

'''
work with local pairs file generated from hickit
has 27 headline start with "#"
'''
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 
def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    #read comment line
    with gzip.open(filename,"rt") as f:
        comment = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comment.append(line)
    #read table format data
    pairs = pd.read_table(filename, header=None, comment="#")
    pairs.attrs["comment"] = comment
    #assign column names
    pairs.columns = name_array[0:pairs.shape[1]]
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def write_pairs(pairs:pd.DataFrame, out_name:str):
    '''
    write dataframe to tab delimited zipped file
    reserve comment lines, no dataframe index and headers
    need to change comment line with data's real column names
    '''
    #sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        f.write(pairs.attrs["comment"])
        pairs.to_csv(f, sep="\t", header=False, index=False, mode="a")