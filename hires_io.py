import time
import sys
import gzip
import pandas as pd

'''
work with local pairs file generated from hickit
has 27 headline start with "#"
'''
def pairs_parser(cell_name:str)->"dataframe":
    '''
    read from 4DN's .pairs format
    '''
    t0 = time.time()
    column_names = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1".split()
    pairs = pd.read_table(cell_name, header=None,skiprows=27)
    pairs.columns = column_names
    sys.stderr.write("pairs_parser parsing time: %.2fs\n"%(time.time()-t0))
    return pairs
def write_pairs(data:"dataframe",in_name:str, out_name:str):
    '''
    write dataframe to tab delimited zipped file
    reserve heads, no dataframe index and headers
    in_file == out_file will replace original file
    '''
    with gzip.open(in_name, "rt") as f:
        #get file head
        head = [next(f) for i in range(0,27)]
        head = "".join(head)
    with gzip.open(out_name,"wt") as f:
        f.write(head)
        data.to_csv(out_name, sep="\t", header=False, index=False, mode="a")