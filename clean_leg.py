import pandas as pd
import numpy as np
import sys
import time
import gzip
from concurrent import futures
from functools import partial

'''
default 4DN .pairs format
READID, CHR1, POS1, CHR2, POS2, STRAND1, STRAND2 = 0,1,2,3,4,5,6
'''
def is_promiscuous(contact:"line", all_legs:dict)->bool:
    '''
    say if either of contact's leg is promiscuous by check leg dictionary
    Too slow! Only for small data set.
    '''
    hit_left = np.abs(all_legs[contact["CHR1"]]["POS"] - contact["POS1"]) < 1000
    hit_right = np.abs(all_legs[contact["CHR2"]]["POS"] - contact["POS2"]) < 1000 
    return (hit_left.astype(np.int8).sum() > 10) or (hit_right.astype(np.int8).sum() > 10)
def clean_promiscuous(contacts:"dataframe", all_legs:dict, thread:int)->"dataframe":
    mask = contacts.apply(is_promiscuous, axis=1, all_legs=all_legs)
    sys.stderr.write("clean_leg: 1/%d chunk done\n"%thread)
    return contacts[~mask]
def clean_leg_main(args):
    cell_name, num_thread, replace, out_name = args.filenames[0], args.thread, args.replace_switch, args.out_name
    t0 = time.time()
    #read data file
    with gzip.open(cell_name, "rt") as f:
        #get file head
        head = [next(f) for i in range(0,27)]
        head = "".join(head)
    cell = pd.read_table(cell_name, header=None, skiprows=27)
    pairs_column = "READID CHR1 POS1 CHR2 POS2 STRAND1 STRAND2 PHASE0 PHASE1".split()
    cell.columns = pairs_column
    sys.stderr.write("clean_leg: parsing data file in %.2fs\n"%(time.time()-t0))
    #merge left and right legs, hash by chromosome_names
    t0 = time.time()
    left, right = cell[["CHR1","POS1"]], cell[["CHR2","POS2"]]
    left.columns, right.columns = ("CHR","POS"),("CHR","POS")
    all_legs = pd.concat((left,right),axis=0,ignore_index=True)
    all_legs = {key:value for key, value in all_legs.groupby("CHR")}
    sys.stderr.write("clean_leg: calculating in %.2fs\n"%(time.time()-t0))
    #multithread filtering
    t0=time.time()
    input = np.array_split(cell, num_thread, axis=0)
    working_func = partial(clean_promiscuous,all_legs=all_legs, thread=num_thread)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input)
    result = pd.concat(res,axis=0)
    sys.stderr.write("clean_leg: cleaning finished in %.2fs\n"%(time.time()-t0))
    if replace == True:
        out_name = cell_name
    with gzip.open(out_name,"wt") as f:
        f.write(head)
    result.to_csv(out_name, sep="\t", \
        header=False, index=False, mode="a")
    