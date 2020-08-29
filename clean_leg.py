import pandas as pd
import numpy as np
import sys
import time
import gzip
from concurrent import futures
from functools import partial

from hires_io import pairs_parser
from hires_io import write_pairs
'''
default 4DN .pairs format
READID, chr1, pos1, chr2, pos2, STRAND1, STRAND2 = 0,1,2,3,4,5,6
'''
def is_promiscuous(contact:"line", all_legs:dict, max_distance, max_count)->bool:
    '''
    say if either of contact's leg is promiscuous by check leg dictionary
    Too slow! Only for small data set.
    '''
    #adjacent legs from the view of leg leg
    hit_left = np.abs(all_legs[contact["chr1"]]["pos"] - contact["pos1"]) < max_distance
    #adjacent legs from the view of right leg
    hit_right = np.abs(all_legs[contact["chr2"]]["pos"] - contact["pos2"]) < max_distance
    return (hit_left.astype(np.int8).sum() > max_count) or (hit_right.astype(np.int8).sum() > max_count)
def clean_promiscuous(contacts:"dataframe", all_legs:dict, thread:int, max_distance, max_count)->"dataframe":
    mask = contacts.apply(is_promiscuous, axis=1, all_legs=all_legs, max_distance=max_distance, max_count=max_count)
    sys.stderr.write("clean_leg: 1/%d chunk done\n"%thread)
    return contacts[~mask]
def clean_leg_main(args):
    cell_name, num_thread, replace, out_name, max_distance, max_count = \
        args.filenames[0], args.thread, args.replace_switch, args.out_name, args.max_distance, args.max_count
    t0 = time.time()
    #read data file
    cell = pairs_parser(cell_name)
    #merge left and right legs, hash by chromosome_names
    t0 = time.time()
    left, right = cell[["chr1","pos1"]], cell[["chr2","pos2"]]
    left.columns, right.columns = ("chr","pos"),("chr","pos")
    all_legs = pd.concat((left,right),axis=0,ignore_index=True)
    all_legs = {key:value for key, value in all_legs.groupby("chr")}
    sys.stderr.write("clean_leg: calculating in %.2fs\n"%(time.time()-t0))
    #multithread filtering
    t0=time.time()
    input = np.array_split(cell, num_thread, axis=0)
    working_func = partial(clean_promiscuous,all_legs=all_legs, thread=num_thread, max_distance=max_distance, max_count=max_count)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input)
    result = pd.concat(res,axis=0)
    sys.stderr.write("clean_leg: remove %d contacts\n"%(len(cell)-len(result)))
    sys.stderr.write("clean_leg: cleaning finished in %.2fs\n"%(time.time()-t0))
    if replace == True:
        out_name = cell_name
    write_pairs(result, cell_name, out_name)
    