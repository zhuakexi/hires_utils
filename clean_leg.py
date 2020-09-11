import pandas as pd
import numpy as np
import sys
import time
import gzip
import bisect
import argparse
from concurrent import futures
from functools import partial
from collections import namedtuple

from classes import Cell
from hires_io import parse_pairs, write_pairs
import booter
'''
default 4DN .pairs format
READID, chr1, pos1, chr2, pos2, STRAND1, STRAND2 = 0,1,2,3,4,5,6
'''
def is_leg_promiscuous(leg, sorted_legs:dict, max_distance, max_count):
    '''
    tell if one leg is promiscuous
    using Tan's phantom leg method
    '''
    this_chromosome = sorted_legs[leg.chr]
    left_end = bisect.bisect_left(this_chromosome["pos"], leg.pos - max_distance)
    right_stretch = left_end + max_count
    if right_stretch >= len(this_chromosome):
        return False
    return this_chromosome.iloc[right_stretch]["pos"] - leg.pos <= max_distance 
def is_promiscuous(contact:"line", sorted_legs:dict, max_distance, max_count)->bool:
    '''
    tell if one contact is promiscuous
    '''
    Leg = namedtuple("Leg", "chr pos")
    leg1, leg2 = Leg(contact["chr1"], contact["pos1"]), Leg(contact["chr2"], contact["pos2"])
    hit = partial(is_leg_promiscuous, sorted_legs=sorted_legs, max_distance=max_distance, max_count=max_count)
    return hit(leg1) or hit(leg2)
def clean_promiscuous(contacts:"dataframe", sorted_legs:dict, thread:int, max_distance, max_count)->"dataframe":
    #wrapper for multi-thread processing
    mask = contacts.apply(is_promiscuous, axis=1, sorted_legs=sorted_legs, max_distance=max_distance, max_count=max_count)
    sys.stderr.write("clean_leg: 1/%d chunk done\n"%thread)
    return contacts[~mask]
def cli(args):
    filenames, num_thread, replace, out_name, max_distance, max_count, batch_switch, parallel_switch = \
        args.filenames, args.thread, args.replace_switch, args.out_name, args.max_distance, args.max_count, args.batch_switch, args.parallel_switch
    working_function = partial(clean_leg, num_thread=num_thread, max_count=max_count, max_distance=max_distance)
    '''
    if parallel_switch == True:
        return booter.parallel(working_function, filenames, out_name, replace)
    if batch_switch == True:
        return booter.batch(clean_leg, filenames, out_name, replace, num_thread)
    '''
    if not parallel_switch and not batch_switch:
        if len(filenames) > 1:
            return booter.multi(working_function, filenames, out_name, replace, num_thread)
        else:
            return booter.single(working_function, filenames, out_name, replace)
def clean_leg(cell:Cell, num_thread:int, max_distance:int, max_count:int):
    #merge left and right legs, hash by chromosome_names
    t0 = time.time()
    left, right = cell.data[["chr1","pos1"]], cell.data[["chr2","pos2"]]
    left.columns, right.columns = ("chr","pos"), ("chr","pos")
    all_legs = pd.concat((left,right), axis=0, ignore_index=True)
    sorted_legs = {key:value.sort_values(by="pos",axis=0,ignore_index=True) for key, value in all_legs.groupby("chr")}
    sys.stderr.write("clean_leg: group sort in %.2fs\n"%(time.time()-t0))
    #multithread filtering
    t0=time.time()
    input = np.array_split(cell.data, num_thread, axis=0)
    working_func = partial(clean_promiscuous,sorted_legs=sorted_legs, 
                           thread=num_thread, max_distance=max_distance, max_count=max_count)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input)
    result = pd.concat(res, axis=0)
    print("clean_leg: remove %d contacts in %s\n"%(len(cell.data)-len(result), cell.name))
    sys.stderr.write("clean_leg: finished in %.2fs\n"%(time.time()-t0))
    cell.data=result
    cell.appendix = ".pairs.gz" #keep old one
    return cell

    