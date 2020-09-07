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

from hires_io import pairs_parser
from hires_io import write_pairs
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
    filenames, num_thread, replace, out_name, max_distance, max_count, batch_switch = \
        args.filenames, args.thread, args.replace_switch, args.out_name, args.max_distance, args.max_count, args.batch_switch
    if len(filenames) > 1:
        for cell_name in filenames:
            if replace == True:
                #--replace will work in multi file input
                the_out_name = cell_name
            else:
                #--out_name will be used as name appendix: xx.appendix.pairs.gz 
                the_out_name = cell_name.split(".")
                the_out_name.insert(1,out_name)
                the_out_name = ".".join(the_out_name)
            return clean_leg_main(cell_name, num_thread, the_out_name, max_distance, max_count)
    if batch_switch == True:
        if len(filenames) > 1:
            raise argparse.ArgumentError(None, "in batch mode only one name list file is permitted.")
        #in batch mode, filenames is parsed a name list file, by line
        with open(filenames[0]) as f:
            filenames = [line.strip() for line in f.readlines()]
        if out_name[-1] != "/":
            #not directory but out name list file
            with open(out_name) as f:
                outlist = [line.strip() for line in f.readlines()]
            if len(outlist) != len(filenames):
                raise argparse.ArgumentError(None, "using out name list now: outlist must be same length with inlist.")
            else:
                for i, cell_name in enumerate(filenames):
                    the_out_name = outlist[i]
                    return clean_leg_main(cell_name, num_thread, the_out_name, max_distance, max_count)
        for cell_name in filenames:
            if replace == True:
                the_out_name = cell_name
            else:
                #using --outname as out directory
                the_out_name = out_name + cell_name.split("/")[-1]
            return clean_leg_main(cell_name, num_thread, the_out_name, max_distance, max_count)
    #neither multi filenames nor batch mode
    if replace == True:
        #--replace in single file case
        the_out_name = cell_name    
    return clean_leg_main(cell_name, num_thread, the_out_name, max_distance, max_count)
def clean_leg_main(cell_name, num_thread, out_name, max_distance, max_count):
    t0 = time.time()
    #read data file
    cell = pairs_parser(cell_name)
    #merge left and right legs, hash by chromosome_names
    t0 = time.time()
    left, right = cell[["chr1","pos1"]], cell[["chr2","pos2"]]
    left.columns, right.columns = ("chr","pos"),("chr","pos")
    all_legs = pd.concat((left,right),axis=0,ignore_index=True)
    sorted_legs = {key:value.sort_values(by="pos",axis=0,ignore_index=True) for key, value in all_legs.groupby("chr")}
    sys.stderr.write("clean_leg: group sort in %.2fs\n"%(time.time()-t0))
    #multithread filtering
    t0=time.time()
    input = np.array_split(cell, num_thread, axis=0)
    working_func = partial(clean_promiscuous,sorted_legs=sorted_legs, 
                           thread=num_thread, max_distance=max_distance, max_count=max_count)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input)
    result = pd.concat(res,axis=0)
    sys.stderr.write("clean_leg: remove %d contacts\n"%(len(cell)-len(result)))
    sys.stderr.write("clean_leg: cleaning finished in %.2fs\n"%(time.time()-t0))
    write_pairs(result, cell_name, out_name)
    return result

    