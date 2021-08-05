# regenerate sample_name
import sys
import numpy as np
import pandas as pd
import rmsd
import os
from itertools import combinations

from .hires_io import parse_3dg, divide_name, gen_record, print_records

def flip_rmsd(struct_a:np.ndarray, struct_b:np.ndarray)->np.float16:
    # calculate median deviation of 2 xyz array
    # aware of mirror effect
    a = rmsd.kabsch_rmsd(struct_a, struct_b, translate=True)
    b = rmsd.kabsch_rmsd(struct_a, - 1.0 * struct_b, translate=True)
    return a if a < b else b
def pick_good(rmsds:dict, threshold:float=2) -> set:
    # pick 3 structures from 3 smallest rmsd pairs
    # input:
    #   {[struct pair] : rmsd of the pair}
    # output:
    #   (final rmsd, name of good_structures)

    # calc smallest 3 mean rmsd
    s_rmsds = sorted(zip(rmsds.values(), rmsds.keys()))
    # for number of good struct >= 3, pick 3 struct
    # for number of good struct < 3, return None
    if s_rmsds[2][0] > threshold:
        good_struct= None
    else:
        pairs = [rec[1] for rec in s_rmsds[0:3]]
        good_struct = set([i for pair in pairs for i in pair])
    # using smallest 3 mean as final rmsd anyway
    final_RMSD = np.mean([rec[0] for rec in s_rmsds[0:3]])
    return final_RMSD, good_struct   
def RMS(data:list):
    return np.sqrt((data ** 2).mean())
def combine_binary(samples:list, data:dict, func)->dict:
    # calc all possible 2-pair's from sample list
    # input:
    #   samples: simple name list
    #   data: using name in *samples* as index, 
    #       store real input for *func*
    #   func: binary function
    # output:
    #   { sample_pair(tuple) : func return }
    res = {}
    for pair in combinations(samples, 2):
        res[pair] = func(data[pair[0]], data[pair[1]])
    return res
def build_records(sample_name, attr, final_RMSD, per_pair_rmsds, good_files):
    # [for pipeline]
    # generate records of chrom_rmsd's result
    # using a 2-layer dict model
    # mix attr in
    # attr is typicaly 3d build binsize
    keys = [attr + "rmsd"]
    values = [final_RMSD]
    for pair in per_pair_rmsds:
        keys.append(attr+"_".join(pair)+"_rmsd")
        values.append(per_pair_rmsds[pair])
    for struct_name in good_files:
        keys.append(attr + struct_name)
        values.append(good_files[struct_name])
    return {sample_name:dict(zip(keys, values))}
def cli(args):
    filenames, result_log, record_dir, sample_name, attr = \
        args.filenames, args.result_log, args.record_dir, args.sample_name, args.attr
    final_RMSD, per_pair_rmsds, good_files = chrom_rmsd(filenames)
    if sample_name == None:
        # infering sample name
        sample_name,_ = divide_name(filename)
    if attr == None:
        attr = ""
    records = build_records(
        sample_name, attr, final_RMSD,
        per_pair_rmsds, good_files)
    if result_log == None:
    # print to stdout
        print_records(records)
    else:
        # print to log file
        with open(result_log,"wt") as f:
            print_records(records, f)
    # [pipeline] record for downstream analysis
    if record_dir != None:
        gen_record(records, record_dir)

def chrom_rmsd(filenames):
    # Input:
    #   filenames
    # Output:
    #   float : smallest-3-mean-r.m.s.deviation of structs
    #   dict : per-pair-r.m.s.dev of sturcts
    #   dict : (if exist) 3 good struct / blank dict
    # load 3dg file
    
    ## inner names/index of different structures
    names = [str(i) for i in range(len(filenames))]
    name_mapping = pd.Series(index=names, data=filenames)
    ## load all structures 
    ss = [parse_3dg(f) for f in name_mapping.values]
    ## get common particles of all structures, 
    ## eliminate NA
    cells = pd.concat(ss, axis=1, join="inner", keys=names, names=["cell","pos"])
    
    # calc all cell-pair's RMSD

    rmsds = combine_binary(names, cells, flip_rmsd)
    # pick good structure and
    # calc filan RMSD between good structures
    final_RMSD, good_struct = pick_good(rmsds)
    if good_struct != None:
        good_files = [os.path.abspath(file) for file in name_mapping[good_struct].values]
        keys = ["g_struct1", "g_struct2", "g_struct3"]
        good_files_d = dict(zip(keys,good_files))
    else:
        good_files_d = {}
    return final_RMSD, rmsds, good_files_d


    