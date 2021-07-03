import sys
import numpy as np
import rmsd
import os
from itertools import combinations

from hires_io import parse_3dg

def flip_rmsd(struct_a:np.ndarray, struct_b:np.ndarray)->np.float16:
    # calculate median deviation of 2 xyz array
    # aware of mirror effect
    a = rmsd.kabsch_rmsd(struct_a, struct_b, translate=True)
    b = rmsd.kabsch_rmsd(struct_a, - 1.0 * struct_b, translate=True)
    return a if a < b else b
def pick_good(rmsds:dict, threshold:float=2) -> set:
    #pick out structure dv bigger than THRESHOLD
    # input:
    #   {[struct pair] : rmsd of the pair}
    # output:
    #    name of good_structures
    ##get all structure name(represent by int index)
    problematic = set()
    for pair in pairs:
        problematic = problematic | set(pair)
        good = problematic | set(pair)
    ##pick bad structure(clustering here may be better)
    for pair in rmsds:
        if rmsds[pair] < threshold:
            problematic -= set(pair)
    ##remaining is the good
    good -= problematic
    return good
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
def cli(args):
    filenames, result_log = args.filenames, args.result_log 
    chrom_rmsd(filenames, result_log)
def chrom_rmsd(filenames, result_log):
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
    good = list(pick_good(rmsds))
    if len(good) < 3:
        final_RMSD = RMS(rmsds)
    else:
        good_rmsds = np.array(list(
                combine_binary(good, 
                        cells, 
                        flip_rmsd).values()
                            ))
        final_RMSD = RMS(good_rmsds)

    # output

    with open(result_log, "wt") as f:
        f.write("#" + str(final_RMSD) + "\n")
        f.writelines("\n".join(good_files))

    