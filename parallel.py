import argparse
import sys
from concurrent import futures
def parallel(handle, filenames, out_name, replace):
    if len(filenames) > 1:
        raise argparse.ArgumentError(None, "in batch mode only one name list file is permitted.") 
    with open(filenames[0]) as f:
        filenames = [line.strip() for line in f.readlines()]
    #case1: out_name is not a directory, use it as name list, read and check length.
    if out_name[-1] != "/":
        with open(out_name) as f:
            outlist = [line.strip() for line in f.readlines()]
        if len(outlist) != len(filenames):
            raise argparse.ArgumentError(None, "using outname list now: outlist must be same length with inlist.")
        else:
            #outlist and inlist good, begin parallel
            with futures.ProcessPoolExecutor() as pool:
                res = pool.map(handle, zip(filenames, outlist))
            return res
    #case2: not outlist, use out_name as prefix or replace old file.
    if replace == True:
        input_data = zip(filenames, filenames)
    else:
        outlist = [out_name + cell_name.split("/")[-1] for cell_name in filenames]
        input_data = zip(filenames, outlist)
    with futures.ProcessPoolExecutor() as pool:
        res = pool.map(handle, input_data)
    return res