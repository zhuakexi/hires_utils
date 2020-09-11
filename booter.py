import argparse
#only for error class
import sys
from concurrent import futures
import os

from hires_io import parse_pairs, write_pairs, queue_read, queue_write
from classes import Cell

def parallel(handle:"one cell argument", filenames:list, out_name:str, replace:bool)->int:
    '''
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
    '''
    pass

def batch(handle:"one cell argument", filenames:list, out_name:str, replace:bool)->int:
    pass
def multi(handle:"one cell argument", filenames:"multi elements list", out_name:str, replace:bool, num_thread=1)->int:
    #call handle in several times
    cells = queue_read(filenames)
    if num_thread == 1:
        #function doesn't have multi thread
        with futures.ProcessPoolExecutor() as pool:
            res = pool.map(handle, cells)
    else:
        res = (handle(cell) for cell in cells)
    queue_write(res, out_name, replace, filenames)
def single(handle:"one cell argument", filenames:"one element list", out_name:str, replace:bool)->int:
    #call handle once
    real_in_name = filenames[0]
    cell = parse_pairs(real_in_name)
    res = handle(cell)
    if replace == True:
        #outname = ""
        write_pairs(res, real_in_name)
    else:
        #not inplace, using out_name
        if os.path.isdir(out_name):
            #out_name is directory
            write_pairs(res, os.path.join(out_name, res.name + res.appendix))
        else:
            write_pairs(res, out_name)