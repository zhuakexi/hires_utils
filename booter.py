import argparse
#only for error class
import sys
from concurrent import futures
import os

from hires_io import parse_pairs, write_pairs, queue_read, queue_write
from classes import Cell

def parallel(handle:"one cell argument", filenames:list, out_name:str, replace:bool)->int:
    pass

def batch(handle:"one cell argument", filenames:list, out_name:str, replace:bool)->int:
    pass
'''
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
'''
def multi(handle:"one cell argument", input:"list_of_cell", pool_size:int=1)->"list_of_cell":
    #call handle in several times
    with futures.ProcessPoolExecutor(pool_size) as pool:
        res = pool.map(handle, input)
    return res
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

booter = {"parallel":parallel, "batch":batch, "multi":multi, "single":single}