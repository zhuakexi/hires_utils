import time
import sys
import gzip
import os
import pandas as pd

from classes import Cell, Data

'''
work with local pairs file generated from hickit
has 27 headline start with "#"
'''
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 
def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's .pairs format
    '''
    #column names are in the last comment line.
    #head are stored in cell, handler may change it.
    #now return Data
    with gzip.open(filename,"rt") as f:
        head = []
        #last_comment = None
        for line in f.readlines():
            if line[0] != "#":
                break
            head.append(line)
            #last_comment= line.strip("#\n")
    # unless you can validate .pairs, don't get col names from last comment
    # column_names = last_comment.split()[1:] 
    pairs = pd.read_table(filename, header=None,comment="#")
    pairs.columns = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1".split()
    sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    pairs_data = Data("pairs","".join(head), pairs, ".pairs", filename)
    name, extend = divide_name(filename)
    return Cell(name, pairs_data)
def write_pairs(cell:Cell, out_name:str):
    #now use data
    '''
    write dataframe to tab delimited zipped file
    reserve heads, no dataframe index and headers
    in_file == out_file will replace original file
    '''
    sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        f.write(cell.get_data("pairs").head)
        cell.get_data("pairs").content.to_csv(f, sep="\t", header=False, index=False, mode="a")
def queue_read(filenames:list)->"list of dataframe":
    '''
    read one by one. return generator
    '''
    #use different parser
    if len(filenames) == 1:
        filename = filenames[0]
        if os.path.isdir(filename):
            #use it as a directory
            return (parse_pairs(name) for name in os.listdir(filename))
        else:
            #use it as a namelist file
            with open(filename) as f:
                real_filenames = [line.strip() for line in f]
            return (parse_pairs(name) for name in real_filenames)
    else:
        #filenames is just list names
        return (parse_pairs(name) for name in filenames)
def queue_write(res:"list of cell", out_name:str, replace:bool, filenames:list):
    '''
    write one by one. out_name to induce real outnames. filenames needed for replace mode.
    '''
    #use different writer
    if replace == True:
        for cell, out_name in zip(res, filenames):
            write_pairs(cell, out_name)
    else:
        if os.path.isdir(out_name):
            #use out_name as out directory, use standard appendix from cell.appendix
            for cell in res:
                write_pairs(cell, os.path.join(out_name, cell.name + cell.appendix))
        else:
            #use outname as namelist
            with open(out_name) as f:
                real_outnames = [line.strip() for line in f]
            for cell, real_outname in zip(res, real_outnames):
                write_pairs(cell, real_outname)          
def parse_i_pairs(filename:str)->"cell":
    '''
    read from hickit's imputed pairs
    '''
    #column names are in the last comment line.
    #head are stored in cell, handler may change it.
    #now return Data
    with gzip.open(filename,"rt") as f:
        head = []
        for line in f.readlines():
            if line[0] != "#":
                break
            head.append(line)
    pairs = pd.read_table(filename, header=None,comment="#")
    pairs.columns = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    pairs_data = Data("pairs","".join(head), pairs, ".pairs", filename)
    name, extend = divide_name(filename)
    return Cell(name, pairs_data)
def write_i_pairs(cell:Cell, out_name:str):
    #now use data
    '''
    write imputed dataframe to tab delimited zipped file
    reserve heads, no dataframe index and headers
    in_file == out_file will replace original file
    '''
    sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        f.write(cell.get_data("i_pairs").head)
        cell.get_data("i_pairs").content.to_csv(f, sep="\t", header=False, index=False, mode="a")