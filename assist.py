from argparse import ArgumentError
from classes import Cell
import os

def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 
def get_cell_name(filename):
    cell_name, _ = divide_name(filename)
    return cell_name
def cli_assist(parallel_switch, batch_switch, filenames, out_name):
    #parse normal arguments for each function's cli
    result = {}
    if parallel_switch == True:
        #not support now
       result["strategy"] = "parallel" 
    if batch_switch == True:
        result["strategy"] = "batch"
        if len(filenames) != 1:
            raise ArgumentError(None, "Please give one name list file as input in batch mode.")
        else:
            filename = filenames
            if os.path.isdir(filename):
                #use it as a directory
                result["real_in_names"] = [name for name in os.listdir(filename)]
            else:
                with open(filename) as f:
                    result["real_in_names"] = [line.strip() for line in f]
    if not parallel_switch and not batch_switch:
        result["strategy"] = "multi"
        if len(filenames) > 1:
            result["real_in_names"] = filenames
        else:
            filename = filenames[0]
            if os.path.isdir(filename):
                #use it as a directory
                result["real_in_names"] = [name for name in os.listdir(filename)]
            else:
                #single filename
                result["real_in_names"] = [filename]
    result["for_init"] = [Cell(get_cell_name(name)) for name in result["real_in_names"]]
    return result