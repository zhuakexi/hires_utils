from argparse import ArgumentError
import os

from assist import cli_assist
from hires_io import parse_pairs, parse_3dg
def cli(args):
    filenames, references, out_name, num_thread, parallel_switch, batch_switch = \
        args.filenames, args.references, args.out_name,args.num_thread, args.parallel_switch, args.batch_switch
    ticket = cli_assist(parallel_switch, batch_switch, filenames, out_name)
    cells = ticket["for_init"]
    in_names = ticket["real_in_names"]
    if len(references) != len(filenames):
        raise ArgumentError(None, "one reference .pairs file for each 3dg")
    else:
        if len(filenames) > 1:
            real_ref_names = references
        else:
            reference = references[0]
            if os.path.isdir(reference):
                #use it as a directory
                real_ref_names = [name for name in os.listdir(reference)]
            else:
                if (batch_switch or parallel_switch) == True:
                    with open(reference) as f:
                        real_ref_names = [line.strip() for line in f]
                else:                    
                    real_ref_names = references  
    assert len(real_ref_names) == len(cells)          
    for cell, in_name, ref_name in zip(cells, in_names, real_ref_names):
        cell.add_data(parse_3dg(in_name))
        cell.add_data(parse_pairs(ref_name)) 
           
def clean3(cell):
    pass