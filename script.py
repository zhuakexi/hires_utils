#script for dip-c pipeline
from argparse import ArgumentError
import os
from subprocess import run
def shear_name(filename:str,sep="."):
    basename = os.path.basename(filename)
    parts = basename.split(sep)
    return [part for part in parts if part != ""]
def is_file_type(filename:str, *key_words:str, sep:str="."):
    parts = shear_name(filename)
    return [part for part in parts if part in key_words] != []
def get_file_under(directory:str, file_sigs:list):
    files = os.listdir(directory)
    return [os.path.abspath(file) for file in files if is_file_type(file,*file_sigs)]
def cli(args):
    filenames, out_name, sub_dir_switch, paired_switch, by_cell_switch = \
        args.filenames, args.out_name, args.sub_dir_switch, args.paired_switch, args.by_cell_switch
    if paired_switch == True:
        pass
    else:
        if os.path.isdir(filenames[0]):
            #use raw/ name1 name2 ...
            #won't check name1/ name2/ ...
            directory, names = filenames[0], filenames[1:]
            [for name in names]  
        if len(filenames)%2 != 0:
            raise ArgumentError(None, "must give paired files")
def split_dna_rna(fastqs:"filename"):
    if len(fastqs) == 1:
        r1 = fastqs[0]
        print(r1)
    elif len(fastqs) == 2:
        r1, r2 = fastqs[0], fastqs[1]
        print(r1,r2)
def no_assigner():
    #run the easy way, ask for resources according to stringe stage
    pass
def assigner():
    pass
