#script for dip-c pipeline
from argparse import ArgumentError
import os
import sys
from subprocess import run
from functools import partial

import booter
from classes import Cell, Data
def shear_name(filename:str,sep="."):
    #shear filename, usefull when working on complex appendix
    basename = os.path.basename(filename)
    parts = basename.split(sep)
    return [part for part in parts if part != ""]
def splitext(filename:str,sep=".") -> tuple:
    #homemade os.path.splitext, original one can't handle complex extend name
    parts = shear_name(filename)
    if len(parts) == 0:
        raise ValueError("blank file name")
    elif len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], ".".join(parts[1:])
def is_file_type(filename:str, *key_words:str, sep:str="."):
    #get file type from key words
    parts = shear_name(filename)
    return [part for part in parts if part in key_words] != []
def get_file_under(directory:str, file_sigs:list)->list:
    #get all file of specific types under the directory
    directory = os.path.abspath(directory)
    files = os.listdir(directory) 
    return [os.path.join(directory, file) for file in files if is_file_type(file,*file_sigs)]
def is_sub_directory(directory, name):
    return os.path.normpath(name) in os.listdir(directory)
def cli(args):
    filenames, out_name, sub_dir_switch, paired_switch, by_cell_switch, num_thread = \
        args.filenames, args.out_name, args.sub_dir_switch, args.paired_switch, args.by_cell_switch, args.num_thread
    #get real fastq filenames and cell names
    nested_fastqs, cell_names = None, None
    if paired_switch == True:
        ##with 1 fastqs for each cell
        if all([os.path.isdir(filename) for filename in filenames]):
            ##cell1/ cell2/ ... 
            cell_names = [os.path.split(name)[1] for name in filenames]
            nested_fastqs = [get_file_under(filename, ["fq","fastq","fasta"]) for filename in filenames]  
        elif os.path.isdir(filenames[0]) and all([is_sub_directory(filenames[0],filename) for filename in filenames[1:]]):
            ##use raw/ name1 name2 ...
            directory, cell_names = filenames[0], filenames[1:]
            nested_fastqs = [get_file_under(os.path.join(directory, name),["fq","fastq","fasta"]) for name in cell_names]
        elif all([os.path.isfile(filename) for filename in filenames]):
            ##use cell1.fastq cell2.fastq ...
            nested_fastqs = [os.path.abspath(name) for name in filenames]
            cell_names = [splitext(name)[0] for name in filenames]
        else:
            sys.stderr.write("Input suggest:\n directory/ name1 name2 ...   OR\n 1.fastq 2.fastq ...   OR\n cell1/ cell2/ ...\n\n\n")
            raise ArgumentError(None, "wrong fastq")
    else:
        ##with 2 fastqs for each cell
        if all([os.path.isdir(filename) for filename in filenames]):    
            cell_names = [os.path.split(name)[1] for name in filenames]
            nested_fastqs = [get_file_under(filename,["fq","fastq","fasta"]) for filename in filenames]
        elif os.path.isdir(filenames[0]) and all([is_sub_directory(filenames[0],filename) for filename in filenames[1:]]):
            directory, cell_names = filenames[0], filenames[1:]
            nested_fastqs = [get_file_under(os.path.join(directory, name),["fq","fastq","fasta"]) for name in cell_names]  
        elif len(filenames)%2 == 1 and all([os.path.isfile(filename) for filename in filenames]):
            ##must give cell_names file. Don't want induction. (can give -f file1 file2)
            cell_names_file, fastq_names = filenames[0], filenames[1:]
            with open(cell_names_file) as f:
                cell_names = [line.strip() for line in f]
            nested_fastqs = list(zip(fastq_names[::2],fastq_names[1::2]))
        else:
            sys.stderr.write("Input suggest:\n raw/ name1 name2 ...   OR\n cell_names 1r1.fastq 1r2.fastq 2r1.fastq 2r2.fastq ...   OR\n cell1/ cell2/ ...\n\n\n")
            raise ArgumentError(None, "wrong fastq")
    #assert len(cell_names) == len(nested_fastqs)
    #print(cell_names)
    #print(nested_fastqs)
    # initiate cells for no_assigner, that is for split_dna_rna in fact
    fastq_data = [Data("fastq","","",".fastq.gz",file_list) for file_list in nested_fastqs]
    cells = [Cell(cell_name, data) for cell_name, data in zip(cell_names,fastq_data)]
    # call no_assinger in simple multiplex
    if len(cells) == 1:
        return no_assigner(cells[0], out_name, num_thread)
    else:
        '''
        workfunction = partial(no_assigner,out_name=out_name, num_thread=num_thread, )
        return booter.multi(workfunction, cells, len(cells)) #simple multiplex 
        '''
        pass

#sketch pipline functions here temperaly, 
#will create independent subcommand for step by step usage    
def split_dna_rna(cell, out_name, num_thread):
    aps = ["dna.R1.fq.gz","dna.R2.fq.gz","rna.R1.fq.gz","rna.R2.fq.gz"]
    subs = ["dna","dna","rna","rna"]
    real_out_names = [os.path.join(out_name, sub,".".join([cell.name, ap])) for sub, ap in zip(subs, aps)]
    '''
    print(
        "cutadapt -G 'XGGTTGAGGTAGTATTGCGCAATG;o=20' -j %d \
        --untrimmed-output %s \
        --untrimmed-paired-output %s \
        -o %s -p %s %s %s \
            " % (num_thread, *real_out_names, *cell.get_data["fastq"]))
    '''
    run(
        "cutadapt -G 'XGGTTGAGGTAGTATTGCGCAATG;o=20' -j %d \
        --untrimmed-output %s \
        --untrimmed-paired-output %s \
        -o %s -p %s %s %s \
            " % (num_thread, *real_out_names, *cell.get_data["fastq"]),
            shell=True)
def no_assigner(cell:Cell, out_name, num_thread):
    #run the easy way, ask for resources according to stringe stage
    #single runing, multiplexed by cli
    #not support by-cell-directory yet
    #print(cell_name, out_name, num_thread, fastq1, fastq2)
    split_dna_rna(cell, out_name, num_thread)
def assigner():
    pass
