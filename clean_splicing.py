import sys
import pickle
import random
import time
import gzip
import argparse
from functools import partial

import pandas as pd


import hires_io
import booter
from classes import Cell

#chromsome contacts and reference only
regular_chromsome_names = ["chr" + str(i) for i in range(1,23)]
regular_chromsome_names.extend(["chrX","chrY"])

def bin_parser(file_name, regular="off"):
    '''
    read bin.bed file
    return dictionary of bins by chromsome
    '''
    time_begin = time.time()
    bins = pd.read_table(file_name, header=None)
    sys.stderr.write("bin_parser parsing time: " + str(time.time() - time_begin) + "\n")
    grouped = bins.groupby(0) #group by chromsome names
    if regular == "on":
        return {key:value.values for key,value in grouped if key in regular_chromsome_names}
        #return [bin for bin in bins.values if bin[0] in regular_chromsome_names]
    return {key:value.values for key,value in grouped}
def filt_in_exon(locus, exons):
    '''
    filter out exons envolope the locus
    '''
    return [exon for exon in exons if exon[0] <= locus <= exon[1]]
def key_to_index(locus, bin_size):
    '''
    in bin_index
    '''
    return locus//bin_size
def in_exon(contact:"line", bin_index:dict, binsize:int)->bool:
    '''
    whether both leg of a contact are in exons from same gene
    '''
    if contact["chr1"] != contact["chr2"]:
        return False
    left_index, right_index = key_to_index(contact["pos1"], binsize), key_to_index(contact["pos2"], binsize)
    left_hit_exons = filt_in_exon(contact["pos1"], bin_index[contact["chr1"]][left_index])
    if left_hit_exons != []:
        right_hit_exons = filt_in_exon(contact["pos2"], bin_index[contact["chr2"]][right_index])
        left_hit_genes = set([exon[2] for exon in left_hit_exons])
        right_hit_genes = set([exon[2] for exon in right_hit_exons])
        return left_hit_genes.intersection(right_hit_genes) != set()
    else:
        return False
def block_search(bin_index:"dict of list", binsize:int, cell:"dataframe")->"data_frame":
    t0 = time.time()
    #vectorize using .pairs, target form
    mask = cell.apply(in_exon, axis=1, bin_index=bin_index, binsize=binsize)
    hit_contacts = cell[mask]
    cleaned_contacts = cell[~mask]
    sys.stderr.write("block_search searching time: %.2fs\n" % (time.time()-t0))
    return hit_contacts, cleaned_contacts

def cli(args)->int:
    BINSIZE, index_name, filenames, out_name, replace, parallel_switch, batch_switch = \
        args.binsize, args.index_file_name, args.filenames, args.out_name, args.replace_switch, args.parallel_switch, args.batch_switch
    working_function = partial(clean_splicing, index_name=index_name, BINSIZE=BINSIZE)
    if parallel_switch == True:
        return booter.parallel(working_function, filenames, out_name, replace)
    if batch_switch == True:
        return booter.batch(clean_splicing, filenames, out_name, replace)
    if not parallel_switch and not batch_switch:
        if len(filenames) > 1:
            return booter.multi(working_function, filenames, out_name, replace)
        else:
            return booter.single(clean_splicing, filenames, out_name, replace)
#def clean_splicing_main(cell_name, out_name, index_name, BINSIZE):
def clean_splicing(cell:Cell, index_name:str, BINSIZE:int)->Cell:
    '''
    clean contacts from splicing
    '''
    # load directly from pickled bin_index
    with open(index_name,"rb") as f:
        bin_index = pickle.load(f)
    # do searching
    hit, cleaned = block_search(bin_index, BINSIZE, cell.data)
    print("clean_splicing: %d contacts removed in %s\n" %(len(hit), cell.name) )
    cell.data = cleaned
    cell.appendix = ".pairs.gz" #keep old one
    return cell