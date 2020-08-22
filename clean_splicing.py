import sys
import pickle
import random
import time
import sys
import gzip
import pandas as pd

#chromsome contacts and reference only
regular_chromsome_names = ["chr" + str(i) for i in range(1,23)]
regular_chromsome_names.extend(["chrX","chrY"])
def pairs_parser(cell_name:str)->"dataframe":
    '''
    read from 4DN's .pairs format
    '''
    t0 = time.time()
    column_names = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1".split()
    pairs = pd.read_table(cell_name, header=None,skiprows=27)
    pairs.columns = column_names
    sys.stderr.write("pairs_parser parsing time: %.2fs\n"%(time.time()-t0))
    return pairs
def bin_parser(file_name, regular="off"):
    '''
    read bin.bed file
    return dictionary of bins by chromsome
    '''
    time_begin = time.time()
    bins = pd.read_table(file_name, header=None)
    print("bin_parser parsing time: " + str(time.time() - time_begin) + "\n")
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
def clean_splicing_main(args):
    '''
    clean contacts from splicing
    '''
    BINSIZE, index_name, cell_name, out_name, replace = \
        args.binsize, args.index_file_name, args.filenames[0], args.out_name, args.replace_switch
    with gzip.open(cell_name, "rt") as f:
        #get file head
        head = [next(f) for i in range(0,27)]
        head = "".join(head)
    #get real data
    cell = pairs_parser(cell_name)

    # load directly from pickled bin_index
    with open(index_name,"rb") as f:
        bin_index = pickle.load(f)
    # do searching
    print("block_search: working...")
    hit, cleaned = block_search(bin_index, BINSIZE, cell)

    if replace == True:
        out_name = cell_name
    with gzip.open(out_name,"wt") as f:
        f.write(head)
    cleaned.to_csv(out_name, sep="\t", \
        header=False, index=False, mode="a")
    sys.stderr.write("block_search total questionable contacts: %d \n" %len(hit) )