import numpy as np
from hires_io import pairs_parser
from hires_io import write_pairs

def L_half(contact1, contact2):
    '''
    caculate L 0.5 norm of contact pair
    in general pandas groupby give chr1 < chr2
    '''
    if ( contact1["chr1"] == contact2["chr1"] ) and ( contact1["chr2"] == contact2["chr2"] ):
        d_x, d_y = abs(contact1["pos1"]-contact2["pos1"]), abs(contact1["pos2"]-contact2["pos2"])
        #print((np.math.sqrt(d_x) + np.math.sqrt(d_y))**2)
        return (np.math.sqrt(d_x) + np.math.sqrt(d_y))**2
    elif ( contact1["chr1"] == contact2["chr2"] ) and ( contact1["chr2"] == contact2["chr1"] ):
        d_x, d_y = abs(contact1["pos1"]-contact2["pos2"]), abs(contact1["pos2"]-contact2["pos1"])
        #print((np.math.sqrt(d_x) + np.math.sqrt(d_y))**2)
        return (np.math.sqrt(d_x) + np.math.sqrt(d_y))**2
    else:
        #print(np.inf,contact1)
        return np.inf

def is_isolated(contact:"line",cell:"dict of dataframe")->bool:
    key = tuple( sorted([contact["chr1"],contact["chr2"]]) )
    count = cell[key].apply(L_half, axis=1, contact2=contact) < 10000000
    print("proxy num:",count.astype(np.int8).sum())
    #print(count.astype(np.int8).sum())
    return count.astype(np.int8).sum()< (5+1)
def clean_contacts(contacts:"dataframe",cell:"dict of dataframe")->"two dataframe":
    '''
    for multiplexing
    '''
    mask = contacts.apply(is_isolated, axis=1, cell=cell)
    return contacts[~mask], contacts[mask]
def clean_isolated_main(args):
    in_name, out_name = args.filenames[0], args.output_file
    cell = pairs_parser(in_name)
    paired_chromosomes = { key:value for key, value in cell.groupby(["chr1","chr2"]) }
    
    '''
    blocks = np.array_split(cell, 1000000, axis=0)
    cleaned, hit = clean_contacts(blocks[0], paired_chromosomes)
    '''
    input_data = cell.iloc[0:3]
    cleaned, hit = clean_contacts(input_data, paired_chromosomes)
    return cleaned