from concurrent import futures
import sys
import time
import bisect
from functools import partial

import pandas as pd
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

def is_isolate(contact, sorted_contacts:"dataframe", up_dense, up_distance)->bool:
    '''
    #two_way_search version, calc in pos2-aixs strip
    #check if contact is isolated, work on one chromosome pair
    #contacts sorted by pos1 without index
    #if two contact pos1 Eu distance > up_distance then L-0.5 distance > up_distance
    '''
    proximity = 0
    index = bisect.bisect_right(sorted_contacts["pos1"], contact["pos1"])
    for _, con in sorted_contacts.iloc[index-1::-1].iterrows():
        if abs(con["pos1"]-contact["pos1"]) > up_distance:
            break
        if L_half(con,contact) < up_distance:
            proximity += 1 
    for _, con in sorted_contacts.iloc[index:].iterrows():
        if abs(con["pos1"]-contact["pos1"]) > up_distance:
            break
        if L_half(con,contact) < up_distance:
            proximity += 1 
    return proximity < up_dense + 1
def clean_contacts_in_pair(contacts:"dataframe", up_dense, up_distance)->"dataframe":
    '''
    #for multiplexing
    '''
    sorted_contacts = contacts.sort_values(by="pos1",axis=0,ignore_index=True)
    mask = contacts.apply(is_isolate, axis=1, sorted_contacts = sorted_contacts,
                          up_dense=up_dense, up_distance=up_distance)
    sys.stderr.write("%s %s %d --> %d" %(contacts.iloc[0]["chr1"], contacts.iloc[0]["chr2"], len(contacts),len(contacts[~mask])) )
    return contacts[~mask]
def clean_isolated_main(args):
    in_name, out_name, num_thread, up_dense, up_distance = \
        args.filenames[0], args.output_file, args.thread, args.dense, args.distance
    cell = pairs_parser(in_name)
    t0 = time.time()
    input_data = ( value for key, value in cell.groupby(["chr1","chr2"]) )
    working_func = partial(clean_contacts_in_pair, up_dense=up_dense, up_distance=up_distance)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input_data)
    cleaned = pd.concat(res,axis=0)
    sys.stderr.write("clean_isolated: finished in %.2fs"%(time.time()-t0))
    write_pairs(cleaned, in_name, out_name)
    return cleaned