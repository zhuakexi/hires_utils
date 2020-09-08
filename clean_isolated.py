from concurrent import futures
import sys
import time
import bisect
from functools import partial

import pandas as pd
import numpy as np

from hires_io import pairs_parser
from hires_io import write_pairs
from batch import batch


def L_half(contact1, contact2):
    '''
    caculate L 0.5 norm of contact pair
    in general pandas groupby give chr1 < chr2
    light version, assume two contacts has same chromosome order
    '''
    d_x, d_y = abs(contact1.pos1-contact2.pos1), abs(contact1.pos2-contact2.pos2)
    return (np.math.sqrt(d_x) + np.math.sqrt(d_y))**2
def is_isolate(contact, sorted_contacts:"dataframe", up_dense, up_distance)->bool:
    '''
    check if contact is isolated, work on one chromosome pair
    contacts sorted by pos1 without index
    if two contact pos1 Eu distance > up_distance then L-0.5 distance > up_distance
    light version
    '''
    proximity = 0
    index = bisect.bisect_right(sorted_contacts["pos1"], contact["pos1"])
    for con in sorted_contacts.iloc[index-1::-1].itertuples(index=False):
        if abs(con.pos1-contact.pos1) > up_distance or proximity >= up_dense+1:
            break
        if L_half(con,contact) <= up_distance:
            proximity += 1 
    for con in sorted_contacts.iloc[index:].itertuples(index=False):
        if abs(con.pos1-contact.pos1) > up_distance or proximity >= up_dense+1:
            break
        if L_half(con,contact) <= up_distance:
            proximity += 1 
    return proximity < up_dense + 1
def clean_contacts_in_pair(contacts:"dataframe", up_dense, up_distance)->"dataframe":
    '''
    #for multiplexing
    '''
    sorted_contacts = contacts.sort_values(by="pos1",axis=0,ignore_index=True)
    mask = contacts.apply(is_isolate, axis=1, sorted_contacts = sorted_contacts,
                          up_dense=up_dense, up_distance=up_distance)
    sys.stderr.write("(%s, %s): %d --> %d\n" %(contacts.iloc[0]["chr1"], contacts.iloc[0]["chr2"], len(contacts),len(contacts[~mask])) )
    return contacts[~mask]
def cli(args):
    filenames, out_name, num_thread, up_dense, up_distance, batch_switch, replace = \
        args.filenames, args.output_file, args.thread, args.dense, args.distance, args.batch_switch, args.replace_switch
    #case1: multi mode. multiple in files begin a loop 
    if len(filenames) > 1:
        for cell_name in filenames:
            if replace == True:
                #--replace will work in multi file input
                the_out_name = cell_name
            else:
                #--out_name will be used as name appendix: xx.appendix.pairs.gz 
                the_out_name = cell_name.split(".")
                the_out_name.insert(1,out_name)
                the_out_name = ".".join(the_out_name)
            clean_isolated_main(cell_name,the_out_name,  num_thread, up_dense, up_distance)
        return 0
    #case2: in batch mode. call batch function to do loop
    if batch_switch == True:
        working_func = partial(clean_isolated_main, num_thread=num_thread, up_dense=up_dense, up_distance=up_distance)
        batch(working_func, filenames, out_name, replace)
    #case3: in single mode. run once.
    cell_name = filenames[0]
    if replace == True:
        the_out_name = cell_name
    else:
        the_out_name = out_name    
    clean_isolated_main(cell_name,the_out_name,  num_thread, up_dense, up_distance)
    return 0
def clean_isolated_main(in_name, out_name, num_thread, up_dense, up_distance):
    cell = pairs_parser(in_name)
    t0 = time.time()
    input_data = ( value for key, value in cell.groupby(["chr1","chr2"]) )
    working_func = partial(clean_contacts_in_pair, up_dense=up_dense, up_distance=up_distance)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input_data)
    cleaned = pd.concat(res,axis=0)
    print("clean_isolated: %d contacts removed in %s" % (len(cell)-len(cleaned),in_name))
    sys.stderr.write("clean_isolated: finished in %.2fs\n"%(time.time()-t0))
    write_pairs(cleaned, in_name, out_name)
    return cleaned