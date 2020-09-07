import argparse
import sys
def batch(handle, filenames, out_name, replace):
    if len(filenames) > 1:
        raise argparse.ArgumentError(None, "in batch mode only one name list file is permitted.")
    #in batch mode, filenames is parsed a name list file, by line
    with open(filenames[0]) as f:
        filenames = [line.strip() for line in f.readlines()]
    if out_name[-1] != "/":
        #not directory but out name list file
        with open(out_name) as f:
            outlist = [line.strip() for line in f.readlines()]
        if len(outlist) != len(filenames):
            raise argparse.ArgumentError(None, "using out name list now: outlist must be same length with inlist.")
        else:
            #outlist and inlist good, begin loop 
            for i, cell_name in enumerate(filenames):
                the_out_name = outlist[i]
                sys.stderr.write("batch clean splicing: working on %s\n"%cell_name)
                handle(cell_name, the_out_name)
            return 0
    #case2:not outlist but out directory/replace, begin loop
    for cell_name in filenames:
        if replace == True:
            the_out_name = cell_name
        else:
            #using --outname as out directory
            the_out_name = out_name + cell_name.split("/")[-1]
        handle(cell_name, the_out_name)
    return 0
