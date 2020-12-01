import pandas as pd
from classes import Cell, Data
from clean_isolated import clean_isolated
from hires_io import parse_i_pairs, write_pairs, write_i_pairs
def cli(args):
    file_name, out_name1, out_name2, num_thread, up_dense, up_distance = \
        args.input_file, args.output_file1, args.output_file2, int(args.num_thread), int(args.dense), int(args.distance)
    cell = parse_i_pairs(file_name)
    cello1, cello2 = sep_clean(cell, num_thread, up_dense, up_distance)
    write_pairs(cello1, out_name1)
    write_i_pairs(cello2, out_name2)

# --------- working module --------
hap_word = {"1":"(mat)", "0":"(pat)"} #1 for maternal, 0 for paternal
def add_ap(data:pd.DataFrame, row_picker:pd.Series, col_picker:list, ap:str):
    # add appendix for any data subset, defined by raw and col picker
    # do inplace
    # in pandas, assign value must be done in single step, don't do chained indexing
    ap = [ap for i in range(0, len(data.loc[row_picker, :]))] # try to figure out a more elegant way
    data.loc[row_picker, col_picker].astype("string")
    data.loc[row_picker, col_picker] = data.loc[row_picker, col_picker].str.cat(ap)
def rm_hap(data:pd.DataFrame):
    data["chr1"] = data["chr1"].astype("string").str.extract(r'(chr[\dXY]+)')
    data["chr2"] = data["chr2"].astype("string").str.extract(r'(chr[\dXY]+)')
    return data
def sep_clean(cell: Cell, num_thread, up_dense, up_distance) -> Cell:
    # sep pairs haplotype and do clean isolated again
    # generate .hap.pairs for 2D analysis as well as .pairs for hickit 3d building
    name, data = cell.name, cell.get_data("pairs")
    head, frame = data.head, data.content
    new_frame = pd.DataFrame()
    # blank DataFrame can concat any DataFrame
    for code in "00 01 10 11".split():
        row_picker = frame["phase_prob" + code] >= 0.75
        frame.loc[row_picker, "phase0"] = code[0]
        frame.loc[row_picker, "phase1"] = code[1]
        add_ap(frame, row_picker, "chr1", hap_word[code[0]])
        add_ap(frame, row_picker, "chr2", hap_word[code[1]])
        new_frame = pd.concat([new_frame, frame.loc[row_picker, :]])
    cell.get_data("pairs").content = new_frame # cell seperated
    cell_c = clean_isolated(cell, num_thread, up_dense, up_distance) # cell cleaned
    cell_cf = Cell(name, Data(dtype="pairs", head=head, content=cell_c.get_data("pairs").content["readID chr1 pos1 chr2 pos2 strand1 strand2".split()])) # filter out redundant columns 
    cell_cm = Cell(name, Data(dtype="i_pairs", head=head, content=rm_hap(cell_c.get_data("pairs").content)) ) # cell mingled
    
    return cell_cf, cell_cm