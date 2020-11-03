import pandas as pd
from classes import Cell, Data




# --------- working module --------
hap_word = {"1":"(mat)", "0":"(pat)"} #1 for maternal, 0 for paternal
def add_ap(data:"DataFrame", row_picker:"boolean series", col_picker:"name of column", ap:"appendix"):
    # assign value must be done in single step, don't do chained indexing
    ap = [ap for i in range(0, len(data.loc[row_picker, :]))] # try to figure out a more elegant way
    data.loc[row_picker, col_picker].astype("string")
    data.loc[row_picker, col_picker] = data.loc[row_picker, col_picker].str.cat(ap)
def main(cell: Cell) -> Cell:
    data = cell.get_data("i_pairs")
    frame = ["content"]
    new_frame = pd.DataFrame()
    # blank DataFrame can concat any DataFrame
    for code in "00 01 10 11".split():
        row_picker = frame["phase_prob" + code] >= 0.75
        frame.loc[row_picker, "phase0"] = code[0]
        frame.loc[row_picker, "phase1"] = code[1]
        add_ap(frame, row_picker, "chr1", hap_word[code[0]])
        add_ap(frame, row_picker, "chr2", hap_word[code[1]])
        new_frame = pd.concat([new_frame, frame.loc[row_picker, :]])
        #print(frame.loc[row_picker, :] )
    #new_data = Data("i_pairs", data.head, new_frame, "hap.pairs.gz", ) 
    cell.get_data("i_pairs").content = new_frame
    return cell