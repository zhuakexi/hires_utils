import sys
import gzip
from .classes import Cell, Data
from .hires_io import parse_pairs
def write_data(cell:Cell, data_type:str, out_name:str):
    #now use data
    '''
    write dataframe to tab delimited zipped file
    reserve heads, no dataframe index and headers
    in_file == out_file will replace original file
    '''
    sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        f.write(cell.get_data(data_type).head)
        cell.get_data(data_type).content.to_csv(f, sep="\t", header=False, index=False, mode="a")
def cli(args):
    input_file, target, output_file = \
        args.input_file, args.target, args.output_file
    #print(input_file, target, output_file)
    main(input_file, target, output_file)
def get_target_pos(target:str)->list:
    target = target.split()
    result = []
    for name in ["chr_a", "cord_a", "chr_b", "cord_b"]:
        result.append(target.index(name))
    return result
def frame_transform(start:"DataFrame", target_pos:list) -> "DataFrame":
    #only support core 4 columns now, that is len(target_pos) = 4
    #assuming start is hires default pairs format
    #chr_a, cord_a, chr_b, cord_b positions are 1, 2, 3, 4(start from python 0), default names are chr1, pos1, chr2, pos2
    #will add format deduction later
    core = start[["chr1", "pos1", "chr2", "pos2"]]
    extra = start.drop(core, axis=1)
    #insert = lambda frame, loc, column : frame.insert(loc=loc, column=column.name, value=column)
    for column, position in zip(core.items(), target_pos):
        extra.insert(loc=position, column=column[0], value=column[1])
    return extra
def main(input_file, target, output_file):
    cell = parse_pairs(input_file)
    origin = cell.get_data("pairs").content
    pos = get_target_pos(target)
    result = Data(type="pairsa",head="", content=frame_transform(origin, pos), appendix="", file_name=output_file)
    cell.add_data(result)
    write_data(cell, "pairsa", output_file)
    


    