from hires_io import parse_pairs, write_pairs
from classes import Cell

def single(handle:"one cell argument", filename:str, outname:str)->int:
    #call handle once
    cell = parse_pairs(filename)
    res = handle(cell)
    write_pairs(res, outname)