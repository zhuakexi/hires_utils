from argparse import ArgumentError
import os

from assist import cli_assist
from hires_io import parse_pairs, parse_3dg
def cli(args):
    structure, pairs, out_filename, clean_quantile, max_clean_distance = \
        args.filename, args.ref_filename, args.output, args.quantile, args.distance
    pass
def clean3(structure, pairs, clean_quantile, max_clean_distance):
    pass