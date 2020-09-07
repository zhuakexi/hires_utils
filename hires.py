import argparse
import clean_leg
import clean_splicing
import clean_isolated
from align import align_main

def cli():
    parser = argparse.ArgumentParser(prog="hires", description="Functions for hires pipline")
    subcommands = parser.add_subparsers(title="These are sub-commands",metavar="command")
    #clean_legs sub command
    clean_leg_arg = subcommands.add_parser(
                            "clean_leg",
                            help="clean promiscuous legs that contacts with multiple legs")
    clean_leg_arg.set_defaults(handle=clean_leg.cli)
    clean_leg_arg.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            nargs="*")
    clean_leg_arg.add_argument(
                            "-t", "--thread",
                            type=int,
                            dest="thread",
                            action="store",
                            default=20,
                            help="set thread number")
    clean_leg_arg.add_argument(
                            "-d","--distance",
                            dest="max_distance",
                            metavar="MAX_DISTANCE",
                            type=int,
                            action="store",
                            default=1000,
                            help="max distance to calculate adjacent legs"
    )
    clean_leg_arg.add_argument(
                            "-n","--count",
                            metavar="MAX_COUNT",
                            dest="max_count",
                            type=int,
                            action="store",
                            default=10,
                            help="number threshold of adjacent legs"
    )     
    clean_leg_arg.add_argument(
                            "-b","--batch",
                            dest="batch_switch",
                            action="store_true",
                            default=False,
                            help="batch mode, -o for output directory, FILENAME for cell list file"
    )              
    clean_leg_arg_out = clean_leg_arg.add_mutually_exclusive_group(required=True)
    clean_leg_arg_out.add_argument("-s", "--replace", 
                            dest="replace_switch", 
                            action="store_true", 
                            default=False,
                            help="do clean in-place and replace input file")
    clean_leg_arg_out.add_argument("-o", "--output", 
                            dest="out_name", action="store",
                            help="set output file name, or output file appendix for multiple file")
    #clean_splicing sub command
    clean_splicing_arg = subcommands.add_parser(
                            "clean_splicing", 
                            help="clean exon splicing from mRNA in contact file")
    clean_splicing_arg.set_defaults(handle=clean_splicing.cli)
    clean_splicing_arg.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs="*")     
    clean_splicing_arg.add_argument(
                            "-r", "--reference", 
                            dest="index_file_name", 
                            help="exon index file, use 'build' sub-command to build from scratch.", 
                            default="bin_10k_FULL_index")
    clean_splicing_arg.add_argument(
                            "-bin", "--binsize",
                            dest="binsize",
                            type=int,
                            default=10000)
    clean_splicing_arg.add_argument(
                            "-b", "--batch",
                            dest="batch_switch",
                            action="store_true",
                            default=False,
                            help="batch mode, give more files or one list file"
    )                           
    clean_splicing_arg_out = clean_splicing_arg.add_mutually_exclusive_group(required=True)
    clean_splicing_arg_out.add_argument(
                            "-s", "--replace", 
                            dest="replace_switch", 
                            help="do clean in-place", 
                            action="store_true", 
                            default=False)
    clean_splicing_arg_out.add_argument(
                            "-o", "--output", 
                            dest="out_name", 
                            help="output file name", 
                            action="store")
    #align subcommand
    align = subcommands.add_parser(
                            "align",
                            help="caculate rmsd between .3dg replicates"
    )
    align.set_defaults(handle=align_main)
    align.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs="*"
    )
    align.add_argument(
                            "-o","--output_prefix",
                            dest="output_prefix",
                            type=str,
                            default=""
    )
    #clean_isolate subcommand
    clean_isolated_arg = subcommands.add_parser(
                            "clean_isolated",
                            help="remove isolated contacts according to L-0.5 distance"
    )
    clean_isolated_arg.set_defaults(handle=clean_isolated.cli)
    clean_isolated_arg.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs="*"
    )
    clean_isolated_arg.add_argument(
                            "-b","--batch",
                            dest="batch_switch",
                            action="store_true",
                            help="batch mode, give more files or one list file",
                            default=False
    )      

    clean_isolated_arg.add_argument(
                            "-t","--thread",
                            dest="thread",
                            type=int,
                            default=23
    )     
    clean_isolated_arg.add_argument(
                            "-m","--dense",
                            dest="dense",
                            type=int,
                            help="number of contacts in proximity",
                            default=5)
    clean_isolated_arg.add_argument(
                           "-d","--distance",
                            dest="distance",
                            type=int,
                            help="check contacts in what L-0.5 range",
                            default=10000000) 
    clean_isolated_arg_out = clean_isolated_arg.add_mutually_exclusive_group(required=True)
    clean_isolated_arg_out.add_argument(
                            "-o","--output",
                            dest="output_file",
                            type=str
    )
    clean_isolated_arg_out.add_argument(
                            "-r","--replace",
                            dest="replace_switch",
                            action="store_true",
                            default=False
    )   

    args = parser.parse_args()
    #print(args.replace_switch)
    #print(args.out_name)
    #print(args.filenames)
    if hasattr(args, "handle"):
        args.handle(args)
    else:
        parser.print_help()
if __name__ == "__main__":
    cli()