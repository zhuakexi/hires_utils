import argparse
from clean_leg import clean_leg_main
from clean_splicing import clean_splicing_main
from align import align_main
from clean_isolated import clean_isolated_main
def cli():
    parser = argparse.ArgumentParser(prog="hires", description="Functions for hires pipline")
    subcommands = parser.add_subparsers(title="These are sub-commands",metavar="command")
    #clean_legs sub command
    clean_legs = subcommands.add_parser(
                            "clean_leg",
                            help="clean promiscuous legs that contacts with multiple legs")
    clean_legs.set_defaults(handle=clean_leg_main)
    clean_legs.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            nargs=1)
    clean_legs.add_argument(
                            "-t", "--thread",
                            type=int,
                            dest="thread",
                            action="store",
                            default=20,
                            help="set thread number")
    clean_legs_out = clean_legs.add_mutually_exclusive_group(required=True)
    clean_legs_out.add_argument("-s", "--replace", 
                            dest="replace_switch", 
                            action="store_true", 
                            default=False,
                            help="do clean in-place and replace input file")
    clean_legs_out.add_argument("-o", "--output", 
                            dest="out_name", action="store",
                            help="set output file name")

    #clean_splicing sub command
    clean_splicing = subcommands.add_parser(
                            "clean_splicing", 
                            help="clean exon splicing from mRNA in contact file")
    clean_splicing.set_defaults(handle=clean_splicing_main)
    clean_splicing.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1)       
    clean_splicing.add_argument(
                            "-r", "--reference", 
                            dest="index_file_name", 
                            help="exon index file, use 'build' sub-command to build from scratch.", 
                            default="bin_10k_FULL_index")
    clean_splicing.add_argument(
                            "-b", "--binsize",
                            dest="binsize",
                            type=int,
                            default=10000)  
    clean_splicing_out = clean_splicing.add_mutually_exclusive_group(required=True)
    clean_splicing_out.add_argument(
                            "-s", "--replace", 
                            dest="replace_switch", 
                            help="do clean in-place", 
                            action="store_true", 
                            default=False)
    clean_splicing_out.add_argument(
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
    clean_isolated = subcommands.add_parser(
                            "clean_isolated",
                            help="remove isolated contacts according to L-0.5 distance"
    )
    clean_isolated.set_defaults(handle=clean_isolated_main)
    clean_isolated.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1
    )      
    clean_isolated.add_argument(
                            "-o","--output",
                            dest="output_file",
                            type=str,
                            required=True
    )
    clean_isolated.add_argument(
                            "-t","--thread",
                            dest="thread",
                            type=int,
                            default=23
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