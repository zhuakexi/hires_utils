import argparse
import clean_leg
import clean_splicing
import clean_isolated
import script
import pairsa
from align import align_main

def cli():
    parser = argparse.ArgumentParser(prog="hires", description="Functions for hires pipline")
    subcommands = parser.add_subparsers(title="These are sub-commands",metavar="command")
#--------- clean_leg sub command ---
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
    clean_leg_arg_out = clean_leg_arg.add_mutually_exclusive_group(required=True)
    clean_leg_arg_out.add_argument("-s", "--replace", 
                            dest="replace_switch", 
                            action="store_true", 
                            default=False,
                            help="do clean in-place and replace input file")
    clean_leg_arg_out.add_argument("-o", "--output", 
                            dest="out_name", action="store",
                            help="set output file name, or output file appendix for multiple file")
    ##parsing different strategy
    clean_leg_arg_strategy = clean_leg_arg.add_mutually_exclusive_group()
    clean_leg_arg_strategy.add_argument(
                            "-p", "--parallel",
                            dest = "parallel_switch",
                            help="do in parallel mode, maximum throghput, give filelist or directory",
                            action="store_true",
                            default=False
    )                                           
    clean_leg_arg_strategy.add_argument(
        "-b", "--batch",
        dest = "batch_switch",
        help="do in batch mode, long and stable at night, give filelist or directory",
        action="store_true",
        default=False
    ) 
 # ---------#clean_splicing sub command ---
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
                            type = str,
                            action="store", 
                            help="exon index file, use 'build' sub-command to build from scratch.", 
                            default="bin_10k_FULL_index")
    clean_splicing_arg.add_argument(
                            "-bin", "--binsize",
                            dest="binsize",
                            type=int,
                            default=10000)
    ##parsing different strategy
    clean_splicing_arg_strategy = clean_splicing_arg.add_mutually_exclusive_group()
    clean_splicing_arg_strategy.add_argument(
                            "-p", "--parallel",
                            dest = "parallel_switch",
                            help="do in parallel mode, maximum throghput, give filelist or directory",
                            action="store_true",
                            default=False
    )                                           
    clean_splicing_arg_strategy.add_argument(
        "-b", "--batch",
        dest = "batch_switch",
        help="do in batch mode, long and stable at night, give filelist or directory",
        action="store_true",
        default=False
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
    
#--------- align subcommand ------
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
                            "-o","--output_dir",
                            dest="output_dir",
                            type=str,
                            help="directory to store aligned 3dg file and rmsd info file",
                            required=True
    )
    align.add_argument(
        "-gd", "--good_dir",
        dest="good_dir",
        help="output directory for good(low rmsd) structure.",
        required=True
    )
    align.add_argument(
        "-bd", "--bad_dir",
        dest="bad_dir",
        help="output directory for bad(causing high rmsd) structure.",
        required=True
    )
#--------- clean_isolate subcommand ------
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
    ##parsing different strategy
    clean_isolated_arg_strategy = clean_isolated_arg.add_mutually_exclusive_group()
    clean_isolated_arg_strategy.add_argument(
        "-p", "--parallel",
        dest = "parallel_switch",
        help="do in parallel mode, maximum throghput, give filelist or directory",
        action="store_true",
        default=False
    )                                           
    clean_isolated_arg_strategy.add_argument(
        "-b", "--batch",
        dest = "batch_switch",
        help="do in batch mode, long and stable at night, give filelist or directory",
        action="store_true",
        default=False   
    )
# --------- pairsa subcommand ------
    pairsa_arg = subcommands.add_parser(
        "pairsa",
        help = "transform .pairs to user-defined pairs-like file formats. e.g. validPairs."
    )
    pairsa_arg.set_defaults(handle=pairsa.cli)
    pairsa_arg.add_argument(
        "--input",
        dest="input_file",
        help="input file path",
        required=True
    )
    pairsa_arg.add_argument(
        "--target",
        dest="target",
        help="a blank seperated string giving out target format\n\
            must have 'chr_a' 'chr_b' 'cord_a' 'cord_b'\n\
                the position of the column name in this 'target string' will be the position of that column\n\
                e.g. 'id chr_a cord_a strand_a chr_b cord_b strand_b restriction_distance'\n\
                    defined the validPairs format from HiCPro and Juicer.",
        action="store",
        required=True
    )
    pairsa_arg.add_argument(
        "--output",
        dest="output_file",
        help="output file path",
        action="store",
        required=True
    )
# --------- script subcommand ------
    script_arg = subcommands.add_parser(
        "script",
        help = "handle from fastq to mmCIF, with multiple cell"
    )
    script_arg.set_defaults(handle=script.cli)
    script_arg.add_argument(
        dest="filenames",
        help="doubled fastq files, or paired single fastq files( check the --paired_in option), or a directory.",
        nargs="+"
    )
    script_arg.add_argument(
        "-o","--output",
        dest="out_name",
        help="output directory, with all mid-files in( check -sub or -cell option for more tidier organization).",
        metavar="DIR/",
        action="store",
        required=True
    )
    ## directory architectual options
    dir_arch = script_arg.add_mutually_exclusive_group()
    dir_arch.add_argument(
        "-sub","--by_type",
        dest="sub_dir_switch",
        help="using subdirectories like seg/ pairs/ 3dg/ to store output, create if not prepared.",
        action="store_true",
        default=False
    )
    dir_arch.add_argument(
        "-cell","--by_cell",
        dest="by_cell_switch",
        help="using cell name as subdirecotries.",
        action="store_true",
        default=False
    )
    script_arg.add_argument(
        "-pr", "--paired_in",
        dest="paired_switch",
        help="use single fastq file with both reads and mate reads in.",
        action="store_true",
        default=False
    )
    script_arg.add_argument(
        "-t", "--num_thread",
        dest="num_thread",
        help="thread used on one cell(Notion:not the total core number allocated).",
        action="store",
        default=4
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