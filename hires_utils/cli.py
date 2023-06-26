import argparse
from . import clean_leg
from . import clean_splicing
from . import clean_isolated
from . import pairsa
from . import sep_clean
from . import chrom_rmsd
from . import clean3
from . import mmcif
from . import clean_leg
from . import gcount
from . import seg_stat
from . import mend_umi
from . import mend_cellbarcode
def cli():
    parser = argparse.ArgumentParser(prog="hires", description="Functions for hires pipline")
    subcommands = parser.add_subparsers(title="These are sub-commands",metavar="command")
#--------- clean_leg sub command ---
    clean_leg_arg = subcommands.add_parser(
                            "clean_leg",
                            help="clean promiscuous legs that contacts with multiple legs")
    clean_leg_arg.set_defaults(handle=clean_leg.cli)
    clean_leg_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            nargs=1)
    clean_leg_arg.add_argument(
                            "-t", "--thread",
                            type=int,
                            dest="thread",
                            action="store",
                            default=4,
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
    clean_leg_arg.add_argument("-o", "--output", 
                            dest="out_name", action="store",
                            metavar="OUTPUT_FILE",
                            required=True,
                            help="set output file name")
 # ---------#clean_splicing sub command ---
    clean_splicing_arg = subcommands.add_parser(
                            "clean_splicing", 
                            help="clean exon splicing from mRNA in contact file")
    clean_splicing_arg.set_defaults(handle=clean_splicing.cli)
    clean_splicing_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1)     
    clean_splicing_arg.add_argument(
                            "-r", "--reference", 
                            dest="gtf_filename",
                            type = str,
                            action="store", 
                            help="annotation gtf file", 
                            required=True)
    clean_splicing_arg.add_argument(
                            "-o", "--output", 
                            dest="out_name",
                            metavar="OUTPUT_FILE",
                            required=True, 
                            help="output file name", 
                            action="store")
    clean_splicing_arg.add_argument(
                            "-t", "--thread",
                            type=int,
                            dest="num_thread",
                            action="store",
                            default=4,
                            help="set thread number")
#--------- chrom_rmsd subcommand ------
    rmsd_arg = subcommands.add_parser(
                            "rmsd",
                            help="calculate rmsd between .3dg replicates, pick good ones"
    )
    rmsd_arg.set_defaults(handle=chrom_rmsd.cli)
    rmsd_arg.add_argument(
                            dest="filenames",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs="*"
    )
    rmsd_arg.add_argument(
                            "-o", "--output",
                            dest="result_log",
                            metavar="RESULT_LOG_FILE",
                            help="a log file to aggregate stat info, useful in workflow",
                            type=str,
                            default=None
    )
    rmsd_arg.add_argument(
                            "-rd", "--record_directory",
                            dest="record_dir",
                            metavar="DIR",
                            help="directory to store key-value json for stat",
                            type=str,
                            default=None
    )
    rmsd_arg.add_argument(
        "-sa","--sample_name",
        dest="sample_name",
        metavar="SAMPLE",
        type=str,
        help="[for pipeline] sample name/ID, used in record",
        default=None
    )
    rmsd_arg.add_argument(
        "-at","--attributes",
        dest="attr",
        metavar="ATTR",
        type=str,
        help="[for pipeline] data key, usually binsize of 3d particle, used in record",
        default=None
    )
#--------- clean3 subcommand ------
    clean3_arg = subcommands.add_parser(
                            "clean3",
                            help="clean 3dg particles poorly supported by contacts"
    )
    clean3_arg.set_defaults(handle=clean3.cli)
    clean3_arg.add_argument(
                            "-i","--input",
                            dest="filename",
                            metavar="STRUCTURE_FILE",
                            help=".3dg/.xyz format structure file to clean",
                            type=str,
                            required=True
    )
    clean3_arg.add_argument(
                            "-r", "--reference",
                            dest="ref_filename",
                            metavar="PAIRS",
                            help=".pairs format contact file",
                            type=str,
                            required=True
    )
    clean3_arg.add_argument(
                            "-o", "--output",
                            dest="output",
                            metavar="CLEANED",
                            help="file name of the output cleaned structure file",
                            type=str,
                            required=True
    )
    clean3_arg.add_argument(
                            "-q", "--quantile",
                            dest="quantile",
                            metavar="QUANTILE",
                            help="quantile of particles to remove",
                            type=int,
                            default=0.06
    )
    clean3_arg.add_argument(
                            "-d", "--distance",
                            dest="distance",
                            metavar="DISTANCE",
                            help="max distance (bp) from a contact leg to a 3D genome particle",
                            type=int,
                            default=500_000
    )
#--------- mmcif subcommand ------
    mmcif_arg = subcommands.add_parser(
                            "mmcif",
                            help="transform 3dg/xyz file to mmcif"
    )
    mmcif_arg.set_defaults(handle=mmcif.cli)
    mmcif_arg.add_argument(
                            "-i", "--input",
                            dest="input_file",
                            metavar="INPUT",
                            help="input xyz/3dg file",
                            type=str,
                            required=True
    )
    mmcif_arg.add_argument(
                            "-o", "--output",
                            dest="output_file",
                            metavar="OUTPUT",
                            help="name of the output mmcif file",
                            type=str,
                            required=True
    )
    mmcif_arg.add_argument(
                            "-b", "--factorB",
                            dest="factorBpath",
                            help="factorB for color your mmcif file",
                            type=str,
                            default=None,
    )
    mmcif_arg.add_argument(
                            "-g", "--maxGap",
                            dest="maxGap",
                            help="max gap for not to have a bond between tow atoms",
                            type=int,
                            default=1000000,
    )
#--------- gcount subcommand ------
    gcount_arg = subcommands.add_parser(
        "gcount",
        help="count number of records \
            in sequencing-related files,\
             eg: fastq, pairs"
    )
    gcount_arg.set_defaults(handle=gcount.cli)
    gcount_arg.add_argument(
        dest="filename",
        metavar="INPUT_FILE",
        help="input filename",
        nargs=1
    )
    gcount_arg.add_argument(
        "-f","--format",
        dest="fformat",
        metavar="FORMAT",
        help="input file format: \
            pe_fastq, se_fastq, merged_fastq, pairs",
        type=str,
        default="pe_fastq"
    )
    gcount_arg.add_argument(
        "-rd","--record_directory",
        dest="record_directory",
        metavar="DIR",
        type=str,
        help="record directory to store key:value results",
        default=None
    )
    gcount_arg.add_argument(
        "-sa","--sample_name",
        dest="sample_name",
        metavar="SAMPLE",
        type=str,
        help="[for pipeline] sample name/ID, used in record",
        default=None
    )
    gcount_arg.add_argument(
        "-at","--attributes",
        dest="attributes",
        metavar="ATTR",
        type=str,
        help="[for pipeline] data key, used in record",
        default=None
    )
#---------seg_stat subcommand ------
    seg_stat_arg = subcommands.add_parser(
        "seg_stat",
        help="summary .seg file"
    )
    seg_stat_arg.set_defaults(handle=seg_stat.cli)
    seg_stat_arg.add_argument(
        dest="filename",
        metavar="INPUT_FILE",
        help="input filename, dip-c .seg file",
        nargs = 1
    )
    seg_stat_arg.add_argument(
        "-o", "--output",
        dest="output",
        metavar="LOGFILE",
        help="output file"
    )
    seg_stat_arg.add_argument(
        "-rd","--record_directory",
        dest="record_directory",
        metavar="DIR",
        type=str,
        help="[for pipeline] record directory to store key:value results",
        default=None
    )
    seg_stat_arg.add_argument(
        "-sa","--sample_name",
        dest="sample_name",
        metavar="SAMPLE",
        type=str,
        help="[for pipeline] sample name/ID, used in record",
        default=None
    )
    seg_stat_arg.add_argument(
        "-dump", "--dump",
        dest = "dump",
        default= True,
        action="store_true",
        help="whether to dump per-chromosome counting, store in additional *dump* dir if rd enabled"
    )
    
#--------- clean_isolate subcommand ------
    clean_isolated_arg = subcommands.add_parser(
                            "clean_isolated",
                            help="remove isolated contacts according to L-0.5 distance"
    )
    clean_isolated_arg.set_defaults(handle=clean_isolated.cli)
    clean_isolated_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1
    )
    clean_isolated_arg.add_argument(
                            "-t","--thread",
                            dest="thread",
                            type=int,
                            help="set thread number",
                            default=4
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
    clean_isolated_arg.add_argument(
                            "-o","--output",
                            dest="output_file",
                            action="store",
                            metavar="OUTPUT_FILE",
                            required=True,
                            help = "output file name",
                            type=str
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
# --------- mend_umi subcommand ------
    mend_umi_arg = subcommands.add_parser(
        "mend_umi",
        help = "adding UB tag with umi stored in readID"
    )
    mend_umi_arg.set_defaults(handle=mend_umi.cli)
    mend_umi_arg.add_argument(
        "--input","-i",
        dest="input_file",
        help="input file path",
        action = "store",
        required=True
    )
    mend_umi_arg.add_argument(
        "--output","-o",
        dest="output_file",
        help="output file path",
        action="store",
        required=True
    )
# --------- mend_cellbarcode subcommand ------
    mend_cellbarcode_arg = subcommands.add_parser(
        "mend_cellbarcode",
        help = "adding CB tag with user-defined cell barcode"
    )
    mend_cellbarcode_arg.set_defaults(handle=mend_cellbarcode.cli)
    mend_cellbarcode_arg.add_argument(
        "--input","-i",
        dest="input_file",
        help="input file path",
        action = "store",
        required=True
    )
    mend_cellbarcode_arg.add_argument(
        "--cellbarcode","-c",
        dest="cellbarcode",
        help="cell barcode",
        action="store",
        required=True
    )
    mend_cellbarcode_arg.add_argument(
        "--output","-o",
        dest="output_file",
        help="output file path",
        action="store",
        required=True
    )
# --------- sep_clean subcommand ---------
    sep_clean_arg = subcommands.add_parser(
                            "sep_clean",
                            help = "seperate homologous chromosome(add (mat)/(pat) for chra chrb colomns), clean isolated contacts again. \
                                generate one more hickit compatible output file.\
                                    works with hickit imputed pairs file"
    )
    sep_clean_arg.set_defaults(handle=sep_clean.cli)
    sep_clean_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1
    )
    sep_clean_arg.add_argument(
                            "-n", "--num_thread",
                            dest="num_thread",
                            help="number of thread use",
                            default="4"
    )
    sep_clean_arg.add_argument(
                            "-o1", "--output1",
                            dest="output_file1",
                            help="output file path for .dip.pairs (the chra(mat) format) file",
                            action="store",
                            required=True
    )
    sep_clean_arg.add_argument(
                            "-o2", "--output2",
                            dest="output_file2",
                            help="output file path for .pairs (the hickit -b default format) file",
                            required=True
    )
    sep_clean_arg.add_argument(
                            "-m","--dense",
                            dest="dense",
                            type=int,
                            help="number of contacts in proximity",
                            default=5
    )
    sep_clean_arg.add_argument(
                            "-d","--distance",
                            dest="distance",
                            type=int,
                            help="check contacts in what L-0.5 range",
                            default=10000000
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