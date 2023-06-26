import pysam
import time
import sys
def cli(args):
    input_file, cellbarcode, output_file = \
        args.input_file, args.cellbarcode, args.output_file
    mend_cellbarcode(input_file, cellbarcode, output_file)
def mend_cellbarcode(filei, barcode, fileo):
    """
    Write provided cell barcode to a seperate CB tag.
    This convert plate-style library (with one file per cell) to dropseq-style library.
    Typically work with uBAM.
    """
    start = time.time()
    with pysam.AlignmentFile(filei,"rb",threads=2, check_sq=False) as infile:
        with pysam.AlignmentFile(fileo,"wb",template=infile,threads=8) as outfile:
            for line in infile:
                tags = line.get_tags()
                tags.extend([("CB",barcode)])
                line.set_tags(tags)
                outfile.write(line)
    print("mend_cellbarcode finished, running time:",(time.time()-start)/60)