import pysam
import time
import sys
def cli(args):
    input_file, output_file = \
        args.input_file, args.output_file
    mend_umi(input_file, output_file)
def mend_umi(filei, fileo):
    start = time.time()
    with pysam.AlignmentFile(filei,"rb",threads=2) as infile:
        with pysam.AlignmentFile(fileo,"wb",template=infile,threads=8) as outfile:
            for line in infile:
                _, umi = line.qname.split("_")
                tags = line.get_tags()
                tags.extend([("UB",umi)])
                line.set_tags(tags)
                outfile.write(line)
    print("mend_umi finished, running time:",(time.time()-start)/60)