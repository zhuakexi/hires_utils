#Calculate Dip-C format of CpG Count
#Usage: snakemake -j10 -s generateCpG.smk

rule all:
    input:
        expand("mm10.CpG.{resolution}.txt",resolution=[10000,20000,50000,100000,200000,500000,1000000])

rule generateCpG:
    input:
        genome = "/share/Data/public/ref_genome/mouse_ref/M23/raw_data/genome.fa",
        index = "/share/Data/public/ref_genome/mouse_ref/M23/raw_data/genome.fa.fai",
    output:
        cpg = "mm10.CpG.{resolution}.txt"
    shell:"""
        set +u;source activate;conda activate py3;set -u
        bedtools getfasta -fi {input.genome} -bed <( bedtools makewindows -g <( cat {input.index} | cut -f1,2 ) -w {wildcards.resolution} | bedtools shift -i stdin -s `expr {wildcards.resolution} / 2` -g {input.index} ) -bedOut |awk '{{print $1"\t"$2"\t"$3"\t"gsub("CG","xx",$4)/length($4)}}'| awk '{{if ($4!=0) print $0}}' | awk '{{print $1"\t"($2+$3)/2"\t"$4}}'> {output.cpg}
        set +u;conda deactivate; set -u
    """

