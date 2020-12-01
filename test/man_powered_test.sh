
#python3 hires.py script -h
#python3 hires.py script -o out -pr raw/Cell1 raw/Cell2
#python3 hires.py script -o out -pr raw/ Cell1 Cell2
#python3 hires.py script -o out -pr raw/cella.fq.gz raw/cellb.fq.gz
#python3 hires.py script -o out  raw/Cell1 raw/Cell2
#python3 hires.py script -o out  raw/ Cell1 Cell2
#python3 hires.py script -o out  raw/cell_list raw/cella.fq.gz raw/cella.r2.fq.gz raw/cellb.fq.gz raw/cellb.r2.fq.gz
#python3 hires.py script -o out  raw/ Cell1
#python3 hires.py script -o out  raw/cell_list raw/cella.fq.gz raw/cella.r2.fq.gz raw/cellb.fq.gz raw/cellb.r2.fq.gz

# suppose you have test.pairs.gz in raw/ folder
## test for argument parsing
#python3 hires.py -h
#python3 hires.py pairsa -h
#python3 hires.py pairsa --input raw/test.pairs.gz --target "x chr_a cord_a x chr_b cord_b" --output out/test.validPairs.gz

# test for sep_clean
python3 hires.py -h
python3 hires.py sep_hap -h
python3 hires.py sep_hap -n 4 -i test.impute.pairs.gz -o1 out/test.hap.pairs.gz -o2 out/for_hickit.pairs.gz