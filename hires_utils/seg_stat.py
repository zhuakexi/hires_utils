import pandas as pd
from .hires_io import parse_seg, gen_record, divide_name, print_records

def cli(args):
    filename, output, record_dir, sample_name = \
        args.filename[0], args.output, args.record_directory, args.sample_name
    hap1_phased, hap2_phased, biasedX_score, hap_score, yp, xp, cp_counts = \
        seg_values(filename)
    assigned = judge(hap_score, yp)
    # generate records
    if sample_name == None:
        sample_name,_ = divide_name(filename)
    records = {
        sample_name:
            {
                "hap1_phased":hap1_phased,
                "hap2_phased":hap2_phased,
                "biasedX_score":biasedX_score,
                "hap_score":hap_score,
                "ypercent":yp,
                "xpercent":xp,
                "cell_state":assigned
            }
    }
    if output == None:
        # print to stdout
        print_records(records)
    else:
        # print to log file
        with open(output,"wt") as f:
            print_records(records, f)
    if record_dir != None:
        gen_record(records, record_dir)
def seg_values(filename:str)->tuple:
    comments, legs = parse_seg(filename)
    df = pd.DataFrame(
        (leg.split("!") for leg in legs),
        columns=["chrom","genome_start","genome_end","strand","phasing","a","b"]
    )
    cp_count = df.value_counts(["chrom","phasing"])
    defaults = pd.DataFrame([(line.split()[1],phase) for line in comments for phase in [".","0","1"]], columns = ["chrom","phasing"])
    defaults.set_index(["chrom","phasing"], inplace=True)
    res = pd.concat([defaults, cp_count], axis=1, join="outer").iloc[:,0]
    res.fillna(0,inplace=True)
    
    u, a, b = res.drop(["chrX","chrY"],level=0).groupby("phasing").sum()
    Xu, Xa, Xb = res["chrX"]
    Y = res["chrY"].sum()
    X = res["chrX"].sum()
    
    try:
        xp = X / (a+b+u) # x percent
        yp = Y / (a+b+u) # y percent
        hap1_phased = a / (a+b+u)
        hap2_phased = b / (a+b+u)
    except ZeroDivisionError:
        yp, hap1_phased, hap2_phased= \
            -1, 0.0, 0.0
    try:
        hap_score = abs(a - b)/(a + b)
    except ZeroDivisionError:
        hap_score = -1
    try:
        biasedX_score = abs(Xa - Xb)/(Xa + Xb)
    except ZeroDivisionError:
        biasedX_score = -1
    return hap1_phased, hap2_phased, biasedX_score, hap_score, yp, xp, res.to_dict()
def judge(hap_score:float, yp:float)->str:
    # hap_score [0,1] 0:dip 1:hap
    # try best to avoid unassigned
    # Tx = 0.75 # 0~0.75 that is, 0.5~0.75 for maternal_X_percent
    Th = 0.5 # for most either <0.2 or >0.9
    Ty = 0.0005
    if hap_score > Th:
        # haploid, using ypercent
        if yp < Ty:
            # in hap, biasedX because has X
            return "hapfem"
        if yp > Ty:
            # in hap, no bias because doesn't have X
            return "hapmal"
    if hap_score < Th:
        # diploid, using both ypercent and biasedX
        # ypercent is the de facto judger
        if yp > Ty:
            # in dip, biasedX becasue has one X,
            # that is, has Y
            return "dipmal"
        if yp < Ty:
            # in dip, no biasedX because has 2 X
            return "dipfem"
    return "unassigned"