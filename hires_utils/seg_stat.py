from .hires_io import parse_seg, gen_record, divide_name, print_records

def cli(args):
    filename, output, record_dir, sample_name = \
        args.filename[0], args.output, args.record_directory, args.sample_name
    hap1_phased, hap2_phased, biasedX_score, hap_score = \
        seg_values(filename)
    assigned = judge(biasedX_score, hap_score)
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
    Xu = Xa = Xb = Y = 0
    u = a = b = 0
    _, segs = parse_seg(filename)
    for seg in segs:
        #print(seg)
        attrs = seg.split("!")
        if attrs[0] == "chrX":
            if attrs[4] == "0":
                Xa += 1
            elif attrs[4] == "1":
                Xb += 1
            else:
                Xu += 1
        elif attrs[0] == "chrY":
            Y += 1
        # autosome 
        elif attrs[4] == "0":
            a += 1
        elif attrs[4] == "1":
            b += 1
        else:
            u += 1
    hap1_phased = a / (a+b+u)
    hap2_phased = b / (a+b+u)
    biasedX_score = abs(Xa - Xb)/(Xa + Xb)
    hap_score = abs(a - b)/(a + b)
    return hap1_phased, hap2_phased, biasedX_score, hap_score
def judge(biasedX_score:float, hap_score:float)->str:
    # singleX_score [0,1] 0:female 1:male
    # hap_score [0,1] 0:dip 1:hap
    # singleX_score threshold: <0.25 or >0.8
    # hap_score threshold: <0.2 or >0.9
    tx1 = 0.25
    tx2 = 0.8
    th1 = 0.2
    th2 = 0.9
    if hap_score > th2:
        # haploid
        if biasedX_score > tx2:
            # in hap, biasedX because has X
            return "hapfem"
        if biasedX_score < tx1:
            # in hap, no bias because doesn't have X
            return "hapmal"
    if hap_score < th1:
        # diploid
        if biasedX_score > tx2:
            # in dip, biasedX becasue has one X,
            # that is, has Y
            return "dipmal"
        if biasedX_score < tx1:
            # in dip, no biasedX because has 2 X
            return "dipfem"
    return "unassigned"