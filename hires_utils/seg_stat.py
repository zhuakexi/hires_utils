from .hires_io import parse_seg, gen_record, divide_name, print_records

def cli(args):
    filename, output, record_dir, sample_name = \
        args.filename[0], args.output, args.record_directory, args.sample_name
    hap1_phased, hap2_phased, biasedX_score, hap_score, yp = \
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
    try:
        yp = Y / (a+b+u) # y percent
        hap1_phased = a / (a+b+u)
        hap2_phased = b / (a+b+u)
        biasedX_score = abs(Xa - Xb)/(Xa + Xb)
        hap_score = abs(a - b)/(a + b)
    except ZeroDivisionError:
        # fully capture the Error
        if a + b == 0:
            return 0.0, 0.0, -1, 1.0, yp
        if Xa + Xb == 0:
            return hap1_phased, hap2_phased, -1, hap_score, yp 
    return hap1_phased, hap2_phased, biasedX_score, hap_score, yp
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