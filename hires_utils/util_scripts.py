# useful scripts to generate ref files

def gen_chr_alias(file_name:str="name_alias.csv"):
    # alias for dip/hap chromosomes
    
    ## generate name list automatically
    base_ele = [str(i) for i in range(1,23)]

    autosomes = ["chr" + i for i in base_ele]
    
    # diploid genome: 22x2 + X(mat) + X(pat) + Y(pat)
    
    # normal chrom name set in hires-utils
    dip_genome0 = []
    for i in autosomes:
        dip_genome0.append(i+"(mat)")
        dip_genome0.append(i+"(pat)")
    dip_genome0.extend(["chrX(mat)", "chrX(pat)", "chrY(pat)"])

    # a variant diploid chromosome name set using by hickit build 3dg
    # in this name set, a:(pat), b:(mat)
    dip_genome1 = []
    for i in autosomes:
        dip_genome1.append(i+"b")
        dip_genome1.append(i+"a")
    dip_genome1.extend(["chrXb", "chrXa", "Y"])
    
    # normal haploid chrom name set in hires-utils
    hap_genome0 = autosomes + ["chrX", "chrY"]
    
    # a variant used in GRCh37/hg19 reference sequence
    hap_genome1 = base_ele + ["X", "Y"]
    
    # hap_genome1 -> hap_genome0 aka 1 -> chr1
    name_map0 = pd.DataFrame({"norm_name":hap_genome0, "alias":hap_genome1})
    # dip_genome1 -> dip_genome0 aka chr1a -> chr1(mat)
    name_map1 = pd.DataFrame({"norm_name":dip_genome0, "alias":dip_genome1})
    # exception1(first capital) aka x -> X
    name_map2 = pd.DataFrame({"norm_name":["chrX","chrY","chrX(mat)","chrX(pat)","chrY(pat)"],\
                              "alias":["chrx","chry", "chrxb","chrxa", "chrya"]})
    # normal -> normal
    name_map3 = pd.DataFrame({"norm_name":dip_genome0, "alias":dip_genome0})
    name_map4 = pd.DataFrame({"norm_name":hap_genome0, "alias":hap_genome0})
    
    # generate table
    res = pd.concat([name_map0, name_map1, name_map2, name_map3, name_map4])
    res.to_csv(file_name, index=None)
    return res