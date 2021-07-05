# useful scripts to generate ref files

def gen_chr_alias(file_name:str="reference/name_alias.csv"):
    # alias for dip/hap chromosomes
    
    ## generate name list automatically
    base_ele = [str(i) for i in range(1,23)]
    base_ele.extend(["X", "Y"])

    hap_name = ["chr"+i for i in base_ele]

    dip_name0 = []
    for i in hap_name:
        dip_name0.append(i+"(mat)")
        dip_name0.append(i+"(pat)")

    dip_name1 = []
    for i in hap_name:
        dip_name1.append(i+"a")
        dip_name1.append(i+"b")
    # 1 -> chr1
    name_map0 = pd.DataFrame({"norm_name":hap_name, "alias":base_ele})
    # chr1a -> chr1(mat)
    name_map1 = pd.DataFrame({"norm_name":dip_name0, "alias":dip_name1})
    # x -> X
    name_map2 = pd.DataFrame({"norm_name":["chrX","chrY","chrX(mat)","chrX(pat)","chrY(mat)","chrY(pat)"],\
                              "alias":["chrx","chry", "chrxa","chrxb","chrya","chryb"]})
    # normal -> normal
    name_map3 = pd.DataFrame({"norm_name":dip_name0, "alias":dip_name0})
    res = pd.concat([name_map0, name_map1, name_map2, name_map3])
    res.to_csv(file_name, index=None)
    return res