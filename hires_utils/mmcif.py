# transform 3dg/xyz to mmcif
# @Date 210723

import sys

import pandas as pd
from .hires_io import parse_3dg, parse_ref

def cli(args):
    input_file, output_file, factorBpath, maxGap= \
        args.input_file, args.output_file, args.factorBpath, args.maxGap
    threedg_to_cif(input_file, output_file, factorBpath, maxGap)

class Bin:
    def __init__(self,id:int,chromosome:str,position:int,x_coord:float,y_coord:float,z_coord:float,factorB:float=None):
        self.id = id
        self.chromosome = chromosome
        self.position = position
        self.nextBin = None
        self.previousBin = None
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord
        self.factorB = factorB
   
    def outputAtomLine(self,atomNum:int):
        """
        """
        return "\t".join(["HETATM",".",str(atomNum),self.chromosome,"003",
                          "1",str(self.position),str(self.x_coord),str(self.y_coord),str(self.z_coord),
                          str(self.factorB if self.factorB else "."),self.chromosome])
    def outputBondLine(self,bondNum:int,maxGap:int):
        """
        return bond of this bin and it's next when they are connected
        """
        if (self.chromosome == self.nextBin.chromosome and abs(self.nextBin.position-self.position) <= maxGap):
            return "\t".join([str(bondNum),"covale",self.chromosome,"003",
                          "1",str(self.position),self.nextBin.chromosome,"003",
                          "1",str(self.nextBin.position)])
        else: return None
def chrom_rm_suffix(chrom:str):
    """
    Remove suffix in chromosome name like (mat), (pat), a, b, A, B
    TODO: this is too specific, should be more general, or recieve a list of suffix to remove
    Input:
        chrom: pd.Series
    Output:
        pd.Series
    """
    return chrom.str.replace(r"\(?[mpatbAB]*\)?","",regex=True)       
def threedg_to_cif(tdgPath:str,outputCifPath:str,factorBpath:str=None,dupref=True,maxGap:int=1000000):
    """
    Transform 3dg/xyz file to mmcif file
    Input:
        tdgPath: path to 3dg/xyz file
        outputCifPath: path to output mmcif file
        factorBpath: path to factorB file
        dupref: whether to draw diploid structure with haploid reference, like CpG frequency
        maxGap: max gap to connect two bins
    """
    positions = parse_3dg(tdgPath).reset_index()
    if(factorBpath!=None):
        factorB = parse_ref(factorBpath,value_name="factorB").reset_index()
        if dupref:
            positions = positions.assign(new_chr=chrom_rm_suffix(positions["chr"]))
            factorB.rename(columns={"chr":"new_chr"},inplace=True)
            positions = pd.merge(
                positions,
                factorB,
                on = ["new_chr","pos"],
                how = "left"
                ).drop(columns=["new_chr"])
        else:
            positions = pd.merge(
                positions,
                factorB,
                on=["chr","pos"],
                how="left"
                )
    grouped = positions.groupby("chr") # split chromosomes
    
    binList = []
    binNum = 0
    for chr_name, coord_per_chrom in grouped:
        for index, series in coord_per_chrom.iterrows():
            # create atom for each bin
            if (factorBpath):
                currentBin = Bin(index,chr_name,series['pos'],series["x"],series["y"],series["z"],series["factorB"])
            else : currentBin = Bin(index,chr_name,series['pos'],series["x"],series["y"],series["z"])
            if(binNum !=0):
                binList[-1].nextBin = currentBin
            binNum += 1
            # if b factor is specificed.
            binList.append(currentBin)
    #print(binList)
    binList[-1].nextBin = binList[0]
    
    with open(outputCifPath,"w") as output_cif:
        output_cif.write("data_"+output_cif.name.replace(".cif","")+"\n")
        output_cif.write("#\n")
        output_cif.write("#\n")
        #write all connnection
        output_cif.write("""loop_\n_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
""")
        bondIndex = 1
        for bin in binList:
            #print(bondIndex)
            temp = bin.outputBondLine(bondIndex,maxGap)
            if temp:
                output_cif.write(temp + "\n")
                bondIndex += 1

        #write all atoms
        output_cif.write("""##
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.id
_atom_site.label_asym_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.label_atom_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.auth_asym_id
""")
        atomIndex = 1
        for bin in binList:
            output_cif.write(bin.outputAtomLine(atomIndex)+"\n")
            atomIndex += 1


if __name__ == '__main__':
    import sys
    sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
    from io import StringIO
    
    import pandas as pd
    from hic_basic.plot.render import clip_b_pymol
    # test diploid
    clip_b_pymol(
        "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg",
        "/share/home/ychi/software/dip-c/color/hg19.cpg.20k.hom.txt",
        "/share/home/ychi/dev/hires_utils/out/GMO1001.clean.20k.4.cpg.dip.png",
        tmpdir = "/share/home/ychi/dev/hires_utils/out/",
        dupref=False
        )
    # test haploid
    clip_b_pymol(
        "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg",
        "/share/home/ychi/software/dip-c/color/hg19.cpg.20k.txt",
        "/share/home/ychi/dev/hires_utils/out/GMO1001.clean.20k.4.cpg.hap.png",
        tmpdir = "/share/home/ychi/dev/hires_utils/out/",
        )
