# transform 3dg/xyz to mmcif
# @Date 210723

import re
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
def chrom_rm_suffix(chrom:pd.Series):
    """
    Remove suffix in chromosome name like (mat), (pat), a, b, A, B
    TODO: this is too specific, should be more general, or recieve a list of suffix to remove
    Input:
        chrom: pd.Series
    Output:
        pd.Series
    """
    if isinstance(chrom,str):
        return re.sub(r"\(?_?[mpatbAB]*\)?_?","",chrom)
    else:
        return chrom.str.replace(r"\(?_?[mpatbAB]*\)?_?","",regex=True)
def chrom_rm_prefix(chrom:str):
    """
    Remove prefix in chromosome name like chr, Chr, CHR
    Input:
        chrom: pd.Series
    Output:
        pd.Series
    """
    return chrom.str.replace(r"^[Cc][Hh][Rr]","",regex=True)
def threedg_to_cif(tdgPath:str,outputCifPath:str,factorBpath:str=None,dupref=True,maxGap:int=1000000,
    flavor="pymol",sample_name:str="hic",resolution:int=None):
    """
    Transform 3dg/xyz file to mmcif file
    Input:
        tdgPath: path to 3dg/xyz file
        outputCifPath: path to output mmcif file
        factorBpath: path to factorB file
        dupref: whether to draw diploid structure with haploid reference, like CpG frequency
        maxGap: max gap to connect two bins
        flavor: "pymol" or "chimera"
        sample_name: name of the sample, only used in chimera flavor
        resolution: resolution of the sample, only used in chimera flavor
    Output:
        None, write output to outputCifPath
    """
    if flavor == "chimera":
        assert resolution is not None, "Please specify resolution when using chimera flavor"
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
    if flavor == "chimera":
        if "factorB" not in positions.columns:
            positions["factorB"] = 0
        positions.columns = ["chrom","pos","x","y","z","factorB"]
        chimera_cif(sample_name, positions, resolution, "factorB", outputCifPath)
    elif flavor == "pymol":
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
    else:
        print("Unsupported flavor")
        sys.exit(1)

def chimera_cif(cellname, tdg, resolution, factor_b, outputpath=None):
    """
    Convert a DataFrame of 3D coordinates to a CIF file.
    Parameters:
        cellname : str
        tdg : pandas.DataFrame containing at least 'chrom', 'pos', 'x', 'y', 'z', 'CpG'
        resolution : int
        factor_b : str, column name for B factor in tdg
        outputpath : str, optional, default is './<cellname>.cif'
    Returns:
        None, writes output to specified outputpath
    """
    assert outputpath is not None, "Please specify an output path"

    file_head_name = f"data_{cellname}_res{int(resolution / 1000)}k"
    cif_str = f"#{file_head_name}\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n_entity_poly.nstd_linkage\n_entity_poly.nstd_monomer\n_entity_poly.pdbx_seq_one_letter_code\n_entity_poly.pdbx_seq_one_letter_code_can\n_entity_poly.pdbx_strand_id\n_entity_poly.pdbx_target_identifier\n"
    
    unique_chroms = tdg['chrom'].unique()
    def sort_chromosomes(chrom):
        num_part = ''.join(filter(str.isdigit, chrom)) 
        num = 100 if 'X' in chrom else 101 if 'Y' in chrom else int(num_part)
        return (chrom[-1], num)

    unique_chroms = sorted(unique_chroms, key=sort_chromosomes)
    tdg["chrom"] = pd.Categorical(tdg["chrom"], categories=unique_chroms)
    for i, chrom in enumerate(unique_chroms, start=1):
        cif_str += f"{i} 'Chromatin' no no ? ? ? ?\n"

    cif_str2 = "#\nloop_\n_entity.id\n_entity.type\n_entity.src_method\n_entity.pdbc_description\n_entity.formula_weight\n_entity.pdbx_number_of_molecules\n_entity.pdbx_ec\n_entity.pdbx_mutation\n_entity.pdbx_fragment\n_entity.details\n"
    for i, chrom in enumerate(unique_chroms, start=1):
        chrom_type = 'maternal' if 'mat' in chrom else 'paternal'
        chrom_num = ''.join(filter(str.isdigit, chrom))
        cif_str2 += f"{i} polymer man 'Chromosome {chrom_num} ({chrom_type})' ? ? ? ? ? ?\n"

    cif_str3 = "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\n"
    chrom_indices = {chrom: 0 for chrom in unique_chroms}
    for index, row in tdg.iterrows():
        chrom_index = chrom_indices[row['chrom']] + 1
        chrom_indices[row['chrom']] = chrom_index
        entity_id = unique_chroms.index(row['chrom']) + 1
        cif_str3 += f"ATOM {index+1} C CA . GLY {row['chrom']} {entity_id} {chrom_index} ? {row['x']} {row['y']} {row['z']} {row[factor_b]}\n"

    with open(outputpath, 'w') as f:
        f.write(cif_str)
        f.write(cif_str2)
        f.write(cif_str3)

    print(f"Done writing CIF file for {cellname}")
