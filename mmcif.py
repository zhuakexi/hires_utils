# transform 3dg/xyz to mmcif
# using gemmi module
import gemmi
import pandas as pd

def cli(args):
    input_file, output_file = \
        args.input_file, args.output_file
    threedg_to_cif(input_file, output_file)
def threedg_to_cif(in_file:str, out_file:str):
    # read hickit 3dg file, write mmCIF file
    positions = pd.read_table(in_file,header=None,names="chr pos x y z".split()) # read in
    grouped = positions.groupby("chr") # split chromosomes
    
    model = gemmi.Model("Nuclear")
    for name, group in grouped:
        name = name.replace("(","_")
        name = name.replace(")","") # brackets not allowed in names
        chain = gemmi.Chain(name)
        for index, series in group.iterrows():
            # create atom for each bin
            atom = gemmi.Atom()
            pos = gemmi.Position(series["x"],series["y"],series["z"])
            atom.pos = pos
            atom.occ = 1
            # create residue for each atom
            residue = gemmi.Residue()
            residue.name = series["chr"] + "_" + str(series["pos"])
            residue.add_atom(atom,0)
            # add residue to chain
            chain.add_residue(residue,index)
        model.add_chain(chain)
    structure = gemmi.Structure()
    structure.name = in_file
    structure.add_model(model)
    structure.make_mmcif_document().write_file(out_file)