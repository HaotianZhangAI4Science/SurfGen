from espsim.electrostatics import GetMolProps, GetEspSim
from rdkit import Chem
from rdkit.Chem import AllChem

if __name__ == 'main':
    ori_mol = Chem.MolFormMol2File('ori_ligand.mol2')
    gen_mol = Chem.MolFormMol2File('gen_ligand.mol2')
    shape_sim = 1-AllChem.ShapeTanimotoDist(ori_mol, gen_mol)
    els_sim = 1-GetEspSim(ori_mol,gen_mol)