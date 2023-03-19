import argparse
import os
import numpy as np
import subprocess
import pymesh
import tempfile, shutil
#import Bio.PDB
from Bio.PDB import PDBParser, PDBIO, Select 
from Bio.PDB import NeighborSearch, Selection
from rdkit import Chem
from scipy.spatial import distance, KDTree
from IPython.utils import io
from joblib import Parallel, delayed
import sys
sys.path.append('../utils/masif')
from compute_normal import compute_normal
from computeAPBS import computeAPBS
from computeCharges import computeCharges, assignChargesToNewMesh
from computeHydrophobicity import computeHydrophobicity
from computeMSMS import computeMSMS
from fixmesh import fix_mesh
from save_ply import save_ply


parser = argparse.ArgumentParser()
parser.add_argument(
    '--pdb_file', action='store',required=False,type=str,default='./1z6e_protein.pdb',
    help='protein file'
)

parser.add_argument(
    '--lig_file', action='store',required=False,type=str,default='./1z6e_ligand.mol2',
    help='protein file'
)

parser.add_argument(
    '--check_software', action='store',required=True,type=str,default='no',
    help='have you ever installed the dependent softwares and manually specify their locations?'
)


args = parser.parse_args()