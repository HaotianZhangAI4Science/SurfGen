#python gen.py --outdir example --ply_file ./example/3cl_pocket_8.0_res_1.5.ply
import os
import argparse
from glob import glob
from easydict import EasyDict
from rdkit import Chem
import torch
from copy import deepcopy
import shutil
import numpy as np
from tqdm.auto import tqdm
from utils.transforms import *
from utils.misc import load_config
from utils.reconstruct import *
from models.surfgen import SurfGen
from utils.chem import read_pkl, write_pkl
from utils.sample import get_init, get_next, logp_to_rank_prob, pdb_to_pocket_data
from utils.sample import STATUS_FINISHED, STATUS_RUNNING
import os.path as osp
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
import warnings

import pickle
def write_pkl(list,file):
    with open(file,'wb') as f:
        pickle.dump(list,f)
        print('pkl file saved at {}'.format(file))
def read_pkl(file):
    with open(file,'rb') as f:
        data = pickle.load(f)
    return data

parser = argparse.ArgumentParser()
parser.add_argument(
    '--config', type=str, default='./configs/sample.yml'
)

parser.add_argument(
    '--outdir', type=str, default='example',
)

parser.add_argument(
    '--device', type=str, default='cuda',
)

parser.add_argument(
    '--check_point',type=str,default='./ckpt/surfgen.pt',
    help='load the parameter'
)

parser.add_argument(
    '--ply_file', action='store',required=False,type=str,default='./example/3cl_pocket_8.0_res_1.5.ply',
    help='surface file specified for generation'
)


args = parser.parse_args()
config = load_config(args.config)

data = pdb_to_pocket_data(args.ply_file)

contrastive_sampler = ContrastiveSample()
ligand_featurizer = FeaturizeLigandAtom()
protein_featurizer = FeaturizeProteinAtom()
transform = Compose([
    RefineData(),
    LigandCountNeighbors(),
    Geodesic_builder(),
    ligand_featurizer,
    protein_featurizer
])


ckpt = torch.load(args.check_point, map_location=args.device)
model = SurfGen(
    ckpt['config'].model, 
    num_classes = contrastive_sampler.num_elements,
    num_bond_types = 3,
    protein_atom_feature_dim = protein_featurizer.feature_dim,
    ligand_atom_feature_dim = ligand_featurizer.feature_dim,
).to(args.device)
model.load_state_dict(ckpt['model'])
print('Num of parameters is {0:.4}M'.format(np.sum([p.numel() for p in model.parameters()]) /100000 ))

mask = LigandMaskAll()
composer = AtomComposer(5, ligand_featurizer.feature_dim, ckpt['config'].model.encoder.knn)
masking = Compose([
    mask, 
    composer
])
def transform_data(data, transform):
    assert data.protein_pos.size(0) > 0
    if transform is not None:
        data = transform(data)
    return data

data = transform(data)
data = transform_data(data, masking)
np.seterr(invalid='ignore') 
pool = EasyDict({
    'queue': [],
    'failed': [],
    'finished': [],
    'duplicate': [],
    'smiles': set(),
})

data = transform_data(deepcopy(data), masking)
init_data_list = get_init(data.to(args.device),   # sample the initial atoms
        model = model,
        transform=composer,
        threshold=config.sample.threshold
)
pool.queue = init_data_list
#rint('Start to generate novel molecules with 3D conformation located in the protein pocket!')
#print('The protein pocket is {}, init length is {}'.format(data.protein_filename, len(init_data_list)))
global_step = 0 
while len(pool.finished) < config.sample.num_samples:
    global_step += 1
    if global_step > config.sample.max_steps:
        break
    queue_size = len(pool.queue)
    # # sample candidate new mols from each parent mol
    queue_tmp = []
    for data in pool.queue:
        nexts = []
        data_next_list = get_next(
            data.to(args.device), 
            model = model,
            transform = composer,
            threshold = config.sample.threshold
        )

        for data_next in data_next_list:
            if data_next.status == STATUS_FINISHED:
                try:
                    rdmol = reconstruct_from_generated_with_edges(data_next)
                    data_next.rdmol = rdmol
                    mol = Chem.MolFromSmiles(Chem.MolToSmiles(rdmol))
                    smiles = Chem.MolToSmiles(mol)
                    data_next.smiles = smiles
                    if smiles in pool.smiles:
                        pool.duplicate.append(data_next)
                    elif '.' in smiles:
                        pool.failed.append(data_next)
                    else:   # Pass checks
                        print('Success: %s' % smiles)
                        pool.finished.append(data_next)
                        pool.smiles.add(smiles)
                except MolReconsError:
                    pool.failed.append(data_next)
            elif data_next.status == STATUS_RUNNING:
                nexts.append(data_next)

        queue_tmp += nexts
    prob = logp_to_rank_prob([p.average_logp[2:] for p in queue_tmp],)  # (logp_focal, logpdf_pos), logp_element, logp_hasatom, logp_bond
    n_tmp = len(queue_tmp)
    if n_tmp == 0:
        print('{}th has filures!'.format(global_step))
        break
    else:
        next_idx = np.random.choice(np.arange(n_tmp), p=prob, size=min(config.sample.beam_size, n_tmp), replace=False)
    pool.queue = [queue_tmp[idx] for idx in next_idx]

# save the generation results
task_name = args.ply_file.split('/')[-1][:-4]
task_dir = osp.join(args.outdir,task_name)
os.makedirs(task_dir,exist_ok=True)
sdf_file = os.path.join(task_dir,f'{task_name}_gen.sdf')
writer = Chem.SDWriter(sdf_file)
for j in range(len(pool['finished'])):
    writer.write(pool['finished'][j].rdmol)
writer.close()

SDF_dir = osp.join(task_dir,'SDF')
os.makedirs(SDF_dir, exist_ok=True)
for j in range(len(pool['finished'])):
    writer = Chem.SDWriter(SDF_dir+f'/{j}.sdf')
    writer.write(pool['finished'][j].rdmol)
    writer.close()

shutil.copy(args.ply_file,task_dir)

print('sucessfully generate at',task_dir)