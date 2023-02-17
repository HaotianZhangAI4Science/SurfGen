import sys
sys.path.append('/home/haotian/molecules_confs/Protein_test/Res2mol/script')
from universe import *
import numpy as np
from rdkit import DataStructs

def compute_sim(ref, gen, source='mol'):
    if source =='mol':
        fp_refmol = Chem.RDKFingerprint(ref)
        fp_genmol = Chem.RDKFingerprint(gen)
        sim = DataStructs.TanimotoSimilarity(fp_refmol, fp_genmol)
    elif source == 'fp':
        sim = DataStructs.TanimotoSimilarity(ref, gen)
    else:
        raise NotImplementedError('Error: you must choose the mol or fp to compute the similariy')
    return sim

def compute_sims(refs, gens, source='mol'):
    '''
    source could be the mol/fp
    '''
    num_refs = len(refs)
    num_gens = len(gens)
    sims = []
    for i in range(num_refs):
        for j in range(num_gens):
            sim = compute_sim(refs[i], gens[j], source=source)
            sims.append(sim)
    return np.array(sims).reshape(num_refs,num_gens)

def compute_fps(mols):
    fps = [Chem.RDKFingerprint(i) for i in mols]
    return fps
def compute_smi_from_fps(i):
    smis = []
    for j in range(len(train_fps)):
        smi = DataStructs.TanimotoSimilarity(fps[i], train_fps[j])
        smis.append(smi)
    return smis
    
#compute the fps
mols = {}
fps = {}
mols['GraphBP'] = read_pkl('./fps/mols/GraphBP_val_allmols.pkl')
mols['pkt2mol'] = read_pkl('./fps/mols/pkt2mol_val_allmols.pkl')
mols['ResGen'] = read_pkl('./fps/mols/resgen_val_allmols.pkl')
methods = list(mols.keys())
for method in methods:
    fps[method] = compute_fps(mols[method])
    write_pkl(fps[method], f'./fps/mols/{method}_valmols_fps.pkl')

# Get the train split and val split
from torch.utils.data import Subset
from utils.datasets.res2mol import Res2MolDataset
import torch 
dataset = Res2MolDataset()
base_path = '/home/haotian/molecules_confs/Protein_test/Res2mol'
split_by_name = torch.load(os.path.join(base_path,'data/split_by_name.pt'))
split = {
    k: [dataset.name2id[n] for n in names if n in dataset.name2id]
    for k, names in split_by_name.items()
}
subsets = {k:Subset(dataset, indices=v) for k, v in split.items()}
train_set, val_set = subsets['train'], subsets['test']

import os.path as osp
from tqdm import tqdm
crossdock= './data/crossdocked_pocket10'
fps = []
for i in tqdm(range(len(train_set))):
    ligand = osp.join(crossdock,train_set[i].ligand_filename)
    train_mol = read_sdf(ligand)[0]
    fp = Chem.RDKFingerprint(train_mol)
    fps.append(fp)
write_pkl(fps,'./fps/fps/trainmols_fps.pkl')


#compute the generated mols vs Train set
fps = {}
train_fps = read_pkl('./fps/fps/trainmols_fps.pkl')
fps['GraphBP'] = read_pkl('./fps/fps/GraphBP_valmols_fps.pkl')
fps['pkt2mol'] = read_pkl('./fps/fps/pkt2mol_valmols_fps.pkl')
fps['ResGen'] = read_pkl('./fps/fps/ResGen_valmols_fps.pkl')
import multiprocessing
from functools import partial 

def compute_sim_for_trian(i,method_fps):
    smis = []
    for j in range(len(train_fps)):
        smi = DataStructs.TanimotoSimilarity(method_fps[i], train_fps[j])
        smis.append(smi)
    return smis

methods = list(fps.keys())
for method in methods:
    pool = multiprocessing.Pool(12)
    func = partial(compute_sim_for_trian, method_fps=fps[method])

    all_sims = []
    for result in tqdm(pool.imap(func, range(len(fps[method]))), total=len(fps[method])):
        all_sims.append(result)
    pool.join()
    pool.close()
    write_pkl(all_sims, f'./fps/{method}_train_sim.pkl')


#compute the generated mols in the Val set mean(max)
def docking_analysis(i, dirs):
    dir_ = dirs[i]
    ori_result = read_pkl(os.path.join(dir_,'ori_docking_results.pkl'))
    ori_aff = ori_result[0][0]['affinity']
    ori_mol = ori_result[0][0]['rdmol']
    gen_results = read_pkl(os.path.join(dir_, 'gen_docking_results.pkl'))
    scores = []
    mols = []
    for i in range(len(gen_results)):
        try:
            scores.append(gen_results[i][0]['affinity'])
            mols.append(gen_results[i][0]['rdmol'])
        except:
            scores.append(0)
            mols.append(0)
    scores_zip = zip(np.sort(scores),np.argsort(scores))
    scores = np.sort(scores)
    percent = len(scores[scores<ori_aff])/len(scores)
    return percent, list(scores_zip), ori_aff, dir_, mols, ori_mol

from glob import glob
gen_mols = {}
ori_mols = {}
methods = ['GraphBP','pkt2mol','resgen_epo57','resgen_epo187']
for method in methods:
    gen_mols[method] = []
    ori_mols[method] = []
    dirs = glob(f'./{method}/*')
    for i in tqdm(range(len(dirs))):
        try:
            perc, sco, ori_a, name, mols, ori_mol = docking_analysis(i, dirs)
            gen_mols[method].extend(mols)
            ori_mols[method].extend(ori_mol)
        except:
            ...
    sims = compute_sims(gen_mols, ori_mols, source='mol')
    write_pkl(sims,f'{method}_val_sim.pkl')


#compute the generated mols in the DUD-E(in signle target)
ori_mols = read_pkl('./ori_drugs.pkl')
mols = {}
mols['GraphBP'] = read_pkl('./GraphBP.pkl')
mols['pkt2mol'] = read_pkl('./pkt2mol_gen.pkl')
mols['ResGen'] = read_pkl('./resgen_gen.pkl')
mols['random'] = read_pkl('./random_drugs.pkl')
methods = list(mols.kesy())
for method in methods:
    sims = compute_sims(mols[method], ori_mols, source='mol')
    write_pkl(sims,f'{method}_ori_sim.pkl')
