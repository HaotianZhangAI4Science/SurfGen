import torch
STATUS_RUNNING = 'running'
STATUS_FINISHED = 'finished'
STATUS_FAILED = 'failed'
FOLLOW_BATCH = []
collate_exclude_keys = ['ligand_nbh_list']
from torch_geometric.data import Batch
from models.sample import get_next_step
import numpy as np
from utils.surface import read_ply, geodesic_matrix, dst2knnedge, read_ply_geom, parse_face
from utils.data import torchify_dict, ProteinLigandData
from torch_geometric.utils import to_undirected


@torch.no_grad()  # for a protein-ligand
def get_init(data, model, transform, threshold):
    batch = Batch.from_data_list([data], follow_batch=FOLLOW_BATCH, exclude_keys = collate_exclude_keys,) #batch only contains one data

    ### Predict next atoms
    model.eval()
    predicitions = model.sample_init(
        compose_feature = batch.compose_feature.float(),
        compose_pos = batch.compose_pos,
        # idx_ligand = batch.idx_ligand_ctx_in_compose,
        idx_protein = batch.idx_protein_in_compose,
        gds_edge_sca = batch.gds_edge_sca,
        gds_knn_edge_index = batch.gds_knn_edge_index, 
        gds_dist = batch.gds_dist, 
        compose_knn_edge_index = batch.compose_knn_edge_index,
        compose_knn_edge_feature = batch.compose_knn_edge_feature,
        n_samples_pos=-1,
        n_samples_atom=5,
    )
    data = data.to('cpu')
    # no frontier
    if not predicitions[0]:
        data.status = STATUS_FINISHED
        return [data]
    # has frontiers
    data.status = STATUS_RUNNING
    (has_frontier, idx_frontier, p_frontier,
    idx_focal_in_compose, p_focal,
    pos_generated, pdf_pos, abs_pos_mu, pos_sigma, pos_pi,
    element_pred, element_prob, has_atom_prob) = [p.cpu() for p in predicitions]

    while True:
        data_next_list = get_next_step(
            data,
            p_focal = p_focal,
            pos_generated = pos_generated,
            pdf_pos = pdf_pos,
            element_pred = element_pred,
            element_prob = element_prob,
            has_atom_prob = has_atom_prob,
            # ind_pred = ind_pred,
            # ind_prob = ind_prob,
            bond_index = torch.empty([2, 0]),
            bond_type = torch.empty([0]),
            bond_prob = torch.empty([0]),
            transform = transform,
            threshold=threshold
        )
        data_next_list = [data for data in data_next_list if data.is_high_prob]
        if len(data_next_list) == 0:
            if torch.all(pdf_pos < threshold.pos_threshold):
                threshold.pos_threshold = threshold.pos_threshold / 2
                print('Positional probability threshold is too high. Change to %f' % threshold.pos_threshold)
            elif torch.all(p_focal < threshold.focal_threshold):
                threshold.focal_threshold = threshold.focal_threshold / 2
                print('Focal probability threshold is too high. Change to %f' % threshold.focal_threshold)
            elif torch.all(element_prob < threshold.element_threshold):
                threshold.element_threshold = threshold.element_threshold / 2
                print('Element probability threshold is too high. Change to %f' % threshold.element_threshold)
            else:
                print('Initialization failed.')
        else:
            break

    return data_next_list


@torch.no_grad()  # for a protein-ligand
def get_next(data, model, transform, threshold, frontier_threshold=0,freeze=None, anchor=None):
    batch = Batch.from_data_list([data], follow_batch=FOLLOW_BATCH) #batch only contains one data
    ### Predict next atoms
    model.eval()
    predicitions = model.sample(
        compose_feature = batch.compose_feature.float(),
        compose_pos = batch.compose_pos,
        idx_ligand = batch.idx_ligand_ctx_in_compose,
        idx_protein = batch.idx_protein_in_compose,
        gds_edge_sca = batch.gds_edge_sca,
        gds_knn_edge_index = batch.gds_knn_edge_index, 
        gds_dist = batch.gds_dist, 
        compose_knn_edge_index = batch.compose_knn_edge_index,
        compose_knn_edge_feature = batch.compose_knn_edge_feature,
        ligand_context_bond_index = batch.ligand_context_bond_index,
        ligand_context_bond_type = batch.ligand_context_bond_type,
        n_samples_pos=-1,
        n_samples_atom=5,
        frontier_threshold=frontier_threshold,
        freeze=freeze,
        anchor = anchor,
    )
    data = data.to('cpu')
    # no frontier
    if not predicitions[0]:
        data.status = STATUS_FINISHED
        return [data]
    else:
        data.status = STATUS_RUNNING

    # has frontiers
    (has_frontier, idx_frontier, p_frontier,
    idx_focal_in_compose, p_focal,
    pos_generated, pdf_pos, abs_pos_mu, pos_sigma, pos_pi,
    element_pred, element_prob, has_atom_prob,
    bond_index, bond_type, bond_prob) = [p.cpu() for p in predicitions]

    data_next_list = get_next_step(
        data,
        p_focal = p_focal,
        pos_generated = pos_generated,
        pdf_pos = pdf_pos,
        element_pred = element_pred,
        element_prob = element_prob,
        has_atom_prob = has_atom_prob,
        bond_index = bond_index,
        bond_type = bond_type,
        bond_prob = bond_prob,
        transform = transform,
        threshold = threshold
    )
    data_next_list = [data for data in data_next_list if data.is_high_prob]

    return data_next_list

def logp_to_rank_prob(logp, weight=1.0):
    
    logp = [list(p) + [-0.2] * (3 - len(p)) if len(p) != 3 else p for p in logp]  #padding
    logp_sum = np.array([np.sum(l) for l in logp])
    prob = np.exp(logp_sum) + 1
    prob = prob * np.array(weight)
    return prob / prob.sum()

def pdb_to_pocket_data(ply_file):
    '''
    use the sdf_file as the center 
    '''
    protein_dict = torchify_dict(read_ply(ply_file))
    
    data = ProteinLigandData.from_protein_ligand_dicts(
        protein_dict = protein_dict,
        ligand_dict = {
            'element': torch.empty([0,], dtype=torch.long),
            'pos': torch.empty([0, 3], dtype=torch.float),
            'atom_feature': torch.empty([0, 8], dtype=torch.float),
            'bond_index': torch.empty([2, 0], dtype=torch.long),
            'bond_type': torch.empty([0,], dtype=torch.long),
        }
    )

    data.face = parse_face(ply_file)
    edge_index = torch.cat([data.face[:2], data.face[1:], data.face[::2]], dim=1)
    dlny_edge_index = to_undirected(edge_index, num_nodes=data.protein_pos.shape[0])
    gds_mat = geodesic_matrix(data.protein_pos, dlny_edge_index)
    gds_knn_edge_index, gds_knn_edge_dist = dst2knnedge(gds_mat, num_knn=16)
    data.gds_knn_edge_index = gds_knn_edge_index
    data.gds_dist = gds_knn_edge_dist
    data.num_nodes = data.ligand_pos.shape[0] + data.protein_pos.shape[0]
    
    return data