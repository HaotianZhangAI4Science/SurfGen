import numpy as np
import networkx as nx
from plyfile import PlyData
import torch
from torch_geometric.data import Data
import torch

def geodesic_matrix(pos,edge_index):
    '''
    edge_index.shape=[2,edges]
    '''
    # In fact, if the max_threshold for geodesic distance computation could be provided,
    # the speed would be dramastically acclearted.  
    if type(pos) == torch.Tensor:
        pos = np.array(pos)
    if type(edge_index) == torch.Tensor:
        edge_index = np.array(edge_index)
    nodes = np.array(range(pos.shape[0]))
    dist = np.linalg.norm(pos[edge_index[0]] - pos[edge_index[1]], axis=-1)
    weighted_list = list(zip(edge_index[0],edge_index[1],dist))
    Graph = nx.Graph()
    Graph.add_nodes_from(nodes)
    Graph.add_weighted_edges_from(weighted_list)
    shortest_path_mat = np.zeros((len(nodes),len(nodes)))
    for i in range(len(nodes)):
        raw = nx.shortest_path_length(Graph,source=i,weight="weight")
        for j in range(len(nodes)):
            shortest_path_mat[i][j] = raw[j]
    return shortest_path_mat


def dst2knnedge(dst_mat, num_knn=24, self_loop=False):
    knn_edge_index_src = []
    knn_edge_index_tgt = []
    knn_edge_dist = []
    num_nodes = dst_mat.shape[0]
    for node_idx in range(num_nodes):
        knn_edge_index_src.extend([node_idx]*num_knn)
        
        if self_loop:
            knn_edge_index_tgt.extend(np.argsort(dst_mat[node_idx])[:num_knn])
            knn_edge_dist.extend(np.sort(dst_mat[node_idx])[:num_knn])
        else:
            knn_edge_index_tgt.extend(np.argsort(dst_mat[node_idx])[1:num_knn+1])
            knn_edge_dist.extend(np.sort(dst_mat[node_idx])[1:num_knn+1])

    return torch.tensor(np.array([knn_edge_index_src,knn_edge_index_tgt])), torch.tensor(np.array(knn_edge_dist,dtype=np.float32))


def read_ply(path, read_face=None):
    with open(path, 'rb') as f:
        data = PlyData.read(f)

    features = ([torch.tensor(data['vertex'][axis.name]) for axis in data['vertex'].properties if axis.name not in ['nx', 'ny', 'nz'] ])
    pos = torch.stack(features[:3], dim=-1)
    features = torch.stack(features[3:], dim=-1)
    if read_face is not None:
        if 'face' in data:
            faces = data['face']['vertex_indices']
            faces = [torch.tensor(fa, dtype=torch.long) for fa in faces]
            face = torch.stack(faces, dim=-1)
            data = {'feature':features,\
                'pos':pos,
                'face':face}
    else:
        data = {'feature':features,\
            'pos':pos}
    return data

def read_ply_geom(path, read_face=None):
    with open(path, 'rb') as f:
        data = PlyData.read(f)

    features = ([torch.tensor(data['vertex'][axis.name]) for axis in data['vertex'].properties if axis.name not in ['nx', 'ny', 'nz'] ])
    pos = torch.stack(features[:3], dim=-1)
    features = torch.stack(features[3:], dim=-1)
    if read_face is not None:
        if 'face' in data:
            faces = data['face']['vertex_indices']
            faces = [torch.tensor(fa, dtype=torch.long) for fa in faces]
            face = torch.stack(faces, dim=-1)
            data = Data(x=features,pos=pos,face=face)
    else:
        data = Data(x=features,pos=pos)
    return data

def parse_face(path, read_face=None):
    with open(path, 'rb') as f:
        data = PlyData.read(f)

    faces = data['face']['vertex_indices']
    faces = [torch.tensor(fa, dtype=torch.long) for fa in faces]
    face = torch.stack(faces, dim=-1)

    return face