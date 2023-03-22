# suppose: 
# node_sca dim=4, edge_sca_dim=6
# node_vec_dim = 3, edge_vec_dim = 3 
from torch import nn
from torch_scatter import scatter_sum
from math import pi 
from models.invariant import VNLinear, GVPerceptronVN
import torch


class Geodesic_GNN(nn.Module):
    def __init__(self, node_sca_dim=4, node_vec_dim=3, edge_sca_dim=6, edge_vec_dim=3, out_sca=16, out_vec=16, cutoff=10.):
        super().__init__()
        # To simplify the model, the out_feats_dim of edges and nodes are the same
        
        self.edge_sca_sca = nn.Linear(edge_sca_dim, out_sca)
        self.node_sca_sca = nn.Linear(node_sca_dim, out_sca)

        self.edge_sca_vec = nn.Linear(edge_sca_dim, out_sca)
        self.node_sca_vec = nn.Linear(node_sca_dim, out_sca)
        self.edge_vec_vec = VNLinear(edge_vec_dim, out_vec)
        self.node_vec_vec = VNLinear(node_vec_dim, out_vec)
        
        self.encoder = GVPerceptronVN(out_sca,out_vec,out_sca,out_sca)
    
    def forward(self, node_feats, edge_feats, edge_index, gds_dist):
        edge_index_raw = edge_index[0]
        num_nodes = node_feats[0].shape[0]
        coeff = 0.5 * (torch.cos(gds_dist * pi / self.cutoff) + 1.0)
        coeff = coeff * (gds_dist <= self.cutoff) * (gds_dist >= 0.0)

        # compute the scalar message
        msg_sca_emb = self.node_sca_sca(node_feats[0])[edge_index_raw] * self.edge_sca_sca(edge_feats[0])
        msg_sca_emb = msg_sca_emb * coeff.view(-1,1)
        
        # compute the vector message
        msg_vec_emb1 = self.node_vec_vec(node_feats[1])[edge_index_raw] * self.edge_sca_vec(edge_feats[0]).unsqueeze(-1)
        msg_vec_emb2 = self.node_sca_vec(node_feats[0])[edge_index_raw].unsqueeze(-1) * self.edge_vec_vec(edge_feats[1])
        msg_vec_emb = msg_vec_emb1 + msg_vec_emb2
        msg_vec_emb = msg_vec_emb * coeff.view(-1,1,1)

        # arrgrate the message 
        aggr_msg_sca = scatter_sum(msg_sca_emb, edge_index_raw, dim=0, dim_size=num_nodes)
        aggr_msg_vec = scatter_sum(msg_vec_emb, edge_index_raw, dim=0, dim_size=num_nodes)
        
        # then encode the geodesic feates 
        node_aggr_sca, node_aggr_vec = self.encoder((aggr_msg_sca,aggr_msg_vec))
        
        return node_aggr_sca, node_aggr_vec