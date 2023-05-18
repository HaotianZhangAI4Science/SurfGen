# suppose: 
# node_sca dim=4, edge_sca_dim=6
# node_vec_dim = 3, edge_vec_dim = 3 
from torch import nn
from torch_scatter import scatter_sum
from math import pi 
from ..invariant import VNLinear, GVPerceptronVN, GVLinear
import torch
from ..model_utils import  GaussianSmearing

class EdgeMapping(nn.Module):
    def __init__(self, edge_channels):
        super().__init__()
        self.nn = nn.Linear(in_features=1, out_features=edge_channels, bias=False)
    
    def forward(self, edge_vector):
        edge_vector = edge_vector / (torch.norm(edge_vector, p=2, dim=1, keepdim=True)+1e-7)
        expansion = self.nn(edge_vector.unsqueeze(-1)).transpose(1, -1)
        return expansion
        
class Geodesic_GNN(nn.Module):
    def __init__(self, node_sca_dim=256, node_vec_dim=64, hid_dim=128, edge_dim=64, num_edge_types=2, \
        out_sca_dim=256, out_vec_dim=64, cutoff = 10.):
        super().__init__()
        self.cutoff = cutoff
        self.edge_expansion = EdgeMapping(edge_dim)
        self.distance_expansion = GaussianSmearing(stop=cutoff, num_gaussians=edge_dim - num_edge_types)

        self.node_mapper = GVLinear(node_sca_dim,node_vec_dim,node_sca_dim,node_vec_dim)
        self.edge_mapper = GVLinear(edge_dim,edge_dim,node_sca_dim,node_vec_dim)

        self.edge_sca_sca = nn.Linear(node_sca_dim, hid_dim)
        self.node_sca_sca = nn.Linear(node_sca_dim, hid_dim)

        self.edge_sca_vec = nn.Linear(node_sca_dim, hid_dim)
        self.node_sca_vec = nn.Linear(node_sca_dim, hid_dim)
        self.edge_vec_vec = VNLinear(node_vec_dim, hid_dim)
        self.node_vec_vec = VNLinear(node_vec_dim, hid_dim)

        self.msg_out = GVLinear(hid_dim, hid_dim, out_sca_dim, out_vec_dim)

        self.resi_connecter = GVLinear(node_sca_dim,node_vec_dim,node_sca_dim,node_vec_dim)
        self.aggr_out = GVLinear(node_sca_dim,node_vec_dim,node_sca_dim,node_vec_dim)
    
    def forward(self, node_feats, edge_feature, edge_vector, edge_index, gds_dist):
        
        num_nodes = node_feats[0].shape[0]
        ## map edge_fetures: original space -> interaction space
        edge_sca_feat = torch.cat([self.distance_expansion(gds_dist), edge_feature], dim=-1)
        edge_vec_feat = self.edge_expansion(edge_vector) 

        # Geodesic Message Passing 
        ## mapping the node and edge features to the same space 
        edge_index_row = edge_index[0]
        node_sca_feats, node_vec_feats = self.node_mapper(node_feats)
        edge_sca_feat, edge_vec_feat = self.edge_mapper([edge_sca_feat, edge_vec_feat])
        node_sca_feats, node_vec_feats = node_sca_feats[edge_index_row], node_vec_feats[edge_index_row]
        ## geodesic coefficient 
        coeff = 0.5 * (torch.cos(gds_dist * pi / self.cutoff) + 1.0)
        coeff = coeff * (gds_dist <= self.cutoff) * (gds_dist >= 0.0)
        ## compute the scalar message
        msg_sca_emb = self.node_sca_sca(node_sca_feats) * self.edge_sca_sca(edge_sca_feat)
        msg_sca_emb = msg_sca_emb * coeff.view(-1,1)

        ## compute the vector message
        msg_vec_emb1 = self.node_vec_vec(node_vec_feats) * self.edge_sca_vec(edge_sca_feat).unsqueeze(-1)
        msg_vec_emb2 = self.node_sca_vec(node_sca_feats).unsqueeze(-1) * self.edge_vec_vec(edge_vec_feat)
        msg_vec_emb = msg_vec_emb1 + msg_vec_emb2
        msg_vec_emb = msg_vec_emb * coeff.view(-1,1,1)
        ## message pssing mapping 
        msg_sca_emb, msg_vec_emb = self.msg_out([msg_sca_emb, msg_vec_emb])

        ## aggregate the message 
        aggr_msg_sca = scatter_sum(msg_sca_emb, edge_index_row, dim=0, dim_size=num_nodes)
        aggr_msg_vec = scatter_sum(msg_vec_emb, edge_index_row, dim=0, dim_size=num_nodes)

        ## residue connection
        resi_sca, resi_vec = self.resi_connecter(node_feats)
        out_sca = resi_sca + aggr_msg_sca
        out_vec = resi_vec + aggr_msg_vec

        ## aggregation mapper
        out_sca, out_vec = self.aggr_out([out_sca, out_vec])
        
        return [out_sca, out_vec]



# class Geodesic_GNN(nn.Module):
#     def __init__(self, node_sca_dim=4, node_vec_dim=3, edge_sca_dim=6, edge_vec_dim=3, out_sca=16, out_vec=16, cutoff=10.):
#         super().__init__()
#         # To simplify the model, the out_feats_dim of edges and nodes are the same
        
#         self.self.edge_sca_sca = nn.Linear(edge_sca_dim, out_sca)
#         self.self.node_sca_sca = nn.Linear(node_sca_dim, out_sca)

#         self.self.edge_sca_vec = nn.Linear(edge_sca_dim, out_sca)
#         self.self.node_sca_vec = nn.Linear(node_sca_dim, out_sca)
#         self.self.edge_vec_vec = VNLinear(edge_vec_dim, out_vec)
#         self.self.node_vec_vec = VNLinear(node_vec_dim, out_vec)
        
#         self.encoder = GVPerceptronVN(out_sca,out_vec,out_sca,out_sca)
    
#     def forward(self, node_feats, edge_feats, edge_index, gds_dist):
#         edge_index_raw = edge_index[0]
#         num_nodes = node_feats[0].shape[0]
#         coeff = 0.5 * (torch.cos(gds_dist * pi / self.cutoff) + 1.0)
#         coeff = coeff * (gds_dist <= self.cutoff) * (gds_dist >= 0.0)

#         # compute the scalar message
#         msg_sca_emb = self.self.node_sca_sca(node_feats[0])[edge_index_raw] * self.self.edge_sca_sca(edge_feats[0])
#         msg_sca_emb = msg_sca_emb * coeff.view(-1,1)
        
#         # compute the vector message
#         msg_vec_emb1 = self.self.node_vec_vec(node_feats[1])[edge_index_raw] * self.self.edge_sca_vec(edge_feats[0]).unsqueeze(-1)
#         msg_vec_emb2 = self.self.node_sca_vec(node_feats[0])[edge_index_raw].unsqueeze(-1) * self.self.edge_vec_vec(edge_feats[1])
#         msg_vec_emb = msg_vec_emb1 + msg_vec_emb2
#         msg_vec_emb = msg_vec_emb * coeff.view(-1,1,1)

#         # arrgrate the message 
#         aggr_msg_sca = scatter_sum(msg_sca_emb, edge_index_raw, dim=0, dim_size=num_nodes)
#         aggr_msg_vec = scatter_sum(msg_vec_emb, edge_index_raw, dim=0, dim_size=num_nodes)
        
#         # then encode the geodesic feates 
#         node_aggr_sca, node_aggr_vec = self.encoder((aggr_msg_sca,aggr_msg_vec))
        
#         return node_aggr_sca, node_aggr_vec