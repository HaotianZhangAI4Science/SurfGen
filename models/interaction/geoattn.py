# supnode_pose: 
# node_sca dim=4, edge_sca_dim=6
# node_vec_dim = 3, edge_vec_dim = 3 
from torch import nn
from torch_scatter import scatter_sum
from ..invariant import VNLinear, GVPerceptronVN
import torch
from torch_scatter import scatter_softmax
from torch.nn import Sigmoid
from ..model_utils import  GaussianSmearing

class EdgeMapping(nn.Module):
    def __init__(self, edge_channels):
        super().__init__()
        self.nn = nn.Linear(in_features=1, out_features=edge_channels, bias=False)
    
    def forward(self, edge_vector):
        edge_vector = edge_vector / (torch.norm(edge_vector, p=2, dim=1, keepdim=True)+1e-7)
        expansion = self.nn(edge_vector.unsqueeze(-1)).transpose(1, -1)
        return expansion

class Geoattn_GNN(nn.Module):
    def __init__(self, node_sca_dim=256,node_vec_dim=64, num_edge_types=4, edge_dim=64, hid_dim=128,\
        out_sca_dim=256, out_vec_dim=64, cutoff=10):
        super().__init__()
        
        self.cutoff = cutoff
        self.edge_expansion = EdgeMapping(edge_dim)
        self.distance_expansion = GaussianSmearing(stop=cutoff, num_gaussians=edge_dim - num_edge_types)
        self.node_mapper = GVLinear(node_sca_dim,node_vec_dim,node_sca_dim,node_vec_dim)
        self.edge_mapper = GVLinear(edge_dim,edge_dim,node_sca_dim,node_vec_dim)

        self.edge_net = nn.Linear(node_sca_dim, hid_dim)
        self.node_net = nn.Linear(node_sca_dim, hid_dim)

        self.edge_sca_net = nn.Linear(node_sca_dim, hid_dim)
        self.node_sca_net = nn.Linear(node_sca_dim, hid_dim)
        self.edge_vec_net = VNLinear(node_vec_dim, hid_dim)
        self.node_vec_net = VNLinear(node_vec_dim, hid_dim)


        self.sca_attn_net = nn.Linear(node_sca_dim*2+1, hid_dim)
        self.vec_attn_net = VNLinear(node_vec_dim, hid_dim)
        self.softmax = scatter_softmax  
        self.sigmoid = nn.Sigmoid()

        self.msg_out = GVLinear(hid_dim, hid_dim, out_sca_dim, out_vec_dim)

        self.resi_connecter = GVLinear(node_sca_dim,node_vec_dim,node_sca_dim,node_vec_dim)
        self.aggr_out = GVPerceptronVN(node_sca_dim,node_vec_dim,node_sca_dim,node_vec_dim)
    
    def forward(self, node_feats, edge_feature, edge_vector, edge_index):
        num_nodes = node_feats[0].shape[0]
        edge_dist = torch.norm(edge_vector, dim=-1)
        edge_index_row = edge_index[0]

        ## map edge_features: original space -> interation space
        edge_dist = torch.norm(edge_vector, dim=-1, p=2)
        edge_sca_feat = torch.cat([self.distance_expansion(edge_dist), edge_feature], dim=-1)
        edge_vec_feat = self.edge_expansion(edge_vector) 

        # message passing framework
        ## extract edge and node features in interaction space
        node_sca_feats, node_vec_feats = self.node_mapper(node_feats)
        edge_sca_feat, edge_vec_feat = self.edge_mapper([edge_sca_feat, edge_vec_feat])
        node_sca_feats, node_vec_feats = node_sca_feats[edge_index_row], node_vec_feats[edge_index_row]

        ## compute the attention score \alpha_ij and A_ij
        alpha_sca = torch.cat([node_sca_feats[edge_index[0]], node_sca_feats[edge_index[1]], edge_dist.unsqueeze(-1)], dim=-1)
        alpha_sca = self.sca_attn_net(alpha_sca)
        alpha_sca = self.softmax(alpha_sca,edge_index_row,dim=0)

        alpha_vec_hid = self.vec_attn_net(node_vec_feats)
        alpha_vec = (alpha_vec_hid[edge_index[0]] * alpha_vec_hid[edge_index[1]]).sum(-1).sum(-1)
        alpha_vec = self.sigmoid(alpha_vec)
    
        ## message: the scalar feats
        node_sca_feat =  self.node_net(node_sca_feats)[edge_index_row] * self.edge_net(edge_sca_feat) 
        ## message: the equivariant interaction between node feature and edge feature
        node_sca_hid = self.node_sca_net(node_sca_feats)[edge_index_row].unsqueeze(-1)
        edge_vec_hid = self.edge_vec_net(edge_vec_feat)
        node_vec_hid = self.node_vec_net(node_vec_feats)[edge_index_row]
        edge_sca_hid =  self.edge_sca_net(edge_sca_feat).unsqueeze(-1)
        msg_sca = node_sca_feat * alpha_sca 
        msg_vec = (node_sca_hid * edge_vec_hid + node_vec_hid*edge_sca_hid)*alpha_vec.unsqueeze(-1).unsqueeze(-1)
        msg_sca,msg_vec = self.msg_out([msg_sca,msg_vec])
        
        ## aggregate the message 
        aggr_msg_sca = scatter_sum(msg_sca, edge_index_row, dim=0, dim_size=num_nodes)
        aggr_msg_vec = scatter_sum(msg_vec, edge_index_row, dim=0, dim_size=num_nodes)

        ## residue connection 
        resi_sca, resi_vec = self.resi_connecter(node_feats)
        out_sca = resi_sca + aggr_msg_sca
        out_vec = resi_vec + aggr_msg_vec

        ## map the aggregated feature
        out_sca, out_vec = self.aggr_out([out_sca, out_vec])

        return [out_sca, out_vec]


# class Geoattn_GNN(nn.Module):
#     def __init__(self, input_node_vec_dim=2, node_vec_dim=3,input_node_sca_dim=13, \
#         input_edge_vec_dim = 1, input_edge_sca_dim=4, out_dim=16, normalize=20.):
#         super().__init__()
#         # To simplify the model, the out_feats_dim of edges and nodes are the same
        
#         ### vector feature mapping 
#         self.node_vec_net = VNLinear(input_node_vec_dim,out_dim)
#         self.node_sca_net = nn.Linear(input_node_sca_dim, out_dim)
#         self.edge_vec_net = VNLinear(input_edge_vec_dim, out_dim)
#         self.edge_sca_net = nn.Linear(input_edge_sca_dim, out_dim)
        
#         ### scalar feature mapping 
#         self.node_net = nn.Linear(input_node_sca_dim, out_dim)
#         self.edge_net = nn.Linear(input_edge_sca_dim, out_dim)

#         ### the input dimension of sca_attn is (node_j||node_i||dist_ij)
#         self.sca_attn_net = nn.Linear(input_node_sca_dim*2+1, out_dim)
#         self.vec_attn_net = nn.Linear(input_node_vec_dim, out_dim)
#         self.softmax = scatter_softmax  
#         self.sigmoid = Sigmoid()

#         self.mapper = GVPerceptronVN(out_dim,out_dim,out_dim,out_dim)
    
#     def forward(self, node_feats, edge_feats, edge_index, node_node_pos):
        
#         edge_dist = torch.norm(node_pos[edge_index[0]]-node_pos[edge_index[1]], dim=-1)
#         node_sca, node_vec = node_feats
#         edge_sca, edge_vec = edge_feats
#         edge_index_raw = edge_index[0]

#         # compute the attention score \alpha_ij and A_ij
#         alpha_sca = torch.cat([node_sca[edge_index[0]], node_sca[edge_index[1]], edge_dist.unsqueeze(-1)], dim=-1)
#         alpha_sca = self.sca_attn_net(alpha_sca)
#         alpha_sca = self.softmax(alpha_sca,edge_index_raw,dim=0)

#         alpha_vec_hid = self.vec_attn_net(node_vec)
#         alpha_vec = (alpha_vec_hid[edge_index[0]] * alpha_vec_hid[edge_index[1]]).sum(-1).sum(-1)
#         alpha_vec = self.sigmoid(alpha_vec)
        
#         # the scalar feats
#         node_sca_feat =  self/node_net(node_sca)[edge_index_raw] * edge_net(edge_sca) * alpha_sca 
#         # the equivariant interaction between node feature and edge feature
#         node_sca_hid = self.node_sca_net(node_sca)[edge_index_raw].unsqueeze(-1)
#         edge_vec_hid = self.edge_vec_net(edge_vec)
#         node_vec_hid = self.node_vec_net(node_vec)[edge_index_raw]
#         edge_sca_hid =  self.edge_sca_net(edge_sca).unsqueeze(-1)
        
#         emb_sca = scatter_sum(node_sca_feat,edge_index_raw, dim=0)
#         emb_vec = scatter_sum((node_sca_hid * edge_vec_hid + node_vec_hid*edge_sca_hid)*alpha_vec.unsqueeze(-1).unsqueeze(-1), edge_index_raw, dim=0)
        
#         ### perform the non-linear transformation between scalar feature and vector feature
#         out = self.mapper([emb_sca,emb_vec])

#         return out
