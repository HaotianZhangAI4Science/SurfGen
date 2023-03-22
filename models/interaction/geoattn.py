# suppose: 
# node_sca dim=4, edge_sca_dim=6
# node_vec_dim = 3, edge_vec_dim = 3 
from torch import nn
from torch_scatter import scatter_sum
from models.invariant import VNLinear, GVPerceptronVN
import torch
from torch_scatter import scatter_softmax
from torch.nn import Sigmoid

class Geoattn_GNN(nn.Module):
    def __init__(self, input_node_vec_dim=2, node_vec_dim=3,input_node_sca_dim=13, \
        input_edge_vec_dim = 1, input_edge_sca_dim=4, out_dim=16, normalize=20.):
        super().__init__()
        # To simplify the model, the out_feats_dim of edges and nodes are the same
        
        ### vector feature mapping 
        self.node_vec_net = VNLinear(input_node_vec_dim,out_dim)
        self.node_sca_net = nn.Linear(input_node_sca_dim, out_dim)
        self.edge_vec_net = VNLinear(input_edge_vec_dim, out_dim)
        self.edge_sca_net = nn.Linear(input_edge_sca_dim, out_dim)
        
        ### scalar feature mapping 
        self.node_net = nn.Linear(input_node_sca_dim, out_dim)
        self.edge_net = nn.Linear(input_edge_sca_dim, out_dim)

        ### the input dimension of sca_attn is (node_j||node_i||dist_ij)
        self.sca_attn_net = nn.Linear(input_node_sca_dim*2+1, out_dim)
        self.vec_attn_net = nn.Linear(input_node_vec_dim, out_dim)
        self.softmax = scatter_softmax  
        self.sigmoid = Sigmoid()

        self.mapper = GVPerceptronVN(out_dim,out_dim,out_dim,out_dim)
    
    def forward(self, node_feats, edge_feats, edge_index, node_pos):
        
        edge_dist = torch.norm(pos[edge_index[0]]-pos[edge_index[1]], dim=-1)
        node_sca, node_vec = node_feats
        edge_sca, edge_vec = edge_feats
        edge_index_raw = edge_index[0]

        # compute the attention score \alpha_ij and A_ij
        alpha_sca = torch.cat([node_sca[edge_index[0]], node_sca[edge_index[1]], edge_dist.unsqueeze(-1)], dim=-1)
        alpha_sca = self.sca_attn_net(alpha_sca)
        alpha_sca = self.softmax(alpha_sca,edge_index_raw,dim=0)

        alpha_vec_hid = self.vec_attn_net(node_vec)
        alpha_vec = (alpha_vec_hid[edge_index[0]] * alpha_vec_hid[edge_index[1]]).sum(-1).sum(-1)
        alpha_vec = self.sigmoid(alpha_vec)
        
        # the scalar feats
        node_sca_feat =  node_net(node_sca)[edge_index_raw] * edge_net(edge_sca) * alpha_sca 
        # the equivariant interaction between node feature and edge feature
        node_sca_hid = self.node_sca_net(node_sca)[edge_index_raw].unsqueeze(-1)
        edge_vec_hid = self.edge_vec_net(edge_vec)
        node_vec_hid = self.node_vec_net(node_vec)[edge_index_raw]
        edge_sca_hid =  self.edge_sca_net(edge_sca).unsqueeze(-1)
        
        emb_sca = scatter_sum(node_sca_feat,edge_index_raw, dim=0)
        emb_vec = scatter_sum((node_sca_hid * edge_vec_hid + node_vec_hid*edge_sca_hid)*alpha_vec.unsqueeze(-1).unsqueeze(-1), edge_index_raw, dim=0)
        
        ### perform the non-linear transformation between scalar feature and vector feature
        out = self.mapper([emb_sca,emb_vec])

        return out


    # class AtomEmbedding(Module):
    # def __init__(self, in_scalar, in_vector,
    #              out_scalar, out_vector, vector_normalizer=20.):
    #     super().__init__()
    #     assert in_vector == 1
    #     self.in_scalar = in_scalar
    #     self.vector_normalizer = vector_normalizer
    #     self.emb_sca = Linear(in_scalar, out_scalar)
    #     self.emb_vec = Linear(in_vector, out_vector)

    # def forward(self, scalar_input, vector_input):
    #     vector_input = vector_input / self.vector_normalizer
    #     assert vector_input.shape[1:] == (3, ), 'Not support. Only one vector can be input'
    #     sca_emb = self.emb_sca(scalar_input[:, :self.in_scalar])  # b, f -> b, f'
    #     vec_emb = vector_input.unsqueeze(-1)  # b, 3 -> b, 3, 1
    #     vec_emb = self.emb_vec(vec_emb).transpose(1, -1)  # b, 1, 3 -> b, f', 3
    #     return sca_emb, vec_emb