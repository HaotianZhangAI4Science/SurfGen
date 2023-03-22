import os
import pickle
import lmdb
import torch
from torch.utils.data import Dataset
from tqdm.auto import tqdm

from ..protein_ligand import PDBProtein, parse_sdf_file
from ..data import ProteinLigandData, torchify_dict


class SurfLigandPairDataset(Dataset):

    def __init__(self, raw_path, transform=None):
        super().__init__()
        self.raw_path = raw_path.rstrip('/')
        self.index_path = os.path.join(self.raw_path, 'index.pkl')
        self.processed_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_processed.lmdb')
        self.name2id_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_name2id.pt')
        self.transform = transform
        self.db = None
        self.keys = None
        if not os.path.exists(self.name2id_path):
            self._precompute_name2id()
        self.name2id = torch.load(self.name2id_path)
        
    def _connect_db(self):
        """
            Establish read-only database connection
        """
        assert self.db is None, 'A connection has already been opened.'
        self.db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=False,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        with self.db.begin() as txn:
            self.keys = list(txn.cursor().iternext(values=False))

    def _precompute_name2id(self):
        name2id = {}
        for i in tqdm(range(self.__len__()), 'Indexing'):
            try:
                data = self.__getitem__(i)
            except AssertionError as e:
                print(i, e)
                continue
            name = (data.protein_filename, data.ligand_filename)
            name2id[name] = i
        torch.save(name2id, self.name2id_path)

    def _close_db(self):
        self.db.close()
        self.db = None
        self.keys = None
    
    def __len__(self):
        if self.db is None:
            self._connect_db()
        return len(self.keys)

    def __getitem__(self, idx):
        if self.db is None:
            self._connect_db()
        key = self.keys[idx]
        data = pickle.loads(self.db.begin().get(key))
        data.id = idx
        assert data.protein_pos.size(0) > 0
        if self.transform is not None:
            data = self.transform(data)
        return data

# class SurfaceLigandPairDataset(Dataset):

#     def __init__(self, raw_path, transform=None):
#         super().__init__()
#         self.raw_path = raw_path
#         self.processed_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_processed.lmdb')
#         self.transform = transform
#         self.db = None
#         self.keys = None
        

#     def _connect_db(self):
#         """
#             Establish read-only database connection
#         """
#         assert self.db is None, 'A connection has already been opened.'
#         self.db = lmdb.open(
#             self.processed_path,
#             map_size=10*(1024*1024*1024),   # 10GB
#             create=False,
#             subdir=False,
#             readonly=True,
#             lock=False,
#             readahead=False,
#             meminit=False,
#         )
#         with self.db.begin() as txn:
#             self.keys = list(txn.cursor().iternext(values=False))

#     def _close_db(self):
#         self.db.close()
#         self.db = None
#         self.keys = None
    
#     def __len__(self):
#         if self.db is None:
#             self._connect_db()
#         return len(self.keys)

#     def __getitem__(self, idx):
#         if self.db is None:
#             self._connect_db()
#         key = self.keys[idx]
#         data = pickle.loads(self.db.begin().get(key))
#         data.id = idx
#         assert data.protein_pos.size(0) > 0
#         if self.transform is not None:
#             data = self.transform(data)
#         return data        

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    args = parser.parse_args()

    PocketLigandPairDataset(args.path)
