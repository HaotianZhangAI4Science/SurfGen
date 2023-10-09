import pickle
import torch
import os
import numpy as np
from torch.utils.data import Subset
from .dataset import SurfLigandPairDataset


def get_dataset(config, *args, **kwargs):
    index_path = config.index_path
    processed_path = config.processed_path
    name2id = config.name2id_path
    dataset = SurfLigandPairDataset(index_path,processed_path,name2id, *args, **kwargs)

    if 'split' in config:
        split_by_name = torch.load(config.split)
        split = {
            k: [dataset.name2id[n] for n in names if n in dataset.name2id]
            for k, names in split_by_name.items()
        }
        subsets = {k:Subset(dataset, indices=v) for k, v in split.items()}
        return dataset, subsets
    else:
        return dataset