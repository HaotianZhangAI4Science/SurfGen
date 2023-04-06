import torch
from torch import  nn
from easydict import  EasyDict
import yaml
config = './configs/train_res.yml'
config = load_config(config)

def load_config(path):
    with open(path, 'r') as f:
        return EasyDict(yaml.safe_load(f))

def get_optimizer(cfg, model):
    if cfg.type == 'adam':
        return torch.optim.Adam(
            model.parameters(),
            lr=cfg.lr,
            weight_decay=cfg.weight_decay,
            betas=(cfg.beta1, cfg.beta2, )
        )
    else:
        raise NotImplementedError('Optimizer not supported: %s' % cfg.type)


def get_scheduler(cfg, optimizer):
    if cfg.type == 'plateau':
        return torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            factor=cfg.factor,
            patience=cfg.patience,
            min_lr=cfg.min_lr
        )
    else:
        raise NotImplementedError('Scheduler not supported: %s' % cfg.type)

model = nn.Linear(12,24)
optimizer = get_optimizer(config.train.optimizer, model)
scheduler = get_scheduler(config.train.scheduler, optimizer)