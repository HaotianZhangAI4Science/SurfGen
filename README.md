# SurfGen

SurfGen: Learning on Topological Surface and Geometric Structure for 3D Molecular Generation



## Data Preparation



### Making surface data on your own. 

#### Approach 1

Although we have prepared the required data for training and evaluation at [data_link](needed to be uploaded). But you may want to apply SurfGen in your own case. So we provide the guideline for creating the surf_maker environment.

`conda create -n surf_maker pymesh2 jupyter scipy joblib biopython rdkit` plyfile

We highly recommend using mamba instead of conda for speeding up. 

`mamba create -n surf_maker pymesh2 jupyter scipy joblib biopython rdkit` plyfile

#### Approach 2

We also provide the .yml file for creating environment

`conda env create -f surf_maker_environment.yml`

#### Install APBS toolkits

When the base python environment was created, then install [APBS-3.0.0](https://github.com/Electrostatics/apbs/releases), [pdb2pqr-2.1.1](https://github.com/Electrostatics/apbs-pdb2pqr/releases) on your computer.Then set the msms_bin, apbs_bin, pdb2pqr_bin, multivalue_bin path directly in your ~/.bashrc, or just set them in the scripts when creating the surface file from the pdb file.  

#### Try Generate surface now !

Now you have deployed all the dependent environments. Please follow the ./data/surf_maker for making surface data. Or run the ./data/surf_maker/surf_maker_test.py for testing whether you have figured out this environment successfully. 

`python ./data/surf_maker/surf_maker_test.py`



