# :loudspeaker: SurfGen: Learning on Topological Surface and Geometric Structure for 3D Molecular Generation

<div align=center>
<img src="./assets/toc.png" width="50%" height="50%" alt="TOC" align=center />
</div>

## Environment

Install via conda .yml file (cuda 11.3)

```python
conda install mamba
mamba env create -f surfgen_environment.yml -n surfgen
conda activate surfgen 
```

If you're reluctant to use mamba: 

```python
conda env create -f surfgen_environment.yml -n surfgen
```

We also provide conda-packed file [here](https://doi.org/10.5281/zenodo.7758282). Download it and then unzip it in your conda/envs/dir. For me, the directory is ~/.conda/envs. Special thanks to the creators and organizers of zenodo, which provides a free platform to store large files for academic use. 

```
mkdir ~/.conda/envs/surfgen
tar -xzvf surfgen.tar.gz -C ~/.conda/envs/surfgen
conda activate surfgen
```

## Data

The main data used for training is CrossDock2020 

#### Download the data from the original source

```python
wget https://bits.csb.pitt.edu/files/crossdock2020/CrossDocked2020_v1.1.tgz -P data/crossdock2020/
tar -C data/crossdock2020/ -xzf data/crossdock2020/CrossDocked2020_v1.1.tgz
wget https://bits.csb.pitt.edu/files/it2_tt_0_lowrmsd_mols_train0_fixed.types -P data/crossdock2020/
wget https://bits.csb.pitt.edu/files/it2_tt_0_lowrmsd_mols_test0_fixed.types -P data/crossdock2020/
```

Then following the guideline to process it.  The train data split is [split_name.pt](https://drive.google.com/file/d/1WUZVNv--gztqDNoA3BEexXdjRfXKHuHn/view?usp=share_link). 

If it's inconvenient for you, we also provided the [processed data](https://drive.google.com/file/d/1WEbPe7XsOCEmDrgMCrD7ZrkBmkFQ7KkK/view?usp=share_link), [processed_data_key](https://drive.google.com/file/d/1Ko_yDF14ck4Y73Tgnbt9jS-R68n0pzOT/view?usp=share_link), [name2id](https://drive.google.com/file/d/1zte17O7HZuNYD56wrvNzu-SSNjmLDE2r/view?usp=share_link). You just need to download them in ./data  and create a ./data/crossdock_pocket10 directory, and put the [index.pkl](https://drive.google.com/file/d/1-YCXOV-MWDOE-p6laQxOKPLPVJRakpL1/view?usp=share_link) in it. 



### (Optional) Making surface data on your own. 

#### Create the base python environment

##### Approach 1

Although we have prepared the required data for training and evaluation above. But you may want to apply SurfGen in your own case. So we provide the guideline for creating the surf_maker environment.

```python
conda create -n surf_maker pymesh2 jupyter scipy joblib biopython rdkit plyfile -c conda-forge
```

We highly recommend using mamba instead of conda for speeding up. 

```python
mamba create -n surf_maker pymesh2 jupyter scipy joblib biopython rdkit plyfile -c conda-forge
```

##### Approach 2

We also provide the .yml file for creating environment

```
conda env create -f surf_maker_environment.yml
```

#### Install APBS Toolkits

When the base python environment was created, then install [APBS-3.0.0](https://github.com/Electrostatics/apbs/releases), [pdb2pqr-2.1.1](https://github.com/Electrostatics/apbs-pdb2pqr/releases) on your computer. Then set the msms_bin, apbs_bin, pdb2pqr_bin, and multivalue_bin path directly in your ~/.bashrc, or just set them in the scripts when creating the surface file from the pdb file.  

#### Try Generate Surface Now !

Now you have deployed all the dependent environments. Please follow the ./data/surf_maker for making surface data. Or run the ./data/surf_maker/surf_maker_test.py for testing whether you have figured out this environment successfully. 

```python
python ./data/surf_maker/surf_maker_test.py
```



## Generation 

To run the generation, run the file. The model's parameters can be downloaded [here](https://drive.google.com/file/d/1tKIib7qRN5IdNXhVB_8v1EBYiWdU12v8/view?usp=share_link). Put it at ./ckpt. 

```python
python gen.py --outdir example --check_point ./ckpt/val_119.pt --ply_file ./example/3cl_pocket_8.0_res_1.5.ply
```



## Training 

```python
python train.py
```

