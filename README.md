# boostedhiggs

## Quickstart
```bash
git clone git@github.com:jeffkrupa-/boostedhiggs.git
cd boostedhiggs
pip install --user --editable .
```


## Add tensorflow capability
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p miniconda/
source miniconda/bin/activate
conda create -n tf tensorflow
conda activate tf
pip install coffea
conda install -c conda-forge xrootd

# To make the environment transportable
conda install -c conda-forge conda-pack
conda-pack
```


