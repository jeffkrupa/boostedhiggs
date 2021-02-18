# boostedhiggs


```wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/nobackup/miniconda.sh -b -p miniconda/
source miniconda/bin/activate
conda create -n env python=3.6
pip install coffea
conda install -c conda-forge xrootd

git clone git@github.com:jeffkrupa-/boostedhiggs.git -b ak1

cd boostedhiggs

pip install --user --editable .
```
## standalone test
```
cd boostedhiggs
python debug_launcher.py
```
## for LPC job submission
```
cd ../
bash envSetup.sh
cd condor
voms-proxy-init -voms cms -valid 192:00
python CoffeaSubmit.py test run_zqq_processor.py <# of files per job (~50)> <1 to remake tar (if processor changes)>
python hadd_coffea_eos.py #point to out dir
```


