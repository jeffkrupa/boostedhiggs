# boostedhiggs



## coffea+tensorflow capability (standalone, transportable conda environment)
```wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/nobackup/miniconda.sh -b -p miniconda/
source miniconda/bin/activate
conda create -n env python=3.6
pip install coffea==0.6.38
conda install -c conda-forge xrootd

git clone git@github.com:jeffkrupa-/boostedhiggs.git -b cleanup

cd boostedhiggs

pip install --user --editable .
```
# standalone test
```
cd boostedhiggs
python debug_launcher.py
```
# for job submission
```
cd ../
bash envSetup.sh
cd condor
python CoffeaSubmit.py test run_zqq_processor.py 15 1
```


