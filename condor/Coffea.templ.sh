#!/bin/sh
#
# Template of the shell script for submitting a CONDOR job
#
# Need to be substituted:
#    - C M S S W B A S E - local release base ($CMSSW_BASE)
#    - D I R E C T O R Y - condor work and output dir
#    - P R E F I X - some generic name for the set of jobs (like ttbar, wjets)
#    - J O B I D  - job number
#
#
#_____ setup the environment ____________________________________________
#

#source /cvmfs/cms.cern.ch/cmsset_default.sh

#get setup coffea area
#source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc9-opt

xrdcp root://cmseos.fnal.gov//store/user/jkrupa/dazsle_coffea.tgz ./dazsle_coffea.tgz
tar -zxf dazsle_coffea.tgz

rm -rf coffeaenv
tar -zxf coffeaenv.tar.gz
ls -ltrh 
source coffeaenv/bin/activate
pip install .
head coffeaenv/bin/activate
python -c "import cloudpickle; print(cloudpickle.__version__)"
which python
echo PYTHONPATH $PYTHONPATH
mkdir test
xrdcp -f root://cmseos.fnal.gov//store/user/jkrupa/SCRIPTNAME test/
python -V
echo python SCRIPTNAME --year YEAR --starti STARTNUM --endi ENDNUM --selsamples SAMPLE --outname ddt_test
python test/SCRIPTNAME --year YEAR --starti STARTNUM --endi ENDNUM --selsamples SAMPLE --outname ddt_test

#move output to eos
env -i xrdcp -f ddt_test.coffea EOSOUT

