import os
import numpy as np
from coffea import processor, util, hist

#%matplotlib inline
import matplotlib.pyplot as plt
# import mplhep
# plt.style.use(mplhep.style.ATLAS)
from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoaod import NanoEvents
from ddt_processor import DDTProcessor
from coffea.nanoaod.methods import Candidate
events = NanoEvents.from_file(
    'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2/200127_174245/0000/nano_mc_2017_1.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISumer19UL17MiniAOD-106X_v6-v2/200127_174245/0000/nano_mc_2017_1.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2/200207_204533/0000/nano_mc_2017_11.root', 
    entrystop=50000, 
    metadata={'dataset': 'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-test'},
    methods={"FatJetPFCands": Candidate}
)
p = DDTProcessor(year='2017')
out = p.process(events)
#ax = hist.plot1d(out['cutflow'].integrate('dataset').integrate('region'))
#plt.savefig('test2.pdf')

gru = out['jet_kin'].integrate('dataset').integrate('region').integrate('jet_pt').integrate('jet_rho')
plt.clf()
ax = hist.plot1d(gru)
plt.savefig('test1.pdf')

