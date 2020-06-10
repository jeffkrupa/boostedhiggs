import os
import numpy as np
from coffea import processor, util, hist
from coffea.util import load, save
#%matplotlib inline
import matplotlib.pyplot as plt
# import mplhep
# plt.style.use(mplhep.style.ATLAS)
from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoaod import NanoEvents
from ddt_processor import DDTProcessor
from zqq_processor import ZQQProcessor
from coffea.nanoaod.methods import Candidate
events = NanoEvents.from_file(
    #"root://cmseos.fnal.gov//store/group/lpcbacon/pancakes/02/2017//WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1/200124_202452/0000/nano_mc_2017_99.root",
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/VectorZPrimeToQQ_M100_pT300_TuneCP5_madgraph_pythia8_13TeV/pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1/200124_195054/0000/nano_mc_2017_2.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/VectorZPrimeToQQ_M100_pT300_TuneCP5_madgraph_pythia8_13TeV/pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1/200124_195054/0000/nano_mc_2017_1.root',
    #'root://cmseos.fnal.gov//store/group/lpcbacon/pancakes/02/2017/tmp-VJets-withPF/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1/200504_201723/0000/nano_mc_2017_102.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2/200127_174245/0000/nano_mc_2017_2.root',
    'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2/200127_173321/0000/nano_mc_2017_99.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISumer19UL17MiniAOD-106X_v6-v2/200127_174245/0000/nano_mc_2017_1.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2/200207_204533/0000/nano_mc_2017_11.root', 
    entrystop=100, 
    metadata={'dataset': 'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8QCD'},#,QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-test'},
    methods={"FatJetPFCands": Candidate}
)
#p = ZQQProcessor(year='2017')
#p = ZQQProcessor(year='2017')
p = DDTProcessor(year='2017')
out = p.process(events)

save(out, 'test.coffea')
#ax = hist.plot1d(out['cutflow'].integrate('dataset').integrate('region'))
#plt.savefig('test2.pdf')
#gru = out['jet_kin'].integrate('dataset').integrate('region').integrate('jet_pt').integrate('jet_rho')
#plt.clf()
#ax = hist.plot1d(gru)
#plt.savefig('test1.pdf')

