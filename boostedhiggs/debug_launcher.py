import os
import numpy as np
from coffea import processor, util, hist
from coffea.util import load, save
#%matplotlib inline
import matplotlib.pyplot as plt
# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoaod import NanoEvents
from hbbprocessor import HbbProcessor
from ddt_processor import DDTProcessor
from zqq_processor import ZQQProcessor
from coffea.nanoaod.methods import Candidate
events = NanoEvents.from_file(
    #'root://cmseos.fnal.gov//store/user/jkrupa/nanopost_process/24Jul20/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_9ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/group/lpcbacon/jkrupa/nanopost_process/6Aug20//WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_9WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/nano_mc_2017_9QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/nano_mc_2017_2QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/nano_mc_2017_15ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/nano_mc_2017_109TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root',
    'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/nano_mc_2017_1TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/nano_mc_2017_1QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_1ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20//ZJetsToQQ_HT800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_9ZJetsToQQ_HT800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20_v2/SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1/nano_data_2017_21SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1.root',
    ##'root://cmseos.fnal.gov//store/group/lpcbacon/jkrupa/nanopost_process/6Aug20/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/nano_mc_2017_99TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/lpcbacon/jkrupa/nanopost_process/6Aug20/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_1WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/jkrupa/nanopost_process/27Jul20_v3/ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_9ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/jkrupa/nanopost_process/27Jul20_v3/ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_98ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8.root',
    #'root://cmseos.fnal.gov//store/user/jkrupa/nanopost_process/22Jul20_v2/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/nano_mc_2017_9_Skim.root', 
    #entrystop=100000, 
    #metadata={'dataset': 'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8'},#,QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-test'},
    #metadata={'dataset': 'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-test'},
    metadata={'dataset':'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'},

    #metadata={'dataset': 'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8'},
    #metadata={'dataset': 'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'},
    #methods={"FatJetPFCands": Candidate}
)
#p = HbbProcessor(year='2017')
p = ZQQProcessor(year='2017',region='signal')
#p = DDTProcessor(year='2017')
out = p.process(events)
print(out)
#print('ZQQProcessor',out['sumw'],out['templates'].sum('dataset','pt','msd','mu_pt','in_v3_ddt','mu_pfRelIso04_all',overflow='allnan').values()) #'in_v3_ddt','Vmatch','gruddt','mu_pt','mu_pfRelIso04_all',).values())
#x = out['templates'].sum('dataset','gruddt','pt','msd','in_v3_ddt','mu_pfRelIso04_all',overflow='allnan').integrate('region','ttbar_muoncontrol')
#x = out['templates'].sum('dataset').integrate('region','ttbar_muoncontrol')
#hist.plot1d(x,overflow='all')
#plt.savefig('test_muonisolation.pdf')
save(out, 'test.coffea')

#p = HbbProcessor(year='2017')
#out = p.process(events)
#print('HbbProcessor',out['sumw'],out['templates'].sum('dataset','pt','msd',).values())
