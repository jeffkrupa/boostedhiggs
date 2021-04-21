# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
NanoAODSchema.warn_missing_crossrefs = True
#from ddt_processor import DDTProcessor
from zqq_processor_systs import ZQQProcessor
from coffea.util import save, load
#from test import ZQQProcessor
import uproot

#filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2016/22Mar21_preUL//QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/nano_mc_2016_189_post-process_hadd.root"
filename="root://cmseos.fnal.gov//store/group/lpcbacon/jkrupa/nanopost_process/6Aug20/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/nano_mc_2017_9TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"
#filename="/uscms/home/jkrupa/nobackup/zprlegacy/CMSSW_11_2_1/src/PhysicsTools/NanoAODTools/scripts/test/nano_mc_2017_99_Skim.root"
#filename="root://cmseos.fnal.gov//store/group/lpcpfnano/jkrupa/nanopost_process//WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_9WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root" #
#filename="test/nano_mc_2017_1TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"
#filename="root://cmseos.fnal.gov//store/group/lpcpfnano/jkrupa/nanopost_process/SingleMuon_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1/nano_data_2017_9SingleMuon_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1.root"
#filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2016/22Mar21_preUL//JetHT/nano_data_2016_363_post-process_pancakes-02_Run2016G-17Jul2018-v1_hadd.root"
#filename="test/nano_mc_2016_36_Skim.root"
f = uproot.open(filename)#, xrootd_handler=uproot.MultithreadedXRootDSource)
events = NanoEventsFactory.from_root(
    filename, 
    #entry_stop=100000,
    metadata={'dataset':'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'},
    #metadata={'dataset':'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root'},
    schemaclass=NanoAODSchema,
).events()
p = ZQQProcessor(year='2017',region=['signal','muonCR','VtaggingCR'])
out = p.process(events)
print(out)
p = ZQQProcessor(year='2017',region='signal')
out = p.process(events)
print(out)
p = ZQQProcessor(year='2017',region='muonCR')
out = p.process(events)
print(out)
p = ZQQProcessor(year='2017',region='VtaggingCR')
out = p.process(events)
print(out)
save(out, 'test.coffea')

