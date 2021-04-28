# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
NanoAODSchema.warn_missing_crossrefs = True
#from ddt_processor import DDTProcessor
from zqq_processor_Yihan import ZQQProcessor
from coffea.util import save, load
#from test import ZQQProcessor
import uproot
#filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2017_preUL/14Apr21_preUL/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/nano_mc_2017_9_post-process_pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1_hadd.root"
#filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2017_preUL/14Apr21_preUL/SingleMuon_pancakes-02_Run2017B-31Mar2018-v1/nano_data_2017_9_post-process__hadd.root"
#filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2016_preUL/14Apr21_preUL/WJetsToQQ_HT-800toInf_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/nano_mc_2016_9_post-process_pancakes-02_RunIISummer16MiniAODv3-PUMoriond17_94X_v3-v1_hadd.root"
#filename="/uscms/home/jkrupa/nobackup/zprlegacy/CMSSW_11_2_1/src/PhysicsTools/NanoAODTools/tmp/nano_mc_2017_99test.root"
filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2016_preUL/14Apr21_preUL/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/nano_mc_2016_9_post-process_pancakes-02_RunIISummer16MiniAODv3-PUMoriond17_94X_v3-v1_hadd.root"
#filename="root://cmseos.fnal.gov//store/user/lpcpfnano/jkrupa/nanopost_process/2016_preUL/14Apr21_preUL/JetHT_pancakes-02_Run2016B-17Jul2018_ver2-v2/nano_data_2016_9_post-process__hadd.root"
f = uproot.open(filename)#, xrootd_handler=uproot.MultithreadedXRootDSource)
events = NanoEventsFactory.from_root(
    filename, 
    metadata={'dataset':filename.split('/')[11]},#TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8'},#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'},
    schemaclass=NanoAODSchema,
).events()
p = ZQQProcessor(year='2017',region=['signal','muonCR','VtaggingCR'])
out = p.process(events)
print(out)
save(out, 'test.coffea')

