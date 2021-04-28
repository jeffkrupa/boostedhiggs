# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
#NanoAODSchema.warn_missing_crossrefs = True
#from ddt_processor import DDTProcessor
from zqq_processor_particleNet import ZQQProcessor
from coffea.util import save, load
#from test import ZQQProcessor
import uproot
filename="root://cmseos.fnal.gov//store/user/lpcpfnano/anovak//v2nano16/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_758fbe82f8c28d1f4ae1816a1c1ab675_0.root"
f = uproot.open(filename)#, xrootd_handler=uproot.MultithreadedXRootDSource)
events = NanoEventsFactory.from_root(
    filename, 
    metadata={'dataset':"ZJetsToQQ_HT-800toInf_qc19_4j_TuneCUETP8M1_13TeV-madgraphMLM-"},#filename.split('/')[11]},#TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8'},#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'},
    schemaclass=NanoAODSchema,
).events()
p = ZQQProcessor(year='2017',region=['signal','muonCR','VtaggingCR'])
out = p.process(events)
print(out)
save(out, 'test.coffea')

