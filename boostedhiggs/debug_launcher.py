# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
#from ddt_processor import DDTProcessor
from zqq_processor import ZQQProcessor
from coffea.util import save, load
#from test import ZQQProcessor
import uproot

filename="test/nano_mc_2017_1TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"
f = uproot.open(filename)#, xrootd_handler=uproot.MultithreadedXRootDSource)
events = NanoEventsFactory.from_root(
    filename, 
    #entry_stop=100000,
    metadata={'dataset':'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'},
    schemaclass=NanoAODSchema,
).events()
p = ZQQProcessor(year='2017',region=['signal','muonCR','VtaggingCR'])
out = p.process(events)
print(out)
p = ZQQProcessor(year='2017',region='muonCR')
out = p.process(events)
print(out)
p = ZQQProcessor(year='2017',region='VtaggingCR')
out = p.process(events)
print(out)
save(out, 'test.coffea')

