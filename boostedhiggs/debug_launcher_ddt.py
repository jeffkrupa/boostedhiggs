# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
NanoAODSchema.warn_missing_crossrefs = True
#from ddt_processor import DDTProcessor
from ddt_processor import DDTProcessor
from coffea.util import save, load
#from test import ZQQProcessor
import uproot

filename="/uscms/home/jkrupa/nobackup/zprlegacy/CMSSW_11_2_1/src/PhysicsTools/NanoAODTools/tmp/nano_mc_2017_99test.root"
f = uproot.open(filename)#, xrootd_handler=uproot.MultithreadedXRootDSource)
events = NanoEventsFactory.from_root(
    filename, 
    metadata={'dataset':'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'},
    schemaclass=NanoAODSchema,
).events()
p = DDTProcessor(year='2017')
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

