# import mplhep
# plt.style.use(mplhep.style.ATLAS)
#from twoProngGRU import *
#from boostedhiggs import HbbProcessor, BTagEfficiency
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
from hbbprocessor import HbbProcessor
from ddt_processor import DDTProcessor
from zqq_processor import ZQQProcessor
#from test import ZQQProcessor
import uproot

filename="test/nano_mc_2017_1TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"
f = uproot.open(filename)#, xrootd_handler=uproot.MultithreadedXRootDSource)
events = NanoEventsFactory.from_root(
    filename, 
    entry_stop=10000,
    metadata={'dataset':'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'},
    schemaclass=NanoAODSchema,
).events()
print(events.fields)
p = ZQQProcessor(year='2017',region='signal')
out = p.process(events)
print(out)
save(out, 'test.coffea')

