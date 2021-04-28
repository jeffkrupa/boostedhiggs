import logging
from functools import partial
import numpy as np
from coffea import processor, hist
np.set_printoptions(threshold=1000)
#from uproot3_methods import TLorentzVectorArray
from coffea.analysis_tools import Weights, PackedSelection 
#from coffea.nanoaod.methods import collection_methods, Candidate
import coffea
print('coffea==',coffea.__version__)
import awkward as ak 
print('awkward==',ak.__version__)
import uproot
print('uproot==', uproot.__version__)
import cloudpickle
print('cloudpickle==', cloudpickle.__version__)
#from coffea.nanoevents.methods import vector
#import coffea
from boostedhiggs.btag import BTagEfficiency, BTagCorrector
from boostedhiggs.corrections_Yihan import (
    corrected_msoftdrop,
    gruddt_shift,
    n2ddt_shift,
    shift,
    add_pileup_weight,
    add_jetTriggerWeight,
    add_VJets_NLOkFactor,
    add_singleMuTriggerWeight,
    jet_factory,
    fatjet_factory,
    add_jec_variables,
    met_factory,
)
from boostedhiggs.common import (
    getBosons,
    bosonFlavor,
)

#ak.behavior.update(vector.behavior)


def VQQgenmatch(events):
   
    dataset = events.metadata['dataset']
    if   'WJetsToQQ'   in dataset: motherId = 24
    elif 'ZJetsToQQ'   in dataset: motherId = 23
    elif 'VectorDiJet' in dataset: motherId = 55
    else: return np.ones(len(events.FatJet.pt),dtype=bool), np.ones(len(events.FatJet.pt),dtype=int) #events.FatJet.pt.ones_like()

    mother = ak.flatten(events.GenPart[(abs(events.GenPart.pdgId) == motherId) & events.GenPart.hasFlags(["isLastCopy","fromHardProcess"])])
    #print(mother)
    #print(mother.children.pdgId)
    #print(mother.children)
    try: 
       q0 = mother.children[:, 0]
       q1 = mother.children[:, 1]
    except:
       q0 = mother.children[:, 0]
       q1 = mother.children[:, 1]
    #print(q0.pdgId,q1.pdgId)
    childid = abs(mother.children.pdgId)
    genflavor = ak.any(childid == 5, axis=-1) * 3 + ak.any(childid == 4, axis=-1) * 2 + ak.any(childid < 4, axis=-1) * 1
    print(childid[genflavor>3])
    print(len(childid[genflavor>3]))
    print(genflavor[genflavor>3])
    leading_jet = ak.firsts(events.FatJet)
    return (leading_jet.delta_r2(q0) < 0.8*0.8) & (leading_jet.delta_r2(q1) < 0.8*0.8), genflavor

def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
    return out    

class ZQQProcessor(processor.ProcessorABC):
    def __init__(self, year='2017',region='signal'):
        self._year = year
        self._region = region

        self._fjpTcut = {"2016" : 500.,
                         "2017" : 525.,
                         "2018" : 500.,
                        }
        self._btagSF = BTagCorrector(year, 'medium')
        self._triggers = {
            '2016': [
                'AK8PFJet360_TrimMass30',
		'AK8PFHT700TrimR0p1PT0p03Mass50',
                'PFJet500',
                'AK8PFJet500', 
                'PFHT800',
                'PFHT900',
                'PFHT650_WideJetMJJ900DEtaJJ1p5',
                'PFHT650_WideJetMJJ950DEtaJJ1p5',
                'AK8DiPFJet300_200_TrimMass30_BTagCSV_p20'
                 ],
            '2017': [
                'PFHT1050',
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFJet500',
                'AK8PFJet500',
                'AK8PFJet550',
                'CaloJet500_NoJetID',
                'CaloJet550_NoJetID', 
                 ]
        }
        self._muontriggers = {
            '2016': [
                'Mu50', 
                 'Mu55',
                 ],
            '2017': [
                'Mu50', 
                 #'Mu55',
                 #'OldMu100',
                 ]
        }
        self._accumulator = processor.dict_accumulator({
            'sumw': processor.defaultdict_accumulator(float),
            'particleNet': hist.Hist(
                'particleNet', hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('msd', 'msd', 25,40,350),
                hist.Bin('genflavor', 'Gen. jet flavor', [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]),
                hist.Bin('particleNetXqqInclusive', 'IN_Sep20_2017 value',  20,0,1),
            ),
            'cutflow' : processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
        
        })
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        print ('###########running test processor!#############')
        dataset = events.metadata['dataset']
        print('process dataset', dataset)
        isRealData = not hasattr(events,'genWeight')
        output = self.accumulator.identity()
             
        if isRealData:
            # Nominal JEC are already applied in data
            output += self.process_shift(events, None)
            return output

        jec_cache = {}
        fatjets = fatjet_factory[f"{self._year}mc"].build(add_jec_variables(events.FatJet,events.FatJet.nearest(events.GenJetAK8), events.fixedGridRhoFastjetAll), jec_cache)
        jets = jet_factory[f"{self._year}mc"].build(add_jec_variables(events.Jet, events.Jet.nearest(events.GenJet), events.fixedGridRhoFastjetAll), jec_cache)
        met = met_factory.build(events.MET, jets, {})

        output += self.process_shift(update(events, {"Jet": jets, "FatJet": fatjets, "MET": met}), None)
        
        return output
    
    def process_shift(self, events, shift_name):

        if shift_name is None:
            systematic = 'nominal'
        else:
            systematic = shift_name
        
        # print ('#####shift name #####')
        # print (systematic)

        def normalize(val, cut):
            return ak.to_numpy(ak.fill_none(val[cut], np.nan)) #val[cut].pad(1, clip=True).fillna(0).flatten()

        def fill(region, cuts, systematic, wmod=None):
            print('filling %s'%region)
            selections = cuts

            sname = systematic

            cut = selection.all(*selections)
            weight = weights.weight()[cut]
            output['particleNet'].fill(
                dataset=dataset,
                region=region,
                genflavor=normalize(candidatejet.genflavor,cut), 
                msd=normalize(candidatejet.msdcorr, cut),
                particleNetXqqInclusive=normalize(candidatejet.particleNetXqqInclusive, cut),
                weight=weight,
            ),
        #common jet kinematics
        #gru = events.GRU
        IN  = events.IN
        fatjets = events.FatJet
        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['qcdrho'] = 2*np.log(fatjets.msdcorr/fatjets.pt)
        dataset = events.metadata['dataset']
        print('process dataset', dataset)
        isRealData = not hasattr(events,'genWeight')
        output = self.accumulator.identity()
        if(len(events) == 0): return output

        selection = PackedSelection('uint64')

        weights_signal = Weights(len(events))
        ##
        if shift_name is None and not isRealData:
            output['sumw'][dataset] += ak.sum(events.genWeight)

        #######################
        if 1: 
            if isRealData:
                trigger_fatjet = np.zeros(len(events), dtype='bool')
                for t in self._triggers[self._year]:
                    try:
                        trigger_fatjet = trigger_fatjet | events.HLT[t]
                    except:
                        print('trigger %s not available'%t)
                        continue

            else:
                trigger_fatjet = np.ones(len(events), dtype='bool')

            fatjets["genMatchFull"], fatjets["genflavor"] = VQQgenmatch(events)

            fatjets["particleNetXqqInclusive"] = ( fatjets.particleNetMD_Xqq + fatjets.particleNetMD_Xcc + fatjets.particleNetMD_Xbb ) /  ( fatjets.particleNetMD_Xqq + fatjets.particleNetMD_Xcc + fatjets.particleNetMD_Xbb + fatjets.particleNetMD_QCD )


            candidatejet = ak.firsts(fatjets)


            goodelectron=(
                (events.Electron.pt > 10.)
                & (abs(events.Electron.eta) < 2.5) 
                & (events.Electron.cutBased >= events.Electron.LOOSE)
                )
            nelectrons = ak.sum(goodelectron, axis=1)

            goodmuon=(
                (events.Muon.pt > 10)
                & (abs(events.Muon.eta) < 2.1)
                & (events.Muon.pfRelIso04_all < 0.4)
                & (events.Muon.looseId)
                )
            nmuons = ak.sum(goodmuon, axis=1)

            ntaus = ak.sum(
                ((events.Tau.pt > 20.)
                & (events.Tau.idDecayMode)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3)
                & ak.all(events.Tau.metric_table(events.Muon[goodmuon]) > 0.4, axis=2)
                & ak.all(events.Tau.metric_table(events.Electron[goodelectron]) > 0.4, axis=2)
                ), 
            axis = 1,
            )

            cuts = { "S_fatjet_trigger" : trigger_fatjet,
                     "S_pt" : candidatejet.pt > self._fjpTcut[self._year],
                     "S_eta" : (abs(candidatejet.eta) < 2.5),
                     "S_msdcorr" : (candidatejet.msdcorr > 40),
                     "S_rho"     : ((candidatejet.qcdrho > -5.5) &(candidatejet.qcdrho < -2.)),
                     "S_jetid"   : (candidatejet.isTight), 
                     "S_VQQgenmatch" : (candidatejet.genMatchFull), 
                     "S_noelectron" : (nelectrons == 0),
                     "S_nomuon"     : (nmuons == 0),
                     "S_notau"      : (ntaus == 0),
                   }

            for name, cut in cuts.items():
                print(name, cut); selection.add(name, cut)

            if isRealData:
                genflavor = 0 #candidatejet.pt.zeros_like().pad(1, clip=True).fillna(-1).flatten()
            if not isRealData: 
                weights_signal.add('genweight', events.genWeight)
                add_pileup_weight(weights_signal, events.Pileup.nPU, self._year, dataset)
                add_jetTriggerWeight(weights_signal, candidatejet.msdcorr, candidatejet.pt, self._year)
                bosons = getBosons(events.GenPart)
                genBosonPt = ak.fill_none(ak.firsts(bosons.pt), 0)
                add_VJets_NLOkFactor(weights_signal, genBosonPt, self._year, dataset)  
                #genflavor = matchedBosonFlavor(candidatejet, bosons).pad(1, clip=True).fillna(-1).flatten()

            allcuts_signal = set()
            output['cutflow_signal'][dataset]['none']+= float(weights_signal.weight().sum())
            for cut in cuts:
                allcuts_signal.add(cut)
                output['cutflow_signal'][dataset][cut] += float(weights_signal.weight()[selection.all(*allcuts_signal)].sum())
            ##
            fill('signal', cuts.keys(), systematic)


        return output

    def postprocess(self, accumulator):
        return accumulator
 
