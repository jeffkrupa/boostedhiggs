import logging
from functools import partial
import numpy as np
import awkward
from coffea import processor, hist
np.set_printoptions(threshold=1000)
from uproot_methods import TLorentzVectorArray
from coffea.analysis_tools import Weights, PackedSelection 
from coffea.nanoaod.methods import collection_methods, Candidate
import awkward1 as ak 
import coffea
print(coffea.__version__)
from boostedhiggs.corrections import (
    corrected_msoftdrop,
    gruddt_shift,
    n2ddt_shift,
    shift,
    add_pileup_weight,
    add_jetTriggerWeight,
    add_VJets_NLOkFactor,
    add_singleMuTriggerWeight,
)
from boostedhiggs.common import (
    getBosons,
    matchedBosonFlavor,
)


def VQQgenmatch(events):
   
    dataset = events.metadata['dataset']
    if   'WJetsToQQ'   in dataset: motherId = 24
    elif 'ZJetsToQQ'   in dataset: motherId = 23
    elif 'VectorDiJet' in dataset: motherId = 55
    else: return events.FatJet.pt.ones_like()

    mother = events.GenPart[(abs(events.GenPart.pdgId) == motherId) & events.GenPart.hasFlags(["isLastCopy","fromHardProcess"])].flatten()

    try: 
       q0 = mother.children[:, 0]
       q1 = mother.children[:, 1]
    except:
       q0 = mother.children[:, 0]
       q1 = mother.children[:, 1]

    return (events.FatJet.delta_r2(q0) < 0.8*0.8) & (events.FatJet.delta_r2(q1) < 0.8*0.8)

    
class ZQQProcessor(processor.ProcessorABC):
    def __init__(self, year='2017',region='signal'):
        self._year = year
        self._region = region
        self._triggers = {
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
            '2017': [
                'Mu50', 
                 #'Mu55',
                 #'OldMu100',
                 ]
        }
        self._accumulator = processor.dict_accumulator({
            'sumw': processor.defaultdict_accumulator(float),
            'event': hist.Hist(
                'Event', hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('MET', r'MET [GeV]', 20,20,500),#[525,575,625,700,800,1500]),#np.arange(525,2000,50)),
                hist.Bin('nJet', r'Number of FatJets', [0.5,1.5,2.5,3.5,4.5,6.5]),#[525,575,625,700,800,1500]),#np.arange(525,2000,50)),
                hist.Bin('nPFConstituents', r'Number of PFCandidates', 30,0,60),#[525,575,625,700,800,1500]),#np.arange(525,2000,50)),
            ),
            'muon': hist.Hist(
                'Muon', hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('mu_pt', 'Leading muon p_{T}', 15,50., 700.),
                hist.Bin('mu_pfRelIso04_all', 'Muon pfRelIso04 isolation', 10,0.,0.25),
            ),
            'in_v3': hist.Hist(
                'in_v3', hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('in_v3', 'IN  value', 25,0,1),
                hist.Bin('gru', 'GRU  value', 25,0,1),
                hist.Bin('n2', 'n2  value', 25,0,0.5),
                hist.Bin('genflavor', 'Gen. jet flavor', [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]),
            ),
            'deepAK8': hist.Hist(
                'deepAK8', hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('deepTagMDWqq', 'DeepTagMDWqq', 25,0,1),
                hist.Bin('deepTagMDZqq', 'DeepTAGMDZqq', 25,0,1),
                hist.Bin('msd', r'Jet $m_{sd}$', 23, 40, 300),
                hist.Bin('genflavor', 'Gen. jet flavor', [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]),
            ),
            'templates': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                #hist.Cat('systematic', 'Systematic'),
                hist.Bin('pt', r'Jet $p_{T}$ [GeV]', 25,200,1000),#[525,575,625,700,800,1500]),#np.arange(525,2000,50)),
                hist.Bin('msd', r'Jet $m_{sd}$', 23, 40, 300),
                hist.Bin('gruddt', 'GRU$^{DDT}$ value',[-2,0,2]),
                hist.Bin('n2ddt', 'N$_2^{DDT}$ value', [-2,0,2]),
                hist.Bin('in_v3_ddt', 'IN$^{DDT}$  value', [-2,0,2]),
            ),
            'cutflow_signal' : processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_muonCR' : processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_VtaggingCR' : processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),

        })
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):

        def normalize(val, cut):
            return ak.to_numpy(ak.fill_none(val[cut], np.nan)) #val[cut].pad(1, clip=True).fillna(0).flatten()

        def fill(region, cuts, systematic=None, wmod=None):
            print('filling %s'%region)
            selections = cuts
            cut = selection.all(*selections)
            weight = weights.weight()[cut]
            output['templates'].fill(
                dataset=dataset,
                region=region,
                pt=normalize(candidatejet.pt, cut),
                msd=normalize(candidatejet.msdcorr, cut),
                n2ddt=normalize(candidatejet.n2ddt, cut),
                gruddt=normalize(candidatejet.gruddt, cut),
                in_v3_ddt=normalize(candidatejet.in_v3_ddt_90pctl, cut),
                weight=weight,
            ),
            output['event'].fill(
                dataset=dataset,
                region=region,
                MET=events.MET.pt[cut],
                nJet=fatjets.counts[cut],
                nPFConstituents=normalize(candidatejet.nPFConstituents,cut),
                weight=weight,
            ),
            output['deepAK8'].fill(
                dataset=dataset,
                region=region,
                deepTagMDWqq=normalize(candidatejet.deepTagMDWqq,cut),
                deepTagMDZqq=normalize(candidatejet.deepTagMDZqq,cut),
                msd=normalize(candidatejet.msdcorr, cut),
                genflavor=genflavor[cut],
                weight=weight,
            ),
            output['in_v3'].fill(
                dataset=dataset,
                region=region,
                genflavor=genflavor[cut], 
                in_v3=normalize(candidatejet.in_v3,cut),
                n2=normalize(candidatejet.n2b1,cut),
                gru=normalize(candidatejet.gru,cut),
                weight=weight,
            ),
            if 'muonCR' in dataset:
                output['muon'].fill(
                    dataset=dataset,
                    region=region,
                    mu_pt=normalize(candidatemuon.pt,cut),
                    mu_pfRelIso04_all=normalize(candidatemuon.pfRelIso04_all,cut),
                    weight=weight,
            ),
        #common jet kinematics
        gru = events.GRU
        IN  = events.IN
        fatjets = events.FatJet
        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['rhocorr'] = 2*np.log(fatjets.msdcorr/fatjets.pt)
        fatjets['gruddt'] = gru.v25 - shift(fatjets,algo='gruddt',year=self._year)
        fatjets['gru'] = gru.v25
        fatjets['in_v3'] = IN.v3 
        fatjets['in_v3_ddt'] = IN.v3 - shift(fatjets,algo='inddt',year=self._year)
        fatjets['in_v3_ddt_90pctl'] = IN.v3 - shift(fatjets,algo='inddt90pctl',year=self._year)
        fatjets['n2ddt'] = fatjets.n2b1 - n2ddt_shift(fatjets,year=self._year)

        dataset = events.metadata['dataset']
        print('process dataset', dataset)
        isRealData = not hasattr(events,'genWeight')
        output = self.accumulator.identity()
        if(len(events) == 0): return output

        selection = PackedSelection()
        weights = Weights(len(events))

        if not isRealData:
            output['sumw'][dataset] += ak.sum(events.genWeight)

        #######################
        if 'signal' in self._region:
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

            selection.add('fatjet_trigger', trigger_fatjet)
            fatjets["genMatchFull"] = VQQgenmatch(events)
            candidatejet = fatjets[:,:1]

            nelectrons = ak.sum(
                (events.Electron.pt > 10.)
                & (abs(events.Electron.eta) < 2.5) 
                & (events.Electron.cutBased >= events.Electron.VETO),
                axis = 1,
            )
            nmuons = ak.sum(
                (events.Muon.pt > 10)
                & (abs(events.Muon.eta) < 2.1)
                & (events.Muon.pfRelIso04_all < 0.4)
                & (events.Muon.looseId).astype(bool),
                axis = 1,
            )
            ntaus = ak.sum(
                (events.Tau.pt > 20.)
                & (events.Tau.idDecayMode).astype(bool)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3),
                axis = 1,
            )

            cuts = { "pt" : (candidatejet.pt > 525),
                     "eta" : (abs(candidatejet.eta) < 2.5),
                     "msdcorr" : (candidatejet.msdcorr > 40),
                     "rho"     : ((candidatejet.rhocorr > -5.5) &(candidatejet.rhocorr < -2.)),
                     "jetid"   : (candidatejet.isTight), 
                     "VQQgenmatch" : (fatjets.genMatchFull), 
                     "noelectron" : (nelectrons == 0),
                     "nomuon"     : (nmuons == 0),
                     "notau"      : (ntaus == 0),
                   }

            for name, cut in cuts.items():
                selection.add(name, cut)

            if isRealData:
                genflavor = 0 #candidatejet.pt.zeros_like().pad(1, clip=True).fillna(-1).flatten()
            if not isRealData: 
                weights.add('genweight', events.genWeight)
                add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
                add_jetTriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)
                bosons = getBosons(events.GenPart)
                genBosonPt = ak.fill_none(ak.firsts(bosons.pt), 0)
                add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)  
                genflavor = matchedBosonFlavor(candidatejet, bosons).pad(1, clip=True).fillna(-1).flatten()

            allcuts_signal = set()
            output['cutflow_signal'][dataset]['none']+= float(weights.weight().sum())
            for cut in cuts:
                allcuts_signal.add(cut)
                output['cutflow_signal'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_signal)].sum())

            fill('signal', cuts.keys())

        #######################
        if 'muonCR' in self._region:
                
            if isRealData:
                trigger_muon = np.zeros(events.size, dtype='bool')
                for t in self._muontriggers[self._year]:
                    trigger_muon = trigger_muon | events.HLT[t]
            else:
                trigger_muon = np.ones(events.size, dtype='bool')
    
    
            candidatejet = fatjets[:,:1]
            candidatemuon = events.Muon[:,:5]

            muon_ak8_pair = candidatemuon.cross(candidatejet,nested=True) 

            jets = events.Jet[
                ((events.Jet.pt > 50.)
                & (abs(events.Jet.eta) < 2.5))
            ][:,:4]

            ak4_ak8_pair = jets.cross(candidatejet, nested=True)
            dr = abs(ak4_ak8_pair.i0.delta_r(ak4_ak8_pair.i1))
            ak4_away = jets[(dr > 0.8).all()]

            nelectrons = (
                ((events.Electron.pt > 10.)
                & (abs(events.Electron.eta) < 2.5) 
                & (events.Electron.cutBased >= events.Electron.VETO))
            ).sum()
            ntaus = (
                ((events.Tau.pt > 20.)
                & (events.Tau.idDecayMode).astype(bool)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3))
            ).sum()

            cuts = { "muon_trigger"       : trigger_muon,
                     "jet_pt"             : (candidatejet.pt > 525).any(),
                     "jet_eta"            : (abs(candidatejet.eta) < 2.5).any(),
                     "jet_rho"            : ((candidatejet.rhocorr > -5.5) & (candidatejet.rhocorr < -2)).any(),
                     "jet_msd"            : (candidatejet.msdcorr > 40).any(),
                     "mu_pt"              : (candidatemuon.pt>55).any(),
                     "mu_eta"             : (abs(candidatemuon.eta)<2.1).any(), 
                     "mu_IDLoose"         : (candidatemuon.looseId).astype(bool).any(),
                     "mu_isolationTight"  : (candidatemuon.pfRelIso04_all < 0.15).any(),
                     "muonDphiAK8"        : (abs(muon_ak8_pair.i0.delta_phi(muon_ak8_pair.i1)) > 2*np.pi/3).all().all(),
                     "ak4btagMedium08"    : (ak4_away.btagCSVV2.max() > 0.8838),
                     "noelectron"         : (nelectrons==0), 
                     "notau"              : (ntaus==0),
                   }
            for name, cut in cuts.items(): 
                selection.add(name, cut)

            if isRealData:
                genflavor = candidatejet.pt.zeros_like().pad(1, clip=True).fillna(-1).flatten()
            if not isRealData: 
                weights.add('genweight', events.genWeight)
                add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
                add_singleMuTriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)
                bosons = getBosons(events)
                genBosonPt = bosons.pt.pad(1, clip=True).fillna(0)
                add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)  
                genflavor = matchedBosonFlavor(candidatejet, bosons).pad(1, clip=True).fillna(-1).flatten()

            allcuts_ttbar_muoncontrol = set()
            output['cutflow_muonCR'][dataset]['none']+= float(weights.weight().sum())
            for cut in cuts:
                allcuts_ttbar_muoncontrol.add(cut)
                output['cutflow_muonCR'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_ttbar_muoncontrol)].sum())
            fill('muonCR', cuts.keys())

        #######################
        if 'VtaggingCR' in self._region:
            if isRealData:
                trigger_muon = np.zeros(events.size, dtype='bool')
                for t in self._muontriggers[self._year]:
                    trigger_muon = trigger_muon | events.HLT[t]
            else:
                trigger_muon = np.ones(events.size, dtype='bool')

            candidatejet = fatjets[:,:1]
            candidatemuon = events.Muon[:,:1]

 
            jets = events.Jet[
                ((events.Jet.pt > 30.)
                & (abs(events.Jet.eta) < 2.4))
            ][:,:4]
            muon_ak8_pair = candidatemuon.cross(candidatejet,nested=True) 
            muon_ak4_pair = jets.cross(candidatemuon,nested=True)
 
            ak4_ak8_pair = jets.cross(candidatejet, nested=True)
            dr_ak4_ak8 = abs(ak4_ak8_pair.i0.delta_r(ak4_ak8_pair.i1))
            dr_muon_ak4 = abs(muon_ak4_pair.i0.delta_r(muon_ak4_pair.i1))
            ak4_away = jets[((dr_ak4_ak8 > 0.8).all() & (dr_muon_ak4 > 0.3).all())]


            mu_p4 = TLorentzVectorArray.from_ptetaphim(candidatemuon.pt.fillna(0),candidatemuon.eta.fillna(0),candidatemuon.phi.fillna(0),candidatemuon.mass.fillna(0))
            met_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in events.MET.pt]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]), awkward.JaggedArray.fromiter([[v] for v in events.MET.phi]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]))
            met_candidatemuon_pair = met_p4.cross(mu_p4)
            Wleptoniccandidate = met_candidatemuon_pair.i0 + met_candidatemuon_pair.i1

            nelectrons = (
                ((events.Electron.pt > 10.)
                & (abs(events.Electron.eta) < 2.5) 
                & (events.Electron.cutBased >= events.Electron.VETO))
            ).sum()
            n_tight_muon = ((events.Muon.pt>53) & (abs(events.Muon.eta) < 2.1) & (events.Muon.tightId).astype(bool))
            n_loose_muon = ((events.Muon.pt>20) & (events.Muon.looseId).astype(bool) & (abs(events.Muon.eta) < 2.4))
            ntaus = (
                ((events.Tau.pt > 20.)
                & (events.Tau.idDecayMode).astype(bool)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3))
            ).sum()

            cuts = { "muon_trigger"        : trigger_muon,
                     "jet_pt"              : (candidatejet.pt > 200).any(),
                     "jet_eta"             : (abs(candidatejet.eta) < 2.5).any(),
                     "jet_msd"             : (candidatejet.msdcorr > 40).any(),
                     "mu_pt"               : (candidatemuon.pt>53).any(),
                     "mu_eta"              : (abs(candidatemuon.eta)<2.1).any(), 
                     "mu_IDTight"          : (candidatemuon.tightId).astype(bool).any(),
                     #"mu_isolationTight"   : (candidatemuon.pfRelIso04_all < 0.15).any(),
                     "one_tightMuon"       : (n_tight_muon.sum()==1),
                     "one_looseMuon"       : (n_loose_muon.sum()==1),
                     "muonDphiAK8"         : (abs(muon_ak8_pair.i0.delta_phi(muon_ak8_pair.i1)) > 2*np.pi/3).all().all(),
                     "ak4btagMedium08"     : (ak4_away.btagCSVV2.max() > 0.8838),
                     "leptonicW"           : (Wleptoniccandidate.pt>200).any(),
                     "MET"                 : (events.MET.pt > 40.),
                     "noelectron"          : (nelectrons==0), 
                     "notau"               : (ntaus==0),
                   }
            for name, cut in cuts.items(): 
                selection.add(name, cut)
    
            #weights.add('metfilter', events.Flag.METFilters) 
            if isRealData:
                genflavor = candidatejet.pt.zeros_like().pad(1, clip=True).fillna(-1).flatten()
            if not isRealData: 
                weights.add('genweight', events.genWeight)
                add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
                add_singleMuTriggerWeight(weights, abs(candidatemuon.eta), candidatemuon.pt, self._year)
                bosons = getBosons(events)
                genBosonPt = bosons.pt.pad(1, clip=True).fillna(0)
                add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)  
                genflavor = matchedBosonFlavor(candidatejet, bosons).pad(1, clip=True).fillna(-1).flatten()
              
                #b-tag weights
            allcuts_vselection = set()
            output['cutflow_VtaggingCR'][dataset]['none']+= float(weights.weight().sum())

            for cut in cuts:
                allcuts_vselection.add(cut)
                output['cutflow_VtaggingCR'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_vselection)].sum())
            fill('VtaggingCR',cuts.keys())

        return output

    def postprocess(self, accumulator):
        return accumulator

