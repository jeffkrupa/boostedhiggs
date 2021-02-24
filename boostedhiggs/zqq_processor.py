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
#from coffea.nanoevents.methods import vector
#import coffea
from boostedhiggs.btag import BTagEfficiency, BTagCorrector
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
    bosonFlavor,
)

#ak.behavior.update(vector.behavior)

def TTsemileptonicmatch(events):

    dataset = events.metadata['dataset']
    if 'TTToSemiLeptonic' not in dataset: 
        return np.zeros(len(events.FatJet.pt),dtype=bool)

    child = events.GenPart[(abs(events.GenPart.pdgId) == 24) & events.GenPart.hasFlags(["isLastCopy","fromHardProcess"])].children
    fatjet = ak.firsts(events.FatJet)
    n_matched_quarks = np.zeros(len(fatjet))

    for ii in [0,1]:
        for jj in [0,1]:
            n_matched_quarks = n_matched_quarks + ak.fill_none( (fatjet.delta_r2(child[:,ii,jj]) < 0.8**2) & (abs(child[:,ii,jj].pdgId) < 6), 0. )
    print(n_matched_quarks) 
    return n_matched_quarks

def VQQgenmatch(events):
   
    dataset = events.metadata['dataset']
    if   'WJetsToQQ'   in dataset: motherId = 24
    elif 'ZJetsToQQ'   in dataset: motherId = 23
    elif 'VectorDiJet' in dataset: motherId = 55
    else: return np.ones(len(events.FatJet.pt),dtype=bool) #events.FatJet.pt.ones_like()

    mother = ak.flatten(events.GenPart[(abs(events.GenPart.pdgId) == motherId) & events.GenPart.hasFlags(["isLastCopy","fromHardProcess"])])
    print(mother)
    print(mother.children.pdgId)
    print(mother.children)
    try: 
       q0 = mother.children[:, 0]
       q1 = mother.children[:, 1]
    except:
       q0 = mother.children[:, 0]
       q1 = mother.children[:, 1]
    print(q0.pdgId,q1.pdgId)
    leading_jet = ak.firsts(events.FatJet)
    return (leading_jet.delta_r2(q0) < 0.8*0.8) & (leading_jet.delta_r2(q1) < 0.8*0.8)

    
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
                #hist.Bin('nJet', r'Number of FatJets', [0.5,1.5,2.5,3.5,4.5,6.5]),#[525,575,625,700,800,1500]),#np.arange(525,2000,50)),
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
                #hist.Bin('genflavor', 'Gen. jet flavor', [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]),
            ),
            'deepAK8': hist.Hist(
                'deepAK8', hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('deepTagMDWqq', 'DeepTagMDWqq', 25,0,1),
                hist.Bin('deepTagMDZqq', 'DeepTAGMDZqq', 25,0,1),
                hist.Bin('msd', r'Jet $m_{sd}$', 41, 40, 200),
                #hist.Bin('genflavor', 'Gen. jet flavor', [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]),
            ),
            'templates': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                #hist.Cat('systematic', 'Systematic'),
                hist.Bin('pt', r'Jet $p_{T}$ [GeV]', 20,200,1000),#[525,575,625,700,800,1500]),#np.arange(525,2000,50)),
                hist.Bin('msd', r'Jet $m_{sd}$', 41, 40, 200),
                hist.Bin('hadW', r'N daughters matched to hadronic W', [-0.5,0.5,1.5,2.5]),#hist.Bin('gruddt', 'GRU$^{DDT}$ value',[-2,0,2]),
                #hist.Bin('n2ddt', 'N$_2^{DDT}$ value', [-2,0,2]),
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
            if 'signal' in region: weight = weights_signal.weight()[cut]
            elif 'muonCR' in region: weight = weights_muonCR.weight()[cut]
            elif 'VtaggingCR' in region: weight = weights_VtaggingCR.weight()[cut]
            output['templates'].fill(
                dataset=dataset,
                region=region,
                pt=normalize(candidatejet.pt, cut),
                msd=normalize(candidatejet.msdcorr, cut),
                #n2ddt=normalize(candidatejet.n2ddt, cut),
                #gruddt=normalize(candidatejet.gruddt, cut),
                in_v3_ddt=normalize(candidatejet.in_v3_ddt, cut),
                hadW=normalize(candidatejet.nmatcheddau,cut),
                weight=weight,
            ),
            output['event'].fill(
                dataset=dataset,
                region=region,
                MET=events.MET.pt[cut],
                #nJet=fatjets.counts[cut],
                nPFConstituents=normalize(candidatejet.nPFConstituents,cut),
                weight=weight,
            ),
            output['deepAK8'].fill(
                dataset=dataset,
                region=region,
                deepTagMDWqq=normalize(candidatejet.deepTagMDWqq,cut),
                deepTagMDZqq=normalize(candidatejet.deepTagMDZqq,cut),
                msd=normalize(candidatejet.msdcorr, cut),
                #genflavor=genflavor[cut],
                weight=weight,
            ),
            output['in_v3'].fill(
                dataset=dataset,
                region=region,
                #genflavor=genflavor[cut], 
                in_v3=normalize(candidatejet.in_v3,cut),
                n2=normalize(candidatejet.n2b1,cut),
                gru=normalize(candidatejet.gru,cut),
                weight=weight,
            ),
            if 'muonCR' in dataset or 'VtaggingCR' in dataset:
                output['muon'].fill(
                    dataset=dataset,
                    region=region,
                    mu_pt=normalize(candidatemuon.pt,cut),
                    mu_eta=normalize(candidatemuon.eta,cut),
                    mu_pfRelIso04_all=normalize(candidatemuon.pfRelIso04_all,cut),
                    weight=weight,
            ),
        #common jet kinematics
        gru = events.GRU
        IN  = events.IN
        fatjets = events.FatJet
        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['qcdrho'] = 2*np.log(fatjets.msdcorr/fatjets.pt)
        fatjets['gruddt'] = gru.v25 - shift(fatjets,algo='gruddt',year=self._year)
        fatjets['gru'] = gru.v25
        fatjets['in_v3'] = IN.v3 
        fatjets['in_v3_ddt'] = IN.v3 - shift(fatjets,algo='inddt',year=self._year)
        fatjets['in_v3_ddt_90pctl'] = IN.v3 - shift(fatjets,algo='inddt90pctl',year=self._year)
        fatjets['n2ddt'] = fatjets.n2b1 - n2ddt_shift(fatjets,year=self._year)
        fatjets['nmatcheddau'] = TTsemileptonicmatch(events)
        dataset = events.metadata['dataset']
        print('process dataset', dataset)
        isRealData = not hasattr(events,'genWeight')
        output = self.accumulator.identity()
        if(len(events) == 0): return output

        selection = PackedSelection('uint64')

        weights_signal = Weights(len(events))
        weights_muonCR = Weights(len(events))
        weights_VtaggingCR = Weights(len(events))

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

            fatjets["genMatchFull"] = VQQgenmatch(events)
            candidatejet = ak.firsts(fatjets)
            candidatejet["genMatchFull"] = VQQgenmatch(events)
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
                & (events.Muon.looseId),
                axis = 1,
            )
            ntaus = ak.sum(
                (events.Tau.pt > 20.)
                & (events.Tau.idDecayMode)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3),
                axis = 1,
            )

            cuts = { "S_fatjet_trigger" : trigger_fatjet,
                     "S_pt" : candidatejet.pt > 525,
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

            fill('signal', cuts.keys())

        #######################
        if 'muonCR' in self._region:
                
            if isRealData:
                trigger_muon = np.zeros(len(events), dtype='bool')
                for t in self._muontriggers[self._year]:
                    trigger_muon = trigger_muon | events.HLT[t]
            else:
                trigger_muon = np.ones(len(events), dtype='bool')
    
    
            candidatejet = ak.firsts(fatjets)
            candidatemuon = events.Muon[:,:5]


            jets = events.Jet[
                ((events.Jet.pt > 50.)
                & (abs(events.Jet.eta) < 2.5)
                & (events.Jet.isTight))
            ][:,:4]

            dphi = abs(jets.delta_phi(candidatejet))

            ak4_away = jets[(dphi > 0.8)]

            nelectrons = ak.sum(
                (events.Electron.pt > 10.)
                & (abs(events.Electron.eta) < 2.5) 
                & (events.Electron.cutBased >= events.Electron.VETO),
                axis = 1,
            )
            nmuons = ak.sum(
                (events.Muon.pt > 10)
                & (abs(events.Muon.eta) < 2.4)
                & (events.Muon.pfRelIso04_all < 0.25)
                & (events.Muon.looseId),
                axis = 1,
            )
            ntaus = ak.sum(
                (events.Tau.pt > 20.)
                & (events.Tau.idDecayMode)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3)
                & (events.Tau.idMVAoldDM2017v1 >= 16),
                axis = 1,
            )


            cuts = { "CR1_muon_trigger"       : trigger_muon,
                     "CR1_jet_pt"             : (candidatejet.pt > 525),
                     "CR1_jet_eta"            : (abs(candidatejet.eta) < 2.5),
                     "CR1_jet_msd"            : (candidatejet.msdcorr > 40),
                     "CR1_jet_rho"            : ((candidatejet.qcdrho > -5.5) & (candidatejet.qcdrho < -2.)),
                     "CR1_mu_pt"              : ak.any(candidatemuon.pt>55,axis=1),
                     "CR1_mu_eta"             : ak.any(abs(candidatemuon.eta)<2.1,axis=1), 
                     "CR1_mu_IDLoose"         : ak.any(candidatemuon.looseId,axis=1),
                     "CR1_mu_isolationTight"  : ak.any(candidatemuon.pfRelIso04_all < 0.15,axis=1),
                     "CR1_muonDphiAK8"        : ak.any(abs(candidatemuon.delta_phi(candidatejet)) > 2*np.pi/3,axis=1),
                     "CR1_ak4btagMedium08"    : (ak.max(ak4_away.btagCSVV2, axis=1, mask_identity=False) > BTagEfficiency.btagWPs[self._year]['medium']), #(ak4_away.btagCSVV2.max() > 0.8838),
                     "CR1_noelectron"         : (nelectrons==0), 
                     "CR1_onemuon"            : (nmuons==1),
                     "CR1_notau"              : (ntaus==0),
                   }
            for name, cut in cuts.items(): 
                selection.add(name, cut)

            if isRealData:
                genflavor = 0 #candidatejet.pt.zeros_like().pad(1, clip=True).fillna(-1).flatten()
            if not isRealData: 
                weights_muonCR.add('genweight', events.genWeight)
                add_pileup_weight(weights_muonCR, events.Pileup.nPU, self._year, dataset)
                #add_singleMuTriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)
                bosons = getBosons(events.GenPart)
                genBosonPt = ak.fill_none(ak.firsts(bosons.pt), 0)
                #add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)  
                #genflavor = matchedBosonFlavor(candidatejet, bosons).pad(1, clip=True).fillna(-1).flatten()

            allcuts_ttbar_muoncontrol = set()
            output['cutflow_muonCR'][dataset]['none']+= float(weights_muonCR.weight().sum())
            for cut in cuts:
                allcuts_ttbar_muoncontrol.add(cut)
                output['cutflow_muonCR'][dataset][cut] += float(weights_muonCR.weight()[selection.all(*allcuts_ttbar_muoncontrol)].sum())
            fill('muonCR', cuts.keys())

        #######################
        if 'VtaggingCR' in self._region:
            if isRealData:
                trigger_muon = np.zeros(len(events), dtype='bool')
                for t in self._muontriggers[self._year]:
                    trigger_muon = trigger_muon | events.HLT[t]
            else:
                trigger_muon = np.ones(len(events), dtype='bool')

            candidatejet = ak.firsts(fatjets)
            candidatemuon = ak.firsts(events.Muon)

 
            jets = events.Jet[
                ((events.Jet.pt > 30.)
                & (abs(events.Jet.eta) < 2.4))
            ][:,:4]

            dr_ak4_ak8 = jets.delta_r(candidatejet)
            dr_ak4_muon = jets.delta_r(candidatemuon)

            ak4_away = jets[(dr_ak4_ak8 > 0.8)]# & (dr_ak4_muon > 0.4)]
            mu_p4 = ak.zip(
                { 
                    "pt"   : ak.fill_none(candidatemuon.pt,0),
                    "eta"  : ak.fill_none(candidatemuon.eta,0),
                    "phi"  : ak.fill_none(candidatemuon.phi,0),
                    "mass" : ak.fill_none(candidatemuon.mass,0),
                },
                with_name="PtEtaPhiMLorentzVector"
            )
    
            met_p4 = ak.zip(
                {
                    "pt"   : ak.from_iter([[v] for v in events.MET.pt]),
                    "eta"  : ak.from_iter([[v] for v in np.zeros(len(events))]),
                    "phi"  : ak.from_iter([[v] for v in events.MET.phi]),
                    "mass" : ak.from_iter([[v] for v in np.zeros(len(events))]),
                },
                with_name="PtEtaPhiMLorentzVector"
            )

            Wleptoniccandidate = mu_p4 + met_p4 

            nelectrons = ak.sum(
                ((events.Electron.pt > 10.)
                & (abs(events.Electron.eta) < 2.5) 
                & (events.Electron.cutBased >= events.Electron.VETO)),
                axis = 1,
            )
            n_tight_muon = ak.sum(
                ((events.Muon.pt>53) 
                & (abs(events.Muon.eta) < 2.1) 
                & (events.Muon.tightId)),
                axis = 1,
            )
            n_loose_muon = ak.sum(
                ((events.Muon.pt>20) 
                & (events.Muon.looseId) 
                & (abs(events.Muon.eta) < 2.4)),
                axis = 1,
            )
            ntaus = ak.sum(
                ((events.Tau.pt > 20.)
                & (events.Tau.idDecayMode)
                & (events.Tau.rawIso < 5)
                & (abs(events.Tau.eta) < 2.3)
                & (events.Tau.idMVAoldDM2017v1 >= 16)),
                axis = 1,   
            )

            cuts = {
                "CR2_muon_trigger"     : trigger_muon,
                "CR2_jet_pt"           : (candidatejet.pt > 200),
                "CR2_jet_eta"          : (abs(candidatejet.eta) < 2.5),
                "CR2_jet_msd"          : (candidatejet.msdcorr > 40),
                "CR2_mu_pt"            : candidatemuon.pt>53,
                "CR2_mu_eta"           : (abs(candidatemuon.eta) < 2.1),
                "CR2_mu_IDTight"       : candidatemuon.tightId,
                "CR2_mu_isolationTight": (candidatemuon.pfRelIso04_all < 0.15),
                "CR2_muonDphiAK8"      : abs(candidatemuon.delta_phi(candidatejet)) > 2*np.pi/3,
                "CR2_ak4btagMedium08"  : (ak.max(ak4_away.btagCSVV2, axis=1, mask_identity=False) > BTagEfficiency.btagWPs[self._year]['medium']),
                "CR2_leptonicW"        : ak.flatten(Wleptoniccandidate.pt>200),
                "CR2_MET"              : (events.MET.pt > 40.),
                "CR2_noelectron"       : (nelectrons==0),
                "CR2_one_tightMuon"    : (n_tight_muon==1),
                "CR2_one_looseMuon"    : (n_loose_muon==1),
                #"CR2_notau"            : (ntaus==0),
            }
    
            for name, cut in cuts.items(): 
                print(name, cut); selection.add(name, cut)
            #weights.add('metfilter', events.Flag.METFilters) 
            if isRealData:
                genflavor = 0 #candidatejet.pt.zeros_like().pad(1, clip=True).fillna(-1).flatten()
            if not isRealData: 
                weights_VtaggingCR.add('genweight', events.genWeight)
                add_pileup_weight(weights_VtaggingCR, events.Pileup.nPU, self._year, dataset)
                #add_singleMuTriggerWeight(weights, abs(candidatemuon.eta), candidatemuon.pt, self._year)
                bosons = getBosons(events.GenPart)
                genBosonPt = ak.fill_none(ak.firsts(bosons.pt), 0)
                #add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)  
                #genflavor = matchedBosonFlavor(candidatejet, bosons).pad(1, clip=True).fillna(-1).flatten()
              
                #b-tag weights
            allcuts_vselection = set()
            output['cutflow_VtaggingCR'][dataset]['none']+= float(weights_VtaggingCR.weight().sum())

            for cut in cuts:
                allcuts_vselection.add(cut)
                output['cutflow_VtaggingCR'][dataset][cut] += float(weights_VtaggingCR.weight()[selection.all(*allcuts_vselection)].sum())
            fill('VtaggingCR',cuts.keys())

        return output

    def postprocess(self, accumulator):
        return accumulator

