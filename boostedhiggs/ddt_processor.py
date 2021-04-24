import logging
import numpy as np
from coffea import processor, hist
#from boostedhiggs.twoProngGRU import *
#from coffea.nanoaod.methods import collection_methods, Candidate
from coffea.analysis_tools import Weights, PackedSelection 
import awkward as ak
import coffea
from functools import partial
from boostedhiggs.corrections import (
    corrected_msoftdrop,
    add_pileup_weight,
    add_jetTriggerWeight,
)

class DDTProcessor(processor.ProcessorABC):
    def __init__(self, year='2017'):
        self._year = year
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
        self._accumulator = processor.dict_accumulator({
            'sumw': processor.defaultdict_accumulator(float),
            'IN_Sep20_2017': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 100, 200, 1500),
                hist.Bin('jet_rho', r'Jet $\rho$', 180, -5.5, -2.),
                hist.Bin('jet_IN_Sep20_2017', 'IN  value', 100, 0.0, 1.0),
            ),
            'IN_Apr21_2017_late': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 100, 200, 1500),
                hist.Bin('jet_rho', r'Jet $\rho$', 180, -5.5, -2.),
                hist.Bin('jet_IN_Apr21_2017_late', 'IN  value', 100, 0.0, 1.0),
            ),
            'IN_Apr21_2017_early': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 100, 200, 1500),
                hist.Bin('jet_rho', r'Jet $\rho$', 180, -5.5, -2.),
                hist.Bin('jet_IN_Apr21_2017_early', 'IN  value', 100, 0.0, 1.0),
            ),
            'cutflow': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('cut', 'Cut index', 11, 0, 11),
            ),
            'cutflow' : processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        def normalize(val, cut):
            return ak.to_numpy(ak.fill_none(val[cut], np.nan)) #val[cut].pad(1, clip=True).fillna(0).flatten()

        def fill(region, cuts, systematic=None, wmod=None):
            selections = cuts
            cut = selection.all(*selections)
            sname = 'nominal' if systematic is None else systematic
            weight = weights.weight()[cut]
            output['IN_Sep20_2017'].fill(
                dataset=dataset,
                region=region,
                jet_pt=normalize(candidatejet.pt, cut),
                jet_rho=normalize(2*np.log(candidatejet.msdcorr/candidatejet.pt), cut),
                jet_IN_Sep20_2017=normalize(candidatejet.IN_Sep20_2017, cut),
            )
            output['IN_Apr21_2017_early'].fill(
                dataset=dataset,
                region=region,
                jet_pt=normalize(candidatejet.pt, cut),
                jet_rho=normalize(2*np.log(candidatejet.msdcorr/candidatejet.pt), cut),
                jet_IN_Apr21_2017_early=normalize(candidatejet.IN_Apr21_2017_early, cut),
            )
            output['IN_Apr21_2017_late'].fill(
                dataset=dataset,
                region=region,
                jet_pt=normalize(candidatejet.pt, cut),
                jet_rho=normalize(2*np.log(candidatejet.msdcorr/candidatejet.pt), cut),
                jet_IN_Apr21_2017_late=normalize(candidatejet.IN_Apr21_2017_late, cut),
            )


        if(len(events) == 0): return output

        dataset = events.metadata['dataset']
        isRealData = False
        selection = PackedSelection()
        weights = Weights(len(events))
        output = self.accumulator.identity()
        output['sumw'][dataset] += ak.sum(events.genWeight)

        fatjets = events.FatJet
        IN  = events.IN

        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['rhocorr'] = 2*np.log(fatjets.msdcorr/fatjets.pt)
        fatjets['IN_Sep20_2017'] = IN.Sep20_2017
        fatjets['IN_Apr21_2017_late']  = IN.Apr21_2017_late
        fatjets['IN_Apr21_2017_early'] = IN.Apr21_2017_early


        candidatejet = ak.firsts(fatjets)

        nmuons = ak.sum(
            (events.Muon.pt > 10)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.pfRelIso04_all < 0.25)
            & (events.Muon.looseId),
            axis=1
        )

        nelectrons = ak.sum(
            (events.Electron.pt > 10)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased >= events.Electron.LOOSE),
            axis=1
        )

        ntaus = ak.sum(
            (events.Tau.pt > 20)
            & (events.Tau.idDecayMode)
            & (events.Tau.rawIso < 5)
            & (abs(events.Tau.eta) < 2.3),
            axis=1
            # bacon iso looser than Nano selection
        )
        add_jetTriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)
        weights.add('genweight', events.genWeight)
        cuts = { 
                     "S_pt" : candidatejet.pt > 200,
                     "S_eta" : (abs(candidatejet.eta) < 2.5),
                     "S_msdcorr" : (candidatejet.msdcorr > 40),
                     "S_rho"     : ((candidatejet.rhocorr > -5.5) &(candidatejet.rhocorr < -2.)),
                     "S_jetid"   : (candidatejet.isTight),
                     "S_noelectron" : (nelectrons == 0),
                     "S_nomuon"     : (nmuons == 0),
                     "S_notau"      : (ntaus == 0),
                   }



        for name, cut in cuts.items():
                print(name, cut); selection.add(name, cut)
        fill('signal',cuts.keys())
        allcuts = set()
        output['cutflow'][dataset]['none']+= float(weights.weight().sum())
        for cut in cuts:
            allcuts.add(cut)
            output['cutflow'][dataset][cut] += float(weights.weight()[selection.all(*allcuts)].sum())

        def normalize(val, cut):
            return val[cut].pad(1, clip=True).fillna(0).flatten()

        return output


    def postprocess(self, accumulator):
        return accumulator

