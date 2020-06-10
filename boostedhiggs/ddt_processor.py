import logging
import numpy as np
from coffea import processor, hist
from boostedhiggs.twoProngGRU import *
from coffea.nanoaod.methods import collection_methods, Candidate
collection_methods["FatJetPFCands"] = Candidate

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
            'jet_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 100, 200, 1200),
                hist.Bin('jet_rho', r'Jet $\rho$', 180, -7., -1.5),
                hist.Bin('jet_twoProngGru', r'Jet GRU score', 100, 0., 1.0),
            ),
            'cutflow': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('cut', 'Cut index', 11, 0, 11),
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):

        dataset = events.metadata['dataset']
        isRealData = False
        selection = processor.PackedSelection()
        weights = processor.Weights(len(events))
        output = self.accumulator.identity()
        output['sumw'][dataset] += events.genWeight.sum()

        trigger_fatjet = np.zeros(events.size, dtype='bool')

        # trigger paths
        for t in self._triggers[self._year]:
            trigger_fatjet = trigger_fatjet | events.HLT[t]
        selection.add('fatjet_trigger', trigger_fatjet)

        # run model on PFCands associated to FatJet (FatJetPFCands)
        events.FatJet.array.content["PFCands"] = type(events.FatJetPFCands.array).fromcounts(events.FatJet.nPFConstituents.flatten(), events.FatJetPFCands.flatten())
        events.FatJet.array.content["twoProngGru"] = run_model(events.FatJet.flatten())


        fatjets = events.FatJet
        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        candidatejet = fatjets[
            # https://github.com/DAZSLE/BaconAnalyzer/blob/master/Analyzer/src/VJetLoader.cc#L269
            (fatjets.pt > 250)
            & (abs(fatjets.eta) < 2.5)
            # & fatjets.isLoose  # not always available

        ][:, 0:1]

        #events.FatJet = candidatejet


        # basic jet selection
        selection.add('minjetkin', (
            (candidatejet.pt >= 450)
            & (candidatejet.msdcorr >= 40.)
            & (abs(candidatejet.eta) < 2.5)
        ).any())

        # lep veto
        nmuons = (
            (events.Muon.pt > 10)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.pfRelIso04_all < 0.25)
            & (events.Muon.looseId).astype(bool)
        ).sum()

        nelectrons = (
            (events.Electron.pt > 10)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased >= events.Electron.LOOSE)
        ).sum()

        ntaus = (
            (events.Tau.pt > 20)
            & (events.Tau.idDecayMode).astype(bool)
            # bacon iso looser than Nano selection
        ).sum()
        selection.add('jetid', candidatejet.isTight.any())

        selection.add('noleptons', (nmuons == 0) & (nelectrons == 0) & (ntaus == 0))
        #weights.add('genweight', events.genWeight)
        #add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
        #add_jetTriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)

        regions = {
           'signal' : ['fatjet_trigger', 'minjetkin','noleptons','jetid']
        }
        for region, cuts in regions.items():
            allcuts = set()
            output['cutflow'].fill(dataset=dataset, region=region, cut=0)
            for i, cut in enumerate(cuts):
                allcuts.add(cut)
                cut = selection.all(*allcuts)
                output['cutflow'].fill(dataset=dataset, region=region, cut=i + 1)# weight=weights.weight()[cut])

        def normalize(val, cut):
            return val[cut].pad(1, clip=True).fillna(0).flatten()

        def fill(region, systematic=None, wmod=None):
            selections = regions[region]
            cut = selection.all(*selections)
            sname = 'nominal' if systematic is None else systematic
            weight = weights.weight()[cut]
            output['jet_kin'].fill(
                dataset=dataset,
                region=region,
                jet_pt=normalize(candidatejet.pt, cut),
                jet_rho=normalize(2*np.log(candidatejet.msdcorr/candidatejet.pt), cut),
                jet_twoProngGru=normalize(candidatejet.twoProngGru, cut),
            )

        for region in regions:
            fill(region)

        return output

    def postprocess(self, accumulator):
        return accumulator

