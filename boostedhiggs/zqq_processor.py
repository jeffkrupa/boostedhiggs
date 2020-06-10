import logging
import numpy as np
from coffea import processor, hist
from boostedhiggs.twoProngGRU import *
from coffea.nanoaod.methods import collection_methods, Candidate
collection_methods["FatJetPFCands"] = Candidate

from boostedhiggs.corrections import (
    corrected_msoftdrop,
    gruddt_shift,
    n2ddt_shift,
    add_pileup_weight,
    add_jetTriggerWeight,
)

class ZQQProcessor(processor.ProcessorABC):
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
            'templates': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                #hist.Cat('systematic', 'Systematic'),
                hist.Bin('pt', r'Jet $p_{T}$ [GeV]', [525,575,625,700,800,1500]),
                hist.Bin('msd', r'Jet $m_{sd}$', 62, 40, 350),
                hist.Bin('gruddt', 'GRU$^{DDT}$ value',100,-1.,1.),
                hist.Bin('gru', 'GRU value',100,0.,1.),
                hist.Bin('rho', 'jet rho', [-6.,-5.5,-5.,-4.5,-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.]),
                #hist.Bin('n2ddt', 'N$_2^{DDT}$ value', 50, -0.4, 0.4),
                #hist.Bin('gru','GRU value',80,0.,1.),
            ),
            #'gruddt' : hist.Hist(
            #    hist.Cat('dataset', 'Dataset'),
            #    hist.Cat('region', 'Region'),
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
        print (dataset)
        if not isRealData:
            output['sumw'][dataset] += events.genWeight.sum()



        trigger_fatjet = np.zeros(events.size, dtype='bool')

        # trigger paths
        for t in self._triggers[self._year]:
            trigger_fatjet = trigger_fatjet | events.HLT[t]
        selection.add('fatjet_trigger', trigger_fatjet)

        # run model on PFCands associated to FatJet (FatJetPFCands)
        events.FatJet.array.content["PFCands"] = type(events.FatJetPFCands.array).fromcounts(events.FatJet.nPFConstituents.flatten(), events.FatJetPFCands.flatten())
        events.FatJet.array.content["twoProngGru"] = run_model(events.FatJet.flatten())
   
        #else:
        #  events.FatJet["genMatchFull"] = np.ones(len(events))
        fatjets = events.FatJet

        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['rho'] = 2*np.log(fatjets.msdcorr/fatjets.pt)

        fatjets['gruddt'] = fatjets.twoProngGru - gruddt_shift(fatjets,year=self._year)
        fatjets['n2ddt'] = fatjets.n2b1 - n2ddt_shift(fatjets,year=self._year)
      
        print(type(dataset))
        if 'QCD' not in dataset: 
          zprime = events.GenPart[((events.GenPart.pdgId==55) | (events.GenPart.pdgId==23) | (abs(events.GenPart.pdgId)==24)) & events.GenPart.hasFlags(["isLastCopy"])].flatten()
          assert len(zprime) == len(events)
          q0 = zprime.children[:, 0]
          q1 = zprime.children[:, 1]
          fatjets["genMatchFull"] = (fatjets.delta_r2(q0) < 0.8*0.8) & (fatjets.delta_r2(q1) < 0.8*0.8)

        print('dataset', dataset)
        candidatejet = fatjets[
            # https://github.com/DAZSLE/BaconAnalyzer/blob/master/Analyzer/src/VJetLoader.cc#L269
            (fatjets.pt > 250)
            & (abs(fatjets.eta) < 2.5)
            & (fatjets.rho >= -5.5)
            & (fatjets.rho <= -2)
            & (fatjets.genMatchFull if 'QCD' not in dataset else fatjets.pt > 0)
            # & fatjets.isLoose  # not always available

        ][:, 0:1]

        # basic jet selection
        selection.add('minjetkin', (
            (candidatejet.pt >= 525)
            & (candidatejet.msdcorr >= 40.)
            & (abs(candidatejet.eta) < 2.5)
        ).any())


        selection.add('jetid', candidatejet.isTight.any())
        #selection.add('gruddt', (candidatejet.n2ddt < 0.).any())


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
            #print(weight[0:100])
            output['templates'].fill(
                dataset=dataset,
                region=region,
                pt=normalize(candidatejet.pt, cut),
                msd=normalize(candidatejet.msdcorr, cut),
                gruddt=normalize(candidatejet.gruddt, cut),
                #n2ddt=normalize(candidatejet.n2ddt, cut),
                gru=normalize(candidatejet.twoProngGru, cut),
                rho=normalize(candidatejet.rho, cut),
            )

        for region in regions:
            fill(region)

        return output

    def postprocess(self, accumulator):
        return accumulator

