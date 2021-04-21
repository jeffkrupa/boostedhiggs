import os
import numpy as np
import awkward as ak
import gzip
import pickle
from coffea.lookup_tools.lookup_base import lookup_base
from coffea.util import load, save
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

compiled = load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'corrections_4.coffea'))
#with gzip.open(os.path.join(os.path.dirname(__file__), 'data', 'corrections.pkl.gz')) as fin:
#    compiled = pickle.load(fin)
#print(compiled)
# hotfix some crazy large weights
compiled['2017_pileupweight']._values = np.minimum(5, compiled['2017_pileupweight']._values)
compiled['2018_pileupweight']._values = np.minimum(5, compiled['2018_pileupweight']._values)


class SoftDropWeight(lookup_base):
    def _evaluate(self, pt, eta):
        gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
        cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
        fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])
        genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
        ptpow = np.power.outer(pt, np.arange(cpar.size))
        cenweight = np.dot(ptpow, cpar)
        forweight = np.dot(ptpow, fpar)
        weight = np.where(np.abs(eta) < 1.3, cenweight, forweight)
        return genw*weight


_softdrop_weight = SoftDropWeight()


def corrected_msoftdrop(fatjets):
    sf = _softdrop_weight(fatjets.pt, fatjets.eta)
    sf = np.maximum(1e-5, sf)
    try:
        # pancakes have the raw value
        dazsle_msd = fatjets.msoftdrop_raw
    except AttributeError:
        # for nanoaod we have to work back to it
        # TODO: this should be ak.sum(..., axis=-1) but not working
        dazsle_msd = (fatjets.subjets * (1 - fatjets.subjets.rawFactor)).sum().mass
    return dazsle_msd * sf

def shift(fatjets, algo, year='2017'):
    fatjets_msdcorr = corrected_msoftdrop(fatjets)
    fatjets_rhocorr = 2*np.log(fatjets_msdcorr/fatjets.pt)
    #2017_inddt90pctl_rho_pt
    return compiled['%s_%s_rho_pt'%(year,algo)](fatjets.pt, fatjets_rhocorr)

def inddt_shift(fatjets, year='2017'):
    fatjets_msdcorr = corrected_msoftdrop(fatjets)
    fatjets_rhocorr = 2*np.log(fatjets_msdcorr/fatjets.pt)

    return compiled[f'2017_gruddt_rho_pt'](fatjets.pt, fatjets_rhocorr)

def gruddt_shift(fatjets, year='2017'):
    fatjets_msdcorr = corrected_msoftdrop(fatjets)
    fatjets_rhocorr = 2*np.log(fatjets_msdcorr/fatjets.pt)

    return compiled[f'2017_gruddt_rho_pt'](fatjets.pt, fatjets_rhocorr)

def n2ddt_shift(fatjets, year='2017'):
    return compiled[f'{year}_n2ddt_rho_pt'](fatjets.qcdrho, fatjets.pt)


def add_pileup_weight(weights, nPU, year='2017', dataset=None):
    if year == '2017' and dataset in compiled['2017_pileupweight_dataset']:
        weights.add(
            'pileup_weight',
            compiled['2017_pileupweight_dataset'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puUp'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puDown'][dataset](nPU),
        )
    else:
        weights.add(
            'pileup_weight',
            compiled[f'{year}_pileupweight'](nPU),
            compiled[f'{year}_pileupweight_puUp'](nPU),
            compiled[f'{year}_pileupweight_puDown'](nPU),
        )


def add_VJets_NLOkFactor(weights, genBosonPt, year, dataset):
    if year == '2017' and 'ZJetsToQQ_HT' in dataset:
        nlo_over_lo_qcd = compiled['2017_Z_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2017' and 'WJetsToQQ_HT' in dataset:
        nlo_over_lo_qcd = compiled['2017_W_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2016' and 'DYJetsToQQ' in dataset:
        nlo_over_lo_qcd = compiled['2016_Z_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2016' and 'WJetsToQQ' in dataset:
        nlo_over_lo_qcd = compiled['2016_W_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    else:
        return
    weights.add('VJets_NLOkFactor', nlo_over_lo_qcd * nlo_over_lo_ewk)

def add_singleMuTriggerWeight(weights, mu_absEta, mu_pt, year):
    mu_absEta = abs(mu_absEta.pad(1, clip=True).fillna(0).flatten()) 
    mu_pt     = mu_pt.pad(1, clip=True).fillna(0).flatten() 
    nom       = compiled[f'{year}_muTrigAbsEta_pt'](mu_absEta, mu_pt)
    weights.add('singleMuonTrigger',nom)

def add_jetTriggerWeight(weights, jet_msd, jet_pt, year):
    nom = compiled[f'{year}_trigweight_msd_pt'](jet_msd, jet_pt)
    up = compiled[f'{year}_trigweight_msd_pt_trigweightUp'](jet_msd, jet_pt)
    down = compiled[f'{year}_trigweight_msd_pt_trigweightDown'](jet_msd, jet_pt)
    weights.add('jet_trigger', nom, up, down)



from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory

jec_name_map = {
    'JetPt': 'pt',
    'JetMass': 'mass',
    'JetEta': 'eta',
    'JetA': 'area',
    'ptGenJet': 'pt_gen',
    'ptRaw': 'pt_raw',
    'massRaw': 'mass_raw',
    'Rho': 'event_rho',
    'METpt': 'pt',
    'METphi': 'phi',
    'JetPhi': 'phi',
    'UnClusteredEnergyDeltaX': 'MetUnclustEnUpDeltaX',
    'UnClusteredEnergyDeltaY': 'MetUnclustEnUpDeltaY',
}

def jet_factory_factory(files):
    ext = extractor()
    ext.add_weight_sets([f"* * {file}" for file in files])
    ext.finalize()

    jec_stack = JECStack(ext.make_evaluator())
    return CorrectedJetsFactory(jec_name_map, jec_stack)


def add_jec_variables(jets, jet_matched, event_rho):
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jet_matched.pt, 0), np.float32)
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    #jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets

jet_factory = { 

    "2017mc": jet_factory_factory(
        files=[
            # https://github.com/cms-jet/JECDatabase/raw/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt.gz"),
            # https://github.com/cms-jet/JECDatabase/raw/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt.gz"),
            # https://raw.githubusercontent.com/cms-jet/JECDatabase/master/textFiles/Fall17_17Nov2017_V32_MC/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt
            os.path.join(DATA_DIR, "RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.junc.txt.gz"),
            # https://github.com/cms-jet/JECDatabase/raw/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.junc.txt.gz"),
            # https://github.com/cms-jet/JRDatabase/raw/master/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_PtResolution_AK4PFchs.txt
            os.path.join(DATA_DIR, "Fall17_V3b_MC_PtResolution_AK4PFchs.jr.txt.gz"),
            # https://github.com/cms-jet/JRDatabase/raw/master/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_SF_AK4PFchs.txt
            os.path.join(DATA_DIR, "Fall17_V3b_MC_SF_AK4PFchs.jersf.txt.gz"),
        ]
    ),
    "2017mcNOJER": jet_factory_factory(
        files=[
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt.gz"),
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt.gz"),
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.junc.txt.gz"),
        ]
    ),
}

fatjet_factory = {

    "2017mc": jet_factory_factory(
        files=[
            # https://github.com/cms-jet/JECDatabase/raw/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFPuppi.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFPuppi.jec.txt.gz"),
            # https://github.com/cms-jet/JECDatabase/raw/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.jec.txt.gz"),
            # https://raw.githubusercontent.com/cms-jet/JECDatabase/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_UncertaintySources_AK8PFPuppi.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_UncertaintySources_AK8PFPuppi.junc.txt.gz"),
            # https://github.com/cms-jet/JECDatabase/raw/master/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.junc.txt.gz"),
            # https://github.com/cms-jet/JRDatabase/raw/master/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_PtResolution_AK8PFPuppi.txt
            os.path.join(DATA_DIR, "Fall17_V3b_MC_PtResolution_AK8PFPuppi.jr.txt.gz"),
            # https://github.com/cms-jet/JRDatabase/raw/master/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_SF_AK8PFPuppi.txt
            os.path.join(DATA_DIR, "Fall17_V3b_MC_SF_AK8PFPuppi.jersf.txt.gz"),
        ]
    ),
    "2017mcNOJER": jet_factory_factory(
        files=[
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFPuppi.jec.txt.gz"),
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.jec.txt.gz"),
            os.path.join(DATA_DIR, "Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.junc.txt.gz"),
        ]
    ),
}
met_factory = CorrectedMETFactory(jec_name_map)
