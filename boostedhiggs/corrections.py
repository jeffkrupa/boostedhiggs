import os
import numpy as np
import awkward as ak
from coffea.util import load
from coffea import hist, lookup_tools

compiled = load(os.path.join(os.path.dirname(__file__), 'data', 'corrections_4.coffea'))

# hotfix some crazy large weights
compiled['2017_pileupweight']._values = np.minimum(5, compiled['2017_pileupweight']._values)
compiled['2018_pileupweight']._values = np.minimum(5, compiled['2018_pileupweight']._values)


def _msoftdrop_weight(pt, eta):
    gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
    cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
    fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])
    genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
    ptpow = np.power.outer(pt, np.arange(cpar.size))
    cenweight = np.dot(ptpow, cpar)
    forweight = np.dot(ptpow, fpar)
    weight = np.where(np.abs(eta) < 1.3, cenweight, forweight)
    return genw*weight


def corrected_msoftdrop(fatjets):
    sf_flat = _msoftdrop_weight(fatjets.pt.flatten(), fatjets.eta.flatten())
    sf_flat = np.maximum(1e-5, sf_flat)
    try:
        # old pancakes
        dazsle_msd = fatjets.msoftdrop_raw
    except AttributeError:
        dazsle_msd = (fatjets.subjets * (1 - fatjets.subjets.rawFactor)).sum().mass
    return dazsle_msd * ak.JaggedArray.fromoffsets(fatjets.array.offsets, sf_flat)

def shift(fatjets, algo, year='2017'):
    fatjets_msdcorr = corrected_msoftdrop(fatjets)
    fatjets_rhocorr = 2*np.log(fatjets_msdcorr/fatjets.pt)

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
    fatjets_msdcorr = corrected_msoftdrop(fatjets)
    fatjets_rhocorr = 2*np.log(fatjets_msdcorr/fatjets.pt)
    return compiled[f'{year}_n2b1_rho_pt'](fatjets.pt, fatjets_rhocorr)


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


def add_jetTriggerWeight(weights, jet_msd, jet_pt, year):
    jet_msd = jet_msd.pad(1, clip=True).fillna(0).flatten()
    jet_pt = jet_pt.pad(1, clip=True).fillna(0).flatten()
    nom = compiled[f'{year}_trigweight_msd_pt'](jet_msd, jet_pt)
    up = compiled[f'{year}_trigweight_msd_pt_trigweightUp'](jet_msd, jet_pt)
    down = compiled[f'{year}_trigweight_msd_pt_trigweightDown'](jet_msd, jet_pt)
    weights.add('jet_trigger', nom, up, down)
