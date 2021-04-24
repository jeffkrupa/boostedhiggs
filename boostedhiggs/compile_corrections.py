#!/usr/bin/env python
import json
import gzip
import uproot
#import numexpr
import numpy as np
from numpy import inf
from coffea import hist, lookup_tools
from coffea.lookup_tools import extractor
from coffea.util import load, save
from coffea.hist import plot
import os
import cloudpickle 
print(cloudpickle.__version__)
corrections = load('data/corrections_4.coffea')

# add gruddt correction derived with 2017 QCD
'''
shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_gru_QCD_debug_5.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_gruddt_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))

shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_n2_QCD_debug_5.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_n2b1_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))

shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_6Aug20_in_ddt.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_inddt_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))
from coffea.lookup_tools import extractor

shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_6Aug20_in_90pctl.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_inddt90pctl_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))
'''

shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_smooth_Apr21_2017_late.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_smooth_Apr21_2017_late_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))

shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_smooth_Apr21_2017_early.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_smooth_Apr21_2017_early_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))

shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_smooth_Sep20_2017.coffea'))
values = shift_hist.values(overflow='none')[()]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_smooth_Sep20_2017_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))

'''
from coffea.lookup_tools import extractor
ext = extractor()
ext.add_weight_sets(["2017_muTrigAbsEta_pt Mu50_PtEtaBins/abseta_pt_ratio data/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root"])
ext.finalize()
corrections['2017_muTrigAbsEta_pt'] = ext.make_evaluator()['2017_muTrigAbsEta_pt']
'''
save(corrections, 'data/corrections_4.coffea')


