#!/usr/bin/env python
import json
import gzip
import uproot
#import numexpr
import numpy as np
from coffea import hist, lookup_tools
from coffea.util import load, save
from coffea.hist import plot
import os

corrections = load('data/corrections.coffea')

# add gruddt correction derived with 2017 QCD
shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_QCD_debug_5.coffea'))
values = list(shift_hist.values(overflow='none').values())[0]
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
corrections['2017_gruddt_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (rho_bins,pt_bins))
print(corrections['2017_gruddt_rho_pt'])
save(corrections, 'data/corrections_2.coffea')



