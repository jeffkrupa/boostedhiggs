#!/usr/bin/env python
import json
import gzip
import uproot
#import numexpr
import numpy as np
from numpy import inf
from coffea import hist, lookup_tools
from coffea.util import load, save
from coffea.hist import plot
import os

corrections = load('data/corrections.coffea')

# add gruddt correction derived with 2017 QCD
shift_hist = load(os.path.join(os.path.dirname(__file__), 'data', 'ddtmap_QCD_debug_6_v1.coffea'))
values = shift_hist.values(overflow='none')[()]
print(values)
#values[values==-inf] = 0.
print(values)
rho_bins = shift_hist.axis("jet_rho").edges(overflow='none')
pt_bins = shift_hist.axis("jet_pt").edges(overflow='none')
print(rho_bins,pt_bins)
corrections['2017_gruddt_rho_pt'] = lookup_tools.dense_lookup.dense_lookup(values, (pt_bins,rho_bins))
print(corrections['2017_gruddt_rho_pt']._evaluate(997.5,-3.170888587071259))

print(corrections['2017_gruddt_rho_pt']._evaluate(881.5, -3.859996923020683))
print(corrections['2017_gruddt_rho_pt']._evaluate(1039.0, -2.655261788438808))
#for pt in range(700, 710, 10):
#  for rho in range(15,75,1):
#    print(pt, float(rho)*-0.1, corrections['2017_gruddt_rho_pt']._evaluate(float(rho)*-0.1,float(pt)))
#print(corrections['2017_gruddt_rho_pt'].values)
save(corrections, 'data/corrections_4.coffea')



