import os
import numpy as np
from coffea import processor, util, hist
import json
#from coffea.nanoaod.methods import Candidate
from coffea.nanoevents import BaseSchema, NanoAODSchema, NanoEventsFactory
NanoAODSchema.warn_missing_crossrefs = True

#from coffea.nanoaod import NanoEvents
from boostedhiggs.ddt_processor import DDTProcessor
import argparse

def run_processor(year,selsamples,starti,endi,outname):
    p = DDTProcessor(year=year)
    
    files = {}
    with open('../data/fileset2017_preUL_QCD.json', 'r') as f:
        newfiles = json.load(f)
        files.update(newfiles)
    
    
    selfiles = {k: files[k][starti:endi] for k in selsamples}
    
    args = { "schema" : NanoAODSchema}  #{'nano': True, 'workers': 1, 'savemetrics': True}


    out = processor.run_uproot_job(selfiles, 'Events', p, processor.iterative_executor, args)
    
    util.save(out, '%s.coffea'%outname)


if __name__ == "__main__":
    #ex. python run_htt.py --year 2018 --starti 0 --endi -1 --selsamples GluGluHToTauTau --outname htt_runtest
    parser = argparse.ArgumentParser()
    parser.add_argument('--year',       dest='year',       default='2017',       help="year",        type=str)
    parser.add_argument('--starti',     dest='starti',     default=0,            help="start index", type=int)
    parser.add_argument('--endi',       dest='endi',       default=-1,           help="end index",   type=int)
    parser.add_argument('--selsamples', dest='selsamples', default=[],           help='selsamples',  nargs='+')
    parser.add_argument('--outname',    dest='outname',    default='ddt_test',   help='outname')
    args = parser.parse_args()

    run_processor(args.year,args.selsamples,args.starti,args.endi,args.outname)

