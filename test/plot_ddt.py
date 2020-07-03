from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.styles.ROOT)
import numpy as np

from coffea import hist
from coffea.util import load

import pickle
import gzip
import math

import argparse
import processmap
from hists_map import *

_rhobins = [-5.5,-5.,-4.5,-4.,-3.5,-3.,-2.5,-2.0,]
plt.rcParams.update({
        'font.size': 22,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        #'text.usetex': True,
        })

fill_opts = {
    'edgecolor': (0,0,0,0.3),
    'alpha': 0.8
    }
err_opts = {
    'label':'Stat. Unc.',
    'hatch':'///',
    'facecolor':'none',
    'edgecolor':(0,0,0,.5),
    'linewidth': 0
    }
def drawSolo(h,sel,vars_cut,savename):

    h = h.sum(*[ax for ax in h.axes() if ax.name not in ['rho','pt','gruddt']]).rebin('gruddt',2) #.values(overflow='allnan')[()]
    for var,val in vars_cut.items():
        #if var=='pt':
        print('integrating ',var,val[0],val[1])
        h = h.integrate(var,slice(val[0],val[1]))
        #else: raise ValueError('only variable is pt')
    print(h)
    qcd_gruddt = {}
    for i in range(len(_rhobins)):
       if i == len(_rhobins) - 1: continue
       qcd_gruddt["%.1f_%.1f"%(_rhobins[i],_rhobins[i+1])] = h.integrate('rho',slice(_rhobins[i],_rhobins[i+1]),overflow='none').values(overflow='none')[()]

    gruaxis=h.axis('gruddt').centers(overflow='none')

    fig,ax = plt.subplots()
    ymaxh = 0.1
    for key, val in qcd_gruddt.items():
       val = val/np.sum(val)

       #plot gruddt distribution  position on each rho curve
       plt.plot(gruaxis,val,label=r"%s < $\rho$ < %s" % (key.split('_')[0], key.split('_')[1] ))
       #plot 95th quantile position on each rho curve
       quantiles = np.cumsum(val) 
       plt.scatter(gruaxis[quantiles.searchsorted(0.95)], val[quantiles.searchsorted(0.95)], marker='.', c='black', s=100)
       ymaxh+=0.1 
    plt.xlabel('$GRU^{DDT}$')
    plt.ylabel('distribution (a.u.)')
    plt.legend(loc='upper right')

    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    addtext = plt.text(0.8, 1., "%i < $p_{T}$ < %i" % (vars_cut['pt'][0] , vars_cut['pt'][1]) ,fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes)
    fig.savefig("gruddt_%s.pdf"%(savename,))

def getPlots(args):
    tag = args.tag
    savename = args.savetag

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    # open hists
    hists_unmapped = load('%s.coffea'%args.hists)
    os.chdir(odir)

    # map to hists
    hists_mapped = {}
    for key, val in hists_unmapped.items():
        if isinstance(val, hist.Hist):
            hists_mapped[key] = processmap.apply(val)
    # properties
    hist_name = args.hist
    vars_cut =  {}
    if (len(args.sel)%3==0):
      for vi in range(int(len(args.sel)/3)):
        vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), float(args.sel[vi*3+2])]
    h = hists_mapped[hist_name]
        
    drawSolo(h,args.hist,vars_cut,savename)

    os.chdir(pwd)

if __name__ == "__main__":
    #python plot_ddt.py --hists ../condor/QCD_debug_6_v2/hists_sum_gru2 --tag QCD_debug_6_v2 --hist templates --sel pt 525 575 --savetag 'gruddt_distribution_525_575'
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",      help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",           help="tag")
    parser.add_argument('--savetag',    dest='savetag',  default="",           help="savetag")
    parser.add_argument('--sel',        dest='sel',      default='',           help='selection',  nargs='+')
    parser.add_argument('--hist',       dest='hist',     default='',           help='histname')
    args = parser.parse_args()

    getPlots(args)

