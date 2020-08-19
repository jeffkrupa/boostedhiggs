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

plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        #'text.usetex': False,
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

def drawSolo(h,sel,var_name,var_label,plottitle,lumifb,vars_cut,regionsel,savename):
    exceptions = ['process', var_name]
    print(h)

    #get GRU distributions
    val_QCD = h.sum('region','jet_mass','jet_twoProngGru').values(overflow='all')[('qcd',)]
    val_QCD_bins = h.sum('region','jet_mass','jet_twoProngGru').axis('jet_n2b1').edges(overflow='all')
    print(val_QCD_bins)

    #make CDF
    qcd_maxval_temp = np.cumsum(val_QCD, axis=0)
    qcd_maxval = qcd_maxval_temp[-1]
    norma = qcd_maxval_temp / np.maximum(1e-10,qcd_maxval[np.newaxis])
    print(norma)
    s_over_sqrtb = [] 
    massrange=(70.,90.)
    sig='wqq'
    for pctl in range(0,59):

      ind_quant = norma.searchsorted(0.01*float(pctl))
      #print(ind_quant)
      quant = norma[ind_quant+1]

      GRUCUT = round(val_QCD_bins[ind_quant],4)
      #print (GRUCUT,quant,ind_quant)
      #print(slice(GRUCUT,1.))
      pass_GRU = h.sum('region','jet_twoProngGru')
      print('grucut',GRUCUT)
      count = pass_GRU.integrate('jet_n2b1',slice(0.,float(GRUCUT)))
      #print(count.axis('msd').edges(overflow='all'))
      print(count.integrate('jet_mass',slice(massrange[0],massrange[1])).values(overflow='all'))
      S = count.integrate('jet_mass',slice(massrange[0],massrange[1])).values(overflow='all')[(sig),]
      B = count.integrate('jet_mass',slice(massrange[0],massrange[1])).values(overflow='all')[('qcd'),]
      #print(S,B)
      s_over_sqrtb.append(S/np.sqrt(B+S))
      print(S , B, S/np.sqrt(B+S),pctl, GRUCUT)
      print("S = %.1f, B = %.1f, S/sqrt(B) = %.3f, percentile = %i, GRUcut = %.2f"%(S , B, S/np.sqrt(B+S),pctl, GRUCUT))

    plt.clf()
    fig,ax = plt.subplots()
    plt.plot(np.arange(0,59), np.array(s_over_sqrtb))
    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    plt.xlabel('N2 quantile')
    plt.ylabel('S/$\sqrt{S+B}$')
    plt.savefig('sensitivity_test.pdf')


def getPlots(args):
    print(args.lumi)
    lumifb = float(args.lumi)
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
    # normalize to lumi
    for h in hists_mapped.values():
        h.scale({p: lumifb for p in h.identifiers('process')}, axis="process")
    
    print(hists_mapped)
    # properties
    hist_name = args.hist
    var_name = args.var
    var_label = r"$%s$"%args.varlabel
    vars_cut =  {}
    #print(args.sel)
    if (len(args.sel)%3==0):
      for vi in range(int(len(args.sel)/3)):
        vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), float(args.sel[vi*3+2])]
    print(vars_cut)
    h = hists_mapped[hist_name]
    print(h)
        
    drawSolo(h,args.hist,var_name,var_label,args.title,lumifb,vars_cut,args.regions,savename)

    os.chdir(pwd)

if __name__ == "__main__":
    #ex. python plot_solo.py --hists htt_test --tag test --var jet_pt --varlabel 'p_{T}(jet)' --title Test --lumi 41.5 --sel lep_pt 20. 200. --regions hadel_signal  --hist trigeff --savetag leppt_20
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",      help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",           help="tag")
    parser.add_argument('--savetag',    dest='savetag',  default="",           help="savetag")
    parser.add_argument('--var',        dest='var',      default="",           help="var")
    parser.add_argument('--varlabel',   dest='varlabel', default="",           help="varlabel")
    parser.add_argument('--title',      dest='title',    default="",           help="title")
    parser.add_argument('--lumi',       dest='lumi',     default=50.,          help="lumi",       type=float)
    parser.add_argument('--sel',        dest='sel',      default='',           help='selection',  nargs='+')
    parser.add_argument('--regions',    dest='regions',  default='',           help='regionsel',  nargs='+')
    parser.add_argument('--hist',       dest='hist',     default='',           help='histname')
    args = parser.parse_args()

    getPlots(args)

