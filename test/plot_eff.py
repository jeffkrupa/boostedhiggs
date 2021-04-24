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
        'font.size': 22,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        #'text.usetex': True,
        })

_ptBins = [500,575,625,700,800,1000]

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
TAGGER = 'in_v3'

mc_processes = ['zqq','wqq','qcd','st','wlnu','tt']#'tttoleptonic','tttosemileptonic','tttohadronic']
data_processes = ['JetHT',]

def integrate(h,processes,ml,mh):
    val = h.integrate('process',processes).integrate('msd',slice(ml,mh)).values()[()]
    unc = np.sqrt(val)
    return val, unc #h.integrate('process',processes).integrate('msd',slice(ml,mh)).values()[()]


def get_effs(hpass,hfail):
    
    fig,ax = plt.subplots()
 
    
    for ip in range(len(_ptBins)):
       if ip == len(_ptBins) - 1: continue


       htmppass = hpass.integrate('pt',slice(_ptBins[ip],_ptBins[ip+1]))
       htmpfail = hfail.integrate('pt',slice(_ptBins[ip],_ptBins[ip+1]))
       msd = htmppass.axis('msd').edges()
       double_ratio = np.array(msd)
       double_ratio_err = np.array(msd)
       for im in range(len(msd)-1):


           pass_data, pass_data_err = integrate(htmppass,data_processes,msd[im],msd[im+1])
           pass_mc,   pass_mc_err   = integrate(htmppass,mc_processes,msd[im],msd[im+1])
           fail_data, fail_data_err = integrate(htmpfail,data_processes,msd[im],msd[im+1]) 
           fail_mc,   fail_mc_err   = integrate(htmpfail,mc_processes,msd[im],msd[im+1])         
 
           double_ratio[im]     = (pass_data/pass_mc) / (fail_data/fail_mc)
           double_ratio_err[im] = np.sqrt( (pass_data_err/pass_data)**2 + (pass_mc_err/pass_mc)**2 + (fail_data_err/fail_data)**2 + (fail_mc_err/fail_mc)**2) 
           #double_ratio[im], double_ratio_err[im] = integrate(htmppass,data_processes,msd[im],msd[im+1]) / integrate(htmppass,mc_processes,msd[im],msd[im+1]) * integrate(htmpfail,mc_processes,msd[im],msd[im+1]) / integrate(htmpfail,data_processes,msd[im],msd[im+1])  
       


       plt.errorbar(msd,double_ratio,yerr=double_ratio_err,label=r'%i < pT < %i'%(_ptBins[ip],_ptBins[ip+1]))

    plt.xlim(40,290)
    plt.ylim(0,2)
    plt.ylabel('Data/MC pass/fail')
    plt.xlabel('Jet mass (GeV)')
    plt.legend(loc='upper right',fontsize=16,title='pT bin',title_fontsize=18)
    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    plt.savefig('effs_by_pTbin.png')
    plt.savefig('effs_by_pTbin.pdf')

    plt.clf()
def drawSolo(h,sel,vars_cut,savename):
    print(h)
    h = h.integrate('region','signal')
    h = h.sum(*[ax for ax in h.axes() if ax.name not in ['process','msd','pt','n2ddt']])#.values(overflow='allnan')[()]
    hpass = h.integrate('n2ddt',slice(-2,0))
    hfail = h.integrate('n2ddt',slice(0,2))

    get_effs(hpass,hfail)

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

