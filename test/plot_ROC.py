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

_decayFlavor = { #"uu":0,
                 #"dd":1,
                 "ll":1,
                 "cc":2,
                 "bb":3,
               }
fill_opts = {
    'edgecolor': (0,0,0,0.3),
    'alpha': 0.8
    }
err_opts = {
    #'label':'Stat. Unc.',
    #'hatch':'///',
    'facecolor':'none',
    'edgecolor':(0,0,0,.5),
    'linewidth': 0
    }
TAGGER = 'gru'
#TAGGER = 'deepTagMDZqq'
def roccer(bkg, sig, ax):
    
    fig,ax = plt.subplots()
    #for p, h in sig.items():
    b = bkg.values()[()]
    s = sig.values()[()]
    b/=np.sum(b)
    s/=np.sum(s)
 
    print(b,s)              
    bp = []
    sp = []
    for i in range(len(b)):
        bp.append(np.sum(b[0:i])) if 'n2' not in TAGGER else bp.append(np.sum(s[0:i]))
        sp.append(np.sum(s[0:i])) if 'n2' not in TAGGER else sp.append(np.sum(b[0:i]))
    plt.plot(bp,sp,)

    plt.xlim([0.,1.])
    plt.ylim([0.,1.])
    plt.xlabel('Signal efficiency')
    plt.ylabel('Background fake rate')
    plt.legend(loc='upper right',fontsize=12,title='decay',title_fontsize=14)
    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    plt.savefig('roc_w_vsQCD_%s.png'%(TAGGER))
    plt.savefig('roc_w_vsQCD_%s.pdf'%(TAGGER))

    plt.clf()
def drawSolo(h,sel,vars_cut,savename):
    print(h)
    h = h.integrate('region','signal')
    h = h.sum(*[ax for ax in h.axes() if ax.name not in ['process',TAGGER,'genflavor']],overflow='allnan')#.values(overflow='allnan')[()]
    print(h)
    for var,val in vars_cut.items():
        #if var=='pt':
        print('integrating ',var,val[0],val[1])
        h = h.integrate(var,slice(val[0],val[1]))
        #else: raise ValueError('only variable is pt')
    print(h)
    zqq = {}
    #for k,v in _decayFlavor.items():
    #   zqq["%s"%k] = h.integrate('genflavor',slice(v-0.5,v+0.5),overflow='none').integrate('process','zqq')
    #   #tmph = h.sum('genflavor',overflow='allnan').values(overflow='allnan')[('zqq',)]
    #   #zqq["%s"%k] = tmph/np.sum(tmph)
    #   #print(zqq["%s"%k])
    #zqq["inc"] = h.integrate('genflavor',slice(0.5,3+0.5),overflow='none').integrate('process','zqq')
    qcd = h.integrate('process','qcd')
    zqq = h.integrate('process','wqq')
    in_axis=h.axis(TAGGER).centers(overflow='none')

    roc = roccer(qcd, zqq, in_axis)
    #roc = roccer(qcd, {"qq":h.sum('genflavor',overflow='allnan').integrate('process','zqq')}, in_axis)

    fig,ax = plt.subplots()
    ymaxh = 0.1
    hist.plot1d(qcd,ax=ax,clear=True,density=True, )#error_opts=err_opts)
    hist.plot1d(h.integrate('process','zqq'),overlay='genflavor',ax=ax,clear=False,density=True)
    #for key, val in zqq.items():

    #   #plot gruddt distribution  position on each rho curve
    #   plt.hist(in_axis,val,label=_decayFlavor[key])
    #   #plot 95th quantile position on each rho curve
    #   #quantiles = np.cumsum(val) 
    #   #plt.scatter(gruaxis[quantiles.searchsorted(0.95)], val[quantiles.searchsorted(0.95)], marker='.', c='black', s=100)
    #   #ymaxh+=0.1 
    #plt.plot(in_axis,qcd,label='qcd')
    plt.xlabel(TAGGER)
    plt.ylabel('distribution (a.u.)')
    plt.legend(loc='upper right')

    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #addtext = plt.text(0.8, 1., "%i < $p_{T}$ < %i" % (vars_cut['pt'][0] , vars_cut['pt'][1]) ,fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes)
    fig.savefig("in_flavor_decay_%s.pdf"%(TAGGER,))

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

