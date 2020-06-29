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


    tmp = h.integrate('pt',slice(525,575),overflow='none')
    qcd_rho = tmp.sum(*[ax for ax in h.axes() if ax.name not in ['rho','pt','gruddt']])#.values(overflow='allnan')[()]
    print(qcd_rho)
    qcd_gruddt = {}
    for i in range(len(_rhobins)):
       if i == len(_rhobins) - 1: continue
       qcd_gruddt["[%.1f,%.1f]"%(_rhobins[i],_rhobins[i+1])] = qcd_rho.integrate('rho',slice(_rhobins[i],_rhobins[i+1]),overflow='none').values(overflow='none')[()]
    print(qcd_gruddt)
    gruaxis=qcd_rho.axis('gruddt').centers(overflow='none')
    fig,ax = plt.subplots()
    for key, val in qcd_gruddt.items():
       val = val/np.sum(val)
       plt.plot(gruaxis,val,label=key.replace('_','< rho <'))

    plt.legend(loc='upper right')

    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    fig.savefig("solo_%s_%s_%s_lumi%i.pdf"%(sel,var_name,savename,lumifb))

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
    #h = hists_unmapped[hist_name]
    print('hi!',hists_unmapped[hist_name],hists_mapped[hist_name])
        
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

