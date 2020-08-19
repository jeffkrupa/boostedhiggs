from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.ticker as mpltick
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.styles.ROOT)
import numpy as np

from coffea import hist, processor
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
line_opts = {
    'color': 'aquamarine',
    'linewidth':2,
    #'marker':'None',
    'linestyle':'dashed',
    'drawstyle':'steps'}

line_opts_nocolor = {
    'linewidth':2,
    #'marker':'None',
    'linestyle':'dashed',
    'drawstyle':'steps'}

special_char_list = ['_','-']

def getindices(s):
    return [i for i, c in enumerate(s) if c.isupper() or c in special_char_list]

def drawCutflow(h,plottitle,lumifb,regionsel):

    cutlist = []
    for c in h[list(h.keys())[0]]:
        cutlist.append(c)

    cutnames = []
    newlinei = 8
    for c in cutlist:
        if (len(c)<newlinei):
            cutnames.append(c)
        else:
            isCut = False
            cuti = getindices(c)
            if not cuti:
                cutnames.append(c)
                isCut = True
            else:
                for ic in cuti:
                    if (ic>=newlinei):
                        if (c[ic] in special_char_list): cutnames.append(c[:ic]+"\n"+c[ic+1:])
                        else: cutnames.append(c[:ic]+"\n"+c[ic:])
                        isCut = True
                        break
            if not isCut:
                cutnames.append(c)
    print(cutlist)
    print(cutnames)

    cuthist_raw = hist.Hist(
        'Events',
        hist.Cat('dataset', 'Dataset'),
        hist.Bin('cut','Cut',len(cutlist),-0.5,float(len(cutlist))-0.5),
    )

    for s in h:
        for ic,c in enumerate(cutlist):
            cuthist_raw.fill(
                dataset=s,
                cut=ic,
                weight=h[s][c],
            )

    x = processmap.apply(cuthist_raw)
    for ih,hkey in enumerate(x.identifiers('process')):
        x.identifiers('process')[ih].label = process_latex[hkey.name]

    mc_processes = ['zqq','wqq','qcd','tt','st','wlnu',]
    data_processes = ['SingleMuon','JetHT',]

    mc = x.remove(data_processes,'process')
    data = x.remove(mc_processes + ['SingleMuon'] if 'signal' in regionsel else mc_processes + ['JetHT'],'process')
    mc.scale({p: lumifb for p in mc.identifiers('process')}, axis="process")
    mc.axis('process').sorting = 'integral'
    fig,ax = plt.subplots()
    hist.plot1d(mc,
                overlay='process',ax=ax,
                clear=False,
                stack=True,
                fill_opts=fill_opts,
                error_opts=err_opts,
                )
    hist.plot1d(data, ax=ax, clear=False, line_opts=line_opts)
    all_bkg = 0
    #for key,val in x[nosig].values().items():
    #    all_bkg+=val.sum()
    #x_nobkg = x[nobkg]
    #x_nobkg.scale(50)

    #all_sig = 0
    #for key,val in x_nobkg.values().items():
    #    all_sig +=val.sum()
    #print('%.4f %.4f %.4f'%(all_bkg,all_sig,all_sig/math.sqrt(all_bkg)))
    #hist.plot1d(x_nobkg,ax=ax,overlay='process',clear=False,line_opts=line_opts)

    ax.autoscale(axis='x', tight=True)
    #ax.set_xlim(20, 200)
    ax.set_ylim(0, None)
    plt.xticks(range(len(cutnames)), cutnames, rotation=60)
    #ax.ticklabel_format(axis='x', style='sci')
    old_handles, old_labels = ax.get_legend_handles_labels()
    new_labels = []
    for xl in old_labels:
        if ('ggH' in xl): xl = xl + " (x 50)"
        new_labels.append(xl)
    leg = ax.legend(handles=old_handles,labels=new_labels,title=r'$%s$'%plottitle,frameon=True,framealpha=1.0,facecolor='white',loc='lower left')
    lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    cmstext = plt.text(0., 1., "CMS",fontsize=19,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    #ax.semilogy()
    ax.set_yscale("log")
    minvals = []
    for xd in x.values():
        minvals.append(min(np.trim_zeros(x.values()[xd])))
    decsplit = str(min(minvals)).split('.')
    if (int(decsplit[0])==0):
        logmin = 0.1**float(len(decsplit[1])-len(decsplit[1].lstrip('0'))+2)
    else:
        logmin = 10.**float(len(decsplit[0])-1)
    ax.set_ylim(logmin/10., None)

    plt.minorticks_on()
    locmaj = mpltick.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
    ax.yaxis.set_major_locator(locmaj)

    locmin = mpltick.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                      numticks=100)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
    fig.savefig("cutflow_%s_lumi%i_logy.pdf"%(regionsel,lumifb))

#----------------------------
    fig,ax = plt.subplots()

    x.scale({(str(xd).split('\''))[1]:1./x.values()[xd][0] for xd in x.values()},'process')

    hist.plot1d(x,
                overlay='process',ax=ax,
                clear=False,
                stack=False,
                #fill_opts=fill_opts,
                #error_opts=err_opts,
                line_opts=line_opts_nocolor
                )

    ax.autoscale(axis='x', tight=True)
    #ax.set_xlim(20, 200)
    ax.set_ylim(0, None)
    plt.xticks(range(len(cutnames)), cutnames, rotation=60)
    #ax.ticklabel_format(axis='x', style='sci')
    old_handles, old_labels = ax.get_legend_handles_labels()
    new_labels = []
    for xl in old_labels:
        new_labels.append(xl)
    leg = ax.legend(handles=old_handles,labels=new_labels,title=r'$%s$'%plottitle,frameon=True,framealpha=1.0,facecolor='white',loc='lower left')
    lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    ax.semilogy()
    minvals = []
    for xd in x.values():
        minvals.append(min(np.trim_zeros(x.values()[xd])))
    ax.set_ylim(10.**float(math.floor(math.log10(min(minvals)))))
    ylo, yhi = ax.get_ylim()
    nticksy = int(math.floor(math.log10(yhi))-math.floor(math.log10(ylo)))+2

    plt.minorticks_on()
    locmaj = mpltick.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
    ax.yaxis.set_major_locator(locmaj)

    locmin = mpltick.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                      numticks=100)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
    fig.savefig("cuteff_%s_lumi%i_logy.pdf"%(regionsel,lumifb))



def getPlots(args):
    print(args.lumi)
    lumifb = float(args.lumi)
    tag = args.tag

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    # open hists
    hists_unmapped = load('%s.coffea'%args.hists)
    os.chdir(odir)

    # map to hists
    dict_mapped = {}
    for key, val in hists_unmapped.items():
        if key in args.regions:
            if isinstance(val, processor.accumulator.defaultdict_accumulator):
                #for proc in processmap.process_map:
                dict_mapped[key] = val
    print(dict_mapped.values())
    # normalize to lumi
    #for h in hists_mapped.values():
    #    h.scale({p: lumifb for p in h.identifiers('process')}, axis="process")
    
    for r in args.regions:
        drawCutflow(dict_mapped[r],args.title,lumifb,r)

    os.chdir(pwd)

if __name__ == "__main__":
    #ex. python plot_cutflow.py --hists htt_test --tag test --title Test --lumi 41.5 --regions cutflow_hadel
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",      help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",           help="tag")
    parser.add_argument('--title',      dest='title',    default="",           help="title")
    parser.add_argument('--lumi',       dest='lumi',     default=50.,          help="lumi",       type=float)
    parser.add_argument('--regions',    dest='regions',  default='',           help='regionsel',  nargs='+')
    args = parser.parse_args()

    getPlots(args)

#        output['cutflow_hadel'][dataset]['none'] += float(weights.weight().sum())
#        output['cutflow_hadmu'][dataset]['none'] += float(weights.weight().sum())
#        output['cutflow_hadel_cr_b'][dataset]['none'] += float(weights.weight().sum())
#        output['cutflow_hadmu_cr_b'][dataset]['none'] += float(weights.weight().sum())
#        output['cutflow_hadhad'][dataset]['none'] += float(weights.weight().sum())

