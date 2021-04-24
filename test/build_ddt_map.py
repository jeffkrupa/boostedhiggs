import os
import numpy as np
from coffea import processor, util, hist
import json
#from coffea.nanoaod.methods import Candidate
from coffea.util import load, save
#from coffea.nanoaod import NanoEvents
#from boostedhiggs.ddt_processor import DDTProcessor
import argparse
import scipy.ndimage as sc
import matplotlib.pyplot as plt
from collections import defaultdict, OrderedDict

DISCRIMINATOR = "jet_IN_Sep20_2017"
#DISCRIMINATOR = "jet_twoProngGru"
HISTNAME = "IN_Sep20_2017"

def plot(template, name):
    plt.clf()

    #ax = hist.plot2d(template, xaxis = "jet_rho",  patch_opts={"vmin":0.5, "vmax":0.99})#,xoverflow='all',yoverflow='all')
    ax = hist.plot2d(template, xaxis = "jet_rho",  patch_opts={"vmin":0.5, "vmax":1.})#,xoverflow='all',yoverflow='all')
    cmstext = plt.text(0.0, 1., "CMS",fontsize=12,horizontalalignment='left',verticalalignment='bottom', fontweight='bold',transform=ax.transAxes)
    addtext = plt.text(0.11, 1., "Simulation Preliminary",fontsize=10,horizontalalignment='left',verticalalignment='bottom', style='italic', transform=ax.transAxes)

    plt.ylim(200,1500)
    plt.xlim(-5.5,-2)
    plt.savefig('plots/'+name+'.pdf')
    plt.savefig('plots/'+name+'.png')
 

def build_ddt_map(coffeafile, percentile, postfix):

    # adapted from https://github.com/SangeonPark/coffeandbacon/blob/master/analysis/DDT_Map_Derivation.ipynb

    hists = load(coffeafile)
    print('imported coffea file:', hists.items())

    #load GRU distribution in each rho, pt bin
    histo = hists[HISTNAME].sum('dataset','region')
    val_QCD = histo.values(overflow='allnan')[()]

    #make CDF 
    qcd_maxval_temp = np.cumsum(val_QCD, axis=2)
    qcd_maxval = qcd_maxval_temp[:,:,-1]
    norma = qcd_maxval_temp / np.maximum(1e-10,qcd_maxval[:,:,np.newaxis])

    #make pt, rho base
    hist_y_QCD = hists[HISTNAME].sum('dataset','region')
    template = hist_y_QCD.sum(DISCRIMINATOR,) #pt, rho base
    hist_y_QCD.clear()
    hist_y_QCD._sumw = {():norma}

    #make an array of X percentile of GRU in (rho, pt) bin
    res = np.apply_along_axis(lambda norma: norma.searchsorted(args.percentile), axis = 2, arr = norma)
    res[res>100]=0

    #evaluation of GRU cut from quantile (trivial if GRU has 100 bins)
    def bineval(a):
        return hist_y_QCD.identifiers(DISCRIMINATOR,overflow='allnan')[a].lo

    binfunc = np.vectorize(bineval)
    
    qmap = binfunc(res)

    template.clear()
    template._sumw = {():qmap}
    template.label = 'IN cut at ' + str(int(100*percentile)) + '%'

    print(template.values()[()])
    save(template, '../boostedhiggs/ddtmap_%s.coffea'%postfix) 
    plot(template, 'ddt_%i_%s'%(int(100*percentile),postfix))
    '''
    #make an array
    print('rho')
    [print(round(rho,3),',', end="") for rho in template.axis('jet_rho').edges()]
    print('pt')
    [print(round(pt,3),',', end="")  for pt  in template.axis('jet_pt') .edges()]
    print('hist')
    print("{",end="")
    for ipt in range(0,len(template.axis('jet_pt').edges())-1):
        print("{", end="")
        for irho in range(0,len(template.axis('jet_rho').edges())-1):

           print(0. if abs(template.values()[()][ipt,irho]) == np.inf else round(template.values()[()][ipt,irho],3),',',end="")
        print("},",end=" ")

    print("}",end="")
    '''
    #smoothing
    #smooth_qmap = sc.filters.gaussian_filter(qmap,1)

    #template.clear()
    #template._sumw = {():smooth_qmap}
    #template.label = 'IN cut at ' + str(100*percentile) + '% (smoothed)'
    values_nonan = template.values()[()]
    print(values_nonan)
    save(template, '../boostedhiggs/ddtmap_smooth_%s.coffea'%postfix) 
    plot(template, 'ddt_%i_smoothed_%s'%(int(100*percentile), postfix))
    import ROOT
    outfile = ROOT.TFile("plots/in_v3_ddtmap.root","recreate")
    outfile.cd()
    print(values_nonan.shape)
    h1 = ROOT.TH2F("h1","h1",100, 200, 1500, 180, -5.5, -2.)
    print(h1.GetNbinsX()-1,h1.GetNbinsY()-1)
    for i in range(h1.GetNbinsX()-1):
        for j in range(h1.GetNbinsY()-1): 
             #print(values_nonan[i][j])
             h1.SetBinContent(i+1,j+1,values_nonan[i][j])

    h1.Write()
    outfile.Close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--coffeafile',       dest='coffeafile',        default='.',       help="coffea file location",        type=str)
    parser.add_argument('--percentile',       dest='percentile',       	default='.',       help="QCD efficiency",              type=float)
    parser.add_argument('--postfix',          dest='postfix',       	default='',        help="file name",                   type=str)
    args = parser.parse_args()


    build_ddt_map(args.coffeafile, args.percentile, args.postfix)
