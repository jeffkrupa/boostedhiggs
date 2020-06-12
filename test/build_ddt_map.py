import os
import numpy as np
from coffea import processor, util, hist
import json
from coffea.nanoaod.methods import Candidate
from coffea.util import load, save
from coffea.nanoaod import NanoEvents
from boostedhiggs.ddt_processor import DDTProcessor
import argparse
import scipy.ndimage as sc
import matplotlib.pyplot as plt

DISCRIMINATOR = "jet_twoProngGru"
HISTNAME = "jet_kin"

def plot(template, name):
    plt.clf()

    z = list(template._sumw.values())[0]
    
    #print(np.amin(z[z>0]), np.amax(z))   

    ax = hist.plot2d(template, xaxis = "jet_rho", patch_opts={"vmin":0.75, "vmax":0.99})
    cmstext = plt.text(0., 1., "CMS",fontsize=12,horizontalalignment='left',verticalalignment='bottom', fontweight='bold',transform=ax.transAxes)
    addtext = plt.text(0.11, 1., "Simulation Preliminary",fontsize=10,horizontalalignment='left',verticalalignment='bottom', style='italic', transform=ax.transAxes)

    plt.ylim(500,1200)
    plt.xlim(-5,-2)
    plt.savefig('plots/'+name+'.pdf')
    plt.savefig('plots/'+name+'.png')
 

def build_ddt_map(coffeafile, percentile, postfix):

    # adapted from https://github.com/SangeonPark/coffeandbacon/blob/master/analysis/DDT_Map_Derivation.ipynb
    hists = load(coffeafile)
    histo = hists[HISTNAME].sum('dataset','region',overflow='none') #hist of pt, rho, gru scores
    val_QCD = histo.values(overflow='allnan')[()]
    #val_QCD = val_QCD[()] 

    #sum hists and normalize
    qcd_maxval_temp = np.cumsum(val_QCD, axis=2)
    qcd_maxval = qcd_maxval_temp[:,:,-1]
    norma = qcd_maxval_temp / np.maximum(1e-3,qcd_maxval[:,:,np.newaxis])
    print(norma)
    print(norma.shape)
    hist_y_QCD = histo
    template = histo.sum(DISCRIMINATOR) #pt, rho base
    hist_y_QCD.clear()
    raw_cut = hist_y_QCD.sum(DISCRIMINATOR)
    hist_y_QCD._sumw = {():norma}


    res = np.apply_along_axis(lambda norma: norma.searchsorted(args.percentile), axis = 2, arr = norma)
    res[res>100]=0
    #print(res)
    print('shape',hist_y_QCD.values(overflow='none')[()].shape)
    print(hist_y_QCD.values()[()])
    def bineval(a):
        return hist_y_QCD.identifiers(DISCRIMINATOR)[a].lo

    binfunc = np.vectorize(bineval)

    qmap = binfunc(res)

    template.clear()
    template._sumw = {():qmap}
    template.label = 'GRU cut at ' + str(int(100*percentile)) + '%'

    plot(template, 'ddt_%i_%s'%(int(100*percentile),postfix))
    save(template, '../boostedhiggs/ddtmap_%s.coffea'%postfix) 

    '''
    #print(template)
    template.scale(-1.)
    print(raw_cut._sumw)
    print(template._sumw)
    flat_cut = raw_cut.add(template)
    template.scale(-1.)
    '''
    smooth_qmap = sc.filters.gaussian_filter(qmap,1)

    template.clear()
    template._sumw = {():smooth_qmap}
    template.label = 'GRU cut at ' + str(100*percentile) + '% (smoothed)'

    values_nonan = template.values()[()]

    plot(template, 'ddt_%i_smoothed_%s'%(int(100*percentile), postfix))
    save(template, '../boostedhiggs/ddtmap_smooth_%s.coffea'%postfix) 
    '''
    import ROOT
    outfile = ROOT.TFile("plots/n2ddtmap_2018bits_GaussianSmoothing1Sigma_CorrectVersion.root","recreate")
    outfile.cd()
    print(values_nonan.shape)
    h1 = ROOT.TH2F("h1","h1",52, -6, -2.1, 100, 300, 1300)
    for i in range(h1.GetNbinsX()):
        for j in range(h1.GetNbinsY()):
            h1.SetBinContent(i+1,j+1,values_nonan[j][i])
    h1.Write()
    outfile.Close()
    '''
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--coffeafile',       dest='coffeafile',        default='.',       help="coffea file location",        type=str)
    parser.add_argument('--percentile',       dest='percentile',       	default='.',       help="QCD efficiency",              type=float)
    parser.add_argument('--postfix',          dest='postfix',       	default='',        help="file name",                   type=str)
    args = parser.parse_args()
    build_ddt_map(args.coffeafile, args.percentile, args.postfix)
