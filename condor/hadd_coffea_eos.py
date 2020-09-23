from coffea import util,hist
import json
import os
import pprint
import numpy as np 

eosdir = "/eos/uscms/store/user/jkrupa/coffea/"
#indir = "6Aug20_in_ddt_all_isoTight_dR08matching_looserTau_v3"
indir = "6Aug20_in_ddt_all_isoTight_dR08matching_looserTau_v3_w10ptcl_withgrun2_v2"
indir = "6Aug20_in_ddt_all_isoTight_looserTau_v3_w10ptclFIX_withgrun2_matchedBoson06" 
indir = "6Aug20_in_ddt_all_isoTight_looserTau_v3_w10ptclFIX_withgrun2_matchedBoson06_v3"
indir = "22Sep20_QCD"
#indir = "6Aug20_in_ddt_all_isoTight_metcut_dR08matching_v3" #"6Aug20_debugQCDmuonCR_17_debug_Jeff_any" #6Aug20_debugQCDmuonCR_14_QCD_"
os.system("mkdir -p %s" % indir)

chunk_size = 5

from os import listdir
from os.path import isfile, join
onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor")]

names = []
for name in onlyfiles:
  #if 'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8_4' in name: continue
  #if 'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_0' in name: continue
  #if 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_17' in name: continue #if 'JetHT' in name or 'SingleMuon' in name: continue
  #if 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_1' in name or 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_2' in name or 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_4' in name or 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_5' in name: continue
  #if 'QCD' in name:
  #if (os.path.isfile("%s%s/%s.coffea" % (eosdir,indir,name))): names.append("%s%s/%s.coffea" % (eosdir,indir,name))
  #if 'TT' in name: continue
  #if 'LNu' in name: continue
  #if 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_14' in name: continue
  #if 'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8_27' in name: continue
  names.append("%s%s/%s.coffea" % (eosdir,indir,name))
from collections import Counter 
os.system('ls %s%s/ | wc -l'%(eosdir,indir))

chunk_names = []
for i in range(0,len(names),chunk_size):
  print('Chunk',i)
  #if i>0: continue
  if (i+chunk_size<len(names)):
    flist = [ util.load(x) for x in names[i:i+chunk_size] ]
  else:
    flist = [ util.load(x) for x in names[i:] ]

  for key in flist[0]:
    if isinstance(key, hist.Hist):
      for fi in range(1,len(flist)):
        flist[0][key].add(flist[fi][key])
    else:
      for fi in range(1,len(flist)):
        flist[0][key] = flist[0][key] + flist[fi][key]
  
  #print(flist[0])
  flist[0]['templates'] = flist[0]['templates'].sum('pt',overflow='all')  
  #flist[0]['templates'].rebin('pt',2)
  #flist[0]['templates'].rebin('mu_pt',2
  #flist[0]['in_v3'] = flist[0]['in_v3'].sum('n2',overflow='allnan')
  util.save(flist[0],'%s/hists_sum_%i.coffea' % (indir,i))

  for f in flist:
    del f

  chunk_names.append('%s/hists_sum_%i.coffea' % (indir,i))

print(chunk_names)

flist = [ util.load(x) for x in chunk_names ]

for key in flist[0]:
  if isinstance(key, hist.Hist):
    for fi in range(1,len(flist)):
      flist[0][key].add(flist[fi][key])
  else:
    for fi in range(1,len(flist)):
      flist[0][key] = flist[0][key] + flist[fi][key]
  
#print('hists',flist[0])
xs = {}
with open('../data/xsec.json', 'r') as f:
   xs = json.load(f)

#flist[0]['templates']   
scale1fb = {k: xs[k] * 1000 / w for k, w in flist[0]['sumw'].items()}
#print(flist[0]['templates'].integrate('region','ttbar_muoncontrol').sum('msd','gruddt','in_v3_ddt','mu_pt','mu_pfRelIso04_all','n2ddt').values())
#for proc in flist[0]['templates'].axis('dataset').identifiers():

# print('process %s sum = %.5f' %(str(proc) , (flist[0]['templates'].integrate('region','ttbar_muoncontrol').sum('msd','gruddt','in_v3_ddt','mu_pt','mu_pfRelIso04_all','n2ddt').values()[(str(proc),)])))
#try:
flist[0]['templates'].scale(scale1fb, 'dataset')
flist[0]['muon'].scale(scale1fb, 'dataset')
flist[0]['event'].scale(scale1fb, 'dataset')
flist[0]['in_v3'].scale(scale1fb, 'dataset')
flist[0]['deepAK8'].scale(scale1fb, 'dataset')
#except:
#  flist[0]['jet_kin'].scale(scale1fb, 'dataset')

#for proc in flist[0]['templates'].axis('dataset').identifiers():

# print('process %s sum = %.5f' %(str(proc) , (flist[0]['templates'].integrate('region','ttbar_muoncontrol').sum('msd','gruddt','in_v3_ddt','mu_pt','mu_pfRelIso04_all','n2ddt').values()[(str(proc),)])))
#print('cutflow pre-scale')
#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(flist[0]['cutflow_ttbar_muoncontrol'])
for proc,flow in flist[0]['cutflow_ttbar_muoncontrol'].items():
    if 'JetHT' in proc or 'SingleMuon' in proc: continue
    for cut, val in flow.items():
   
       flist[0]['cutflow_ttbar_muoncontrol'][proc][cut] *= scale1fb[proc]

for proc,flow in flist[0]['cutflow_signal'].items():
    if 'JetHT' in proc or 'SingleMuon' in proc: continue
    for cut, val in flow.items():
   
       flist[0]['cutflow_signal'][proc][cut] *= scale1fb[proc]
print('cutflow post-scale')
pp = pprint.PrettyPrinter(indent=2)
#pp.pprint(flist[0]['cutflow_ttbar_muoncontrol'])
print('signal')
pp.pprint(flist[0]['cutflow_signal'])

'''for key,val in scale1fb.items():
  for keyp,valp in flist[0]['cutflow_ttbar_muoncontrol'].items():
     if key==keyp: 
        print(flist[0]['cutflow_ttbar_muoncontrol'][keyp])
        print (key,val); 
        for x,y in flist[0]['cutflow_ttbar_muoncontrol'][keyp].items():
            flist[0]['cutflow_ttbar_muoncontrol'][keyp][y] =  val*flist[0]['cutflow_ttbar_muoncontrol'][keyp][y]
  for keyp,valp in flist[0]['cutflow_ttbar_muoncontrol'].items():
     if key==keyp: 
        print (key,val); flist[0]['cutflow_ttbar_muoncontrol'][keyp] = val*flist[0]['cutflow_ttbar_muoncontrol'][keyp]
'''
#flist[0]['templates'].scale({k:0.63 for k,_ in flist[0]['sumw'].items() if 'QCD' in k}, 'dataset')
util.save(flist[0],'%s/hists_sum_gru2.coffea' % (indir))
for i,x in enumerate(chunk_names):
  os.system("rm %s/hists_sum_%i.coffea" % (indir,i*chunk_size))

