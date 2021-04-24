from coffea import util,hist
import json
import os
import pprint
import numpy as np 

eosdir = "/eos/uscms/store/user/jkrupa/coffea_ak1/"
indir="14Apr21_2017_QCD_2"
os.system("mkdir -p %s" % indir)

chunk_size = 5

from os import listdir
from os.path import isfile, join
onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor")]

names = []
for name in onlyfiles:
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
  #flist[0]['templates'] = flist[0]['templates'].sum('genflavor',overflow='allnan')  
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
flist[0]['IN_Sep20_2017'].scale(scale1fb, 'dataset')
flist[0]['IN_Apr21_2017_late'].scale(scale1fb, 'dataset')
flist[0]['IN_Apr21_2017_early'].scale(scale1fb, 'dataset')
#flist[0]['event'].scale(scale1fb, 'dataset')
#flist[0]['in_v3'].scale(scale1fb, 'dataset')
#flist[0]['deepAK8'].scale(scale1fb, 'dataset')
#except:
#  flist[0]['jet_kin'].scale(scale1fb, 'dataset')

#for proc in flist[0]['templates'].axis('dataset').identifiers():

# print('process %s sum = %.5f' %(str(proc) , (flist[0]['templates'].integrate('region','ttbar_muoncontrol').sum('msd','gruddt','in_v3_ddt','mu_pt','mu_pfRelIso04_all','n2ddt').values()[(str(proc),)])))
#print('cutflow pre-scale')
#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(flist[0]['cutflow_ttbar_muoncontrol'])
#print(flist[0])
'''
for proc,flow in flist[0]['cutflow_muonCR'].items():
    if 'JetHT' in proc or 'SingleMuon' in proc: continue
    for cut, val in flow.items():
   
       flist[0]['cutflow_muonCR'][proc][cut] *= scale1fb[proc]

for proc,flow in flist[0]['cutflow_signal'].items():
    if 'JetHT' in proc or 'SingleMuon' in proc: continue
    for cut, val in flow.items():
   
       flist[0]['cutflow_signal'][proc][cut] *= scale1fb[proc]
for proc,flow in flist[0]['cutflow_VtaggingCR'].items():
    if 'JetHT' in proc or 'SingleMuon' in proc: continue
    for cut, val in flow.items():
   
       flist[0]['cutflow_VtaggingCR'][proc][cut] *= scale1fb[proc]
'''
#print('cutflow post-scale')
#pp = pprint.PrettyPrinter(indent=2)
#pp.pprint(flist[0]['cutflow_ttbar_muoncontrol'])
#print('VtaggingCR')
#pp.pprint(flist[0]['cutflow_VtaggingCR'])

#flist[0]['templates'].scale({k:0.63 for k,_ in flist[0]['sumw'].items() if 'QCD' in k}, 'dataset')
util.save(flist[0],'%s/hists_sum_gru2.coffea' % (indir))
for i,x in enumerate(chunk_names):
  os.system("rm %s/hists_sum_%i.coffea" % (indir,i*chunk_size))

