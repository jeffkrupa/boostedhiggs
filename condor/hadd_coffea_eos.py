from coffea import util,hist
import json
import os
import pprint


eosdir = "/eos/uscms/store/user/jkrupa/coffea/"
indir = "6Aug20_debugQCDmuonCR_4"
os.system("mkdir -p %s" % indir)

chunk_size = 10

from os import listdir
from os.path import isfile, join
onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor")]

names = []
for name in onlyfiles:
  #if (os.path.isfile("%s%s/%s.coffea" % (eosdir,indir,name))): names.append("%s%s/%s.coffea" % (eosdir,indir,name))
  #if 'TT' in name: continue
  #if 'LNu' in name: continue
  names.append("%s%s/%s.coffea" % (eosdir,indir,name))
 
print(len(names))
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
  
  print(flist[0])
  flist[0]['templates'] = flist[0]['templates'].sum('pt',overflow='allnan')  
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
  
print('hists',flist[0])
xs = {}
with open('../data/xsec.json', 'r') as f:
   xs = json.load(f)

#flist[0]['templates']   
scale1fb = {k: xs[k] * 1000 / w for k, w in flist[0]['sumw'].items()}

try:
  flist[0]['templates'].scale(scale1fb, 'dataset')
except:
  flist[0]['jet_kin'].scale(scale1fb, 'dataset')

print('cutflow pre-scale')
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(flist[0]['cutflow_ttbar_muoncontrol'])
for proc,flow in flist[0]['cutflow_ttbar_muoncontrol'].items():
    if 'JetHT' in proc or 'SingleMuon' in proc: continue
    for cut, val in flow.items():
   
       flist[0]['cutflow_ttbar_muoncontrol'][proc][cut] *= scale1fb[proc]

print('cutflow post-scale')
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(flist[0]['cutflow_ttbar_muoncontrol'])

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

