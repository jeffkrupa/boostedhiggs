from coffea import util,hist
import json
import os


eosdir = "/eos/uscms/store/user/jkrupa/coffea/"
indir = "6Aug20_debugQCDmuonCR_2"
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
  flist[0]['templates'] = flist[0]['templates'].sum('msd',overflow='allnan')  
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

print('hists-prescale', flist[0]) 
try:
  flist[0]['templates'].scale(scale1fb, 'dataset')
except:
  flist[0]['jet_kin'].scale(scale1fb, 'dataset')
#flist[0]['templates'].scale({k:0.63 for k,_ in flist[0]['sumw'].items() if 'QCD' in k}, 'dataset')
print('hists-postscale', flist[0]) 
util.save(flist[0],'%s/hists_sum_gru2.coffea' % (indir))
for i,x in enumerate(chunk_names):
  os.system("rm %s/hists_sum_%i.coffea" % (indir,i*chunk_size))

