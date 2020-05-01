from coffea import util,hist
import json
import os

eosdir = "/eos/uscms/store/user/jkrupa/coffea/"
indir = "Apr30_26"

os.system("mkdir -p %s" % indir)

chunk_size = 40

from os import listdir
from os.path import isfile, join
onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor")]

names = []
for name in onlyfiles:
  #if (os.path.isfile("%s%s/%s.coffea" % (eosdir,indir,name))): names.append("%s%s/%s.coffea" % (eosdir,indir,name))
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
  
print(flist[0])
  
util.save(flist[0],'%s/hists_sum.coffea' % (indir))
for i,x in enumerate(chunk_names):
  os.system("rm %s/hists_sum_%i.coffea" % (indir,i*chunk_size))

