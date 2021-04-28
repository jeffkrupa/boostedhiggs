#!/usr/bin/python

import argparse
import os
import re
import fileinput

import json
import glob

#python CoffeaSubmit.py Apr23 run_htt.py 50 1
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('settings', metavar='S', type=str, nargs='+',
                   help='label scriptname (re-tar)')
args = parser.parse_args()

if (not (len(args.settings) is 5) or (len(args.settings) is 4)):
    print("Wrong number of arguments (must be 3 or 4, found",len(args.settings),")")
    sys.exit()

label = args.settings[0]
script = args.settings[1]
files_per_job = int(args.settings[2])

loc_base = os.environ['PWD']
logdir = label
outdir = '/store/user/jkrupa/coffea_ak1/'+logdir+'/'

samplelist = {

        "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "WJetsToQQ_HT-800toInf_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "WJetsToQQ_HT400to600_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "WJetsToQQ_HT600to800_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "ZJetsToQQ_HT-800toInf_qc19_4j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "ZJetsToQQ_HT400to600_qc19_4j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "ZJetsToQQ_HT600to800_qc19_4j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016": "2016",
        "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-2016": "2016",
        "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-2016": "2016",
        "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016": "2016",
        "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016": "2016",
        "JetHT_pancakes-02_Run2016B-17Jul2018_ver2-v2-2016": "2016",
        "JetHT_pancakes-02_Run2016C-17Jul2018-v1-2016": "2016",
        "JetHT_pancakes-02_Run2016D-17Jul2018-v1-2016": "2016",
        "JetHT_pancakes-02_Run2016E-17Jul2018-v1-2016": "2016",
        "JetHT_pancakes-02_Run2016F-17Jul2018-v1-2016": "2016",
        "JetHT_pancakes-02_Run2016G-17Jul2018-v1-2016": "2016",
        "JetHT_pancakes-02_Run2016H-17Jul2018-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016B-17Jul2018_ver2-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016C-17Jul2018-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016D-17Jul2018-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016E-17Jul2018-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016F-17Jul2018-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016G-17Jul2018-v1-2016": "2016",
        "SingleMuon_pancakes-02_Run2016H-17Jul2018-v1-2016": "2016",
        # add signal, rest of ST, WJetsToLNu, DYJets

        # 2017
        #'QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8-2017': '2017',
        'QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8-2017': '2017',
        'QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8-2017': '2017',
        'QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8-2017': '2017',
        'QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8-2017': '2017',
        'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-2017' : '2017',
        'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-2017' : '2017',
        'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-2017' : '2017',
        'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-2017' : '2017',
        'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-2017' : '2017',               
        'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-2017' : '2017',
        'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8-2017' : '2017',
        'TTToHadronic_TuneCP5_13TeV-powheg-pythia8-2017' : '2017',
        'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8-2017' : '2017',
        'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017' : '2017',
        'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017' : '2017',
        'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017' : '2017',
        'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017' : '2017',
        'JetHT_pancakes-02_Run2017B-31Mar2018-v1-2017': '2017',
        'JetHT_pancakes-02_Run2017C-31Mar2018-v1-2017' : '2017',
        'JetHT_pancakes-02_Run2017D-31Mar2018-v1-2017' : '2017',
        'JetHT_pancakes-02_Run2017E-31Mar2018-v1-2017' : '2017',
        'JetHT_pancakes-02_Run2017F-31Mar2018-v1-2017' : '2017',
        'SingleMuon_pancakes-02_Run2017B-31Mar2018-v1-2017' : '2017',
        'SingleMuon_pancakes-02_Run2017C-31Mar2018-v1-2017' : '2017',
        'SingleMuon_pancakes-02_Run2017D-31Mar2018-v1-2017' : '2017',
        'SingleMuon_pancakes-02_Run2017E-31Mar2018-v1-2017' : '2017',
        'SingleMuon_pancakes-02_Run2017F-31Mar2018-v1-2017' : '2017'
        # add signal, WJetsToLNu, DYJets
}

#################################################

os.chdir('..')
os.system('xrdcp -f test/%s root://cmseos.fnal.gov//store/user/jkrupa/'%script)
if (args.settings[3]):
    #os.system('tar -vzcf dazsle_coffea.tgz . --exclude="*.root" --exclude="*.pdf" --exclude="*.pyc" --exclude="nobackup" --exclude="nobackup1" --exclude=tmp --exclude="*.tgz" --exclude="*std*" --exclude="*sum*" --exclude-vcs --exclude-caches-all --exclude="*.condor" --exclude="*.sh" --exclude=./tf/ --exclude="*tgz" ')
    os.system('tar -zvcf dazsle_coffea.tgz . --exclude="*.root" --exclude="*.pdf" --exclude="*.pyc" --exclude=tmp --exclude="*.tgz" --exclude-vcs --exclude-caches-all')
    os.system('xrdcp -f dazsle_coffea.tgz root://cmseos.fnal.gov//store/user/jkrupa/dazsle_coffea.tgz')
os.chdir(loc_base)
year=None
if (len(args.settings) == 5): year = args.settings[4]
#################################################
### Names to give to your output root files
#################################################

prefix = label

################################################
### Where is your list of root files to run over
################################################
print(label,' ',outdir)

print(str(files_per_job)+' files per job...')

#make local directory
locdir = logdir
os.system('mkdir -p  %s' %locdir)

print('CONDOR work dir: '+outdir)
#os.system('rm -rf '+outdir+label)
os.system('eos root://cmseos.fnal.gov// mkdir -p' +outdir)

totfiles = {}

#with open('../data/fileset2016.json', 'r') as f:
#    newfiles = json.load(f)
#    totfiles.update(newfiles)

#with open('../data/fileset2017.json', 'r') as f:
#    newfiles = json.load(f)
#    totfiles.update(newfiles)
with open('../data/fileset2016_preUL.json', 'r') as f:
    newfiles = json.load(f)
    totfiles.update(newfiles)

with open('../data/fileset2017_preUL.json', 'r') as f:
    newfiles = json.load(f)
    totfiles.update(newfiles)


#with open('../data/fileset2017VJets.json', 'r') as f:
#    newfiles = json.load(f)
#    totfiles.update(newfiles)

#with open('../data/fileset2017ULhadd.json', 'r') as f:
#    newfiles = json.load(f)
#    totfiles.update(newfiles)

for sample in samplelist:
    totfiles[sample] = len(totfiles[sample])

nsubmit = 0


for sample in samplelist:
    if year not in sample and year is not None: continue
    if 'QCD' in sample: continue
    prefix = sample
    print('Submitting '+prefix)
    
    #if 'JetHT' in sample: files_per_job=1
    #if 'TTToHadronic' in sample: files_per_job=10

    njobs = int(totfiles[sample]/files_per_job)+1
    remainder = totfiles[sample]-int(files_per_job*(njobs-1))

    for j in range(njobs):
        #if j>0: break
        condor_templ_file = open(loc_base+"/Coffea.templ.condor")
        sh_templ_file    = open(loc_base+"/Coffea.templ.sh")
    
        localcondor = locdir+'/'+prefix+"_"+str(j)+".condor"
        condor_file = open(localcondor,"w")
        for line in condor_templ_file:
            line=line.replace('DIRECTORY',locdir)
            line=line.replace('PREFIX',prefix)
            line=line.replace('JOBID',str(j))
            condor_file.write(line)
        condor_file.close()
    
        #copy local to eos
        #os.system('xrdcp -f %s %s' % (localcondor,eoscondor))
        #remove local copy
        #os.system('rm %s' % localcondor)
    
        localsh=locdir+'/'+prefix+"_"+str(j)+".sh"
        eosoutput="root://cmseos.fnal.gov/"+outdir+"/"+prefix+'_'+str(j)+'.coffea'
        sh_file = open(localsh,"w")
        for line in sh_templ_file:
            line=line.replace('SCRIPTNAME',script)
            line=line.replace('FILENUM',str(j))
            line=line.replace('YEAR',samplelist[sample])
            line=line.replace('SAMPLE',sample)
            line=line.replace('STARTNUM',str(j*files_per_job))
            line=line.replace('ENDNUM',str((j+1)*files_per_job))
            line=line.replace('EOSOUT',eosoutput)
            sh_file.write(line)
        sh_file.close()

        os.system('chmod u+x '+locdir+'/'+prefix+'_'+str(j)+'.sh')
        #print('condor file is: '+localcondor)
        if (os.path.exists('%s.log'  % localcondor)):
            os.system('rm %s.log' % localcondor)
        os.system('condor_submit %s' % localcondor)

        condor_templ_file.close()
        sh_templ_file.close()

        nsubmit = nsubmit + 1

print(nsubmit,"jobs submitted.")

