import os
FR = open('flat-MC-cfg_miniAOD.py')

content=FR.read()
path = '/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/BstarToJJ/'
PATH = ['1000/submit_20180615_193226/','2000/submit_20180615_193241/','3000/submit_20180615_193257/','4000/submit_20180615_193314/','5000/submit_20180615_193329/','6000/submit_20180615_193345/','7000/submit_20180615_193401/','8000/submit_20180615_193417/','9000/submit_20180615_193433/','500/submit_20180615_193210/']

F1 = ''
F2 = ''
F3 = ''

for i in PATH:
  Rp = content
  F1 = ''
  F2 = ''
  F3 = ''
  rpath = path +i
  File =  os.listdir(rpath)
  index = 0
  for j in File :
    if not('step3' in j and 'root' in j):
      continue
    if 'step3_0' in j and i.split('/')[0]=='1000':
      continue
    if 'step3_103' in j and i.split('/')[0]=='4000':
      continue
    index += 1
    if index < 100:
      F1 += '\'file:'+rpath+j+'\',\n'
    if 200 > index > 100:
      F2 += '\'file:'+rpath+j+'\',\n'
    if index > 200:
      F3 += '\'file:'+rpath+j+'\',\n'
  Rp = Rp.replace('File111',F1).replace('File222',F2).replace('File333',F3).replace('THISROOTFILE','\'test_'+i.split('/')[0]+'.root\'')
  Rp  = Rp.replace('TROOTFILE','\'/tmp/bStar_'+i.split('/')[0]+'GeV.root\'')
  Fout=open('flat-MC-cfg_miniAOD_'+i.split('/')[0]+'.py','w+')
  Fout.write(Rp)
  Fout.close()


  con = '''#!/bin/sh

cd /afs/cern.ch/work/z/zhixing/private/CMSSW_9_4_0/src/CMSDIJET/DijetRootTreeMaker/prod/

eval `scramv1 runtime -sh`

cmsRun '''

  con += 'flat-MC-cfg_miniAOD_'+i.split('/')[0]+'.py\n\n'
  con += 'mv /tmp/bStar_'+i.split('/')[0]+'GeV.root /eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/TylerW/2017JetHT/Bstar_'+i.split('/')[0]+'GeV/'
  Fbatch = open('batch_'+i.split('/')[0]+'.sh','w+')
  Fbatch.write(con)
  Fbatch.close()
  os.system('chmod 751 batch_'+i.split('/')[0]+'.sh')
  if i.split('/')[0]==1000: 
    continue
  os.system('bsub -q 8nh -J '+i.split('/')[0]+'GeV batch_'+i.split('/')[0]+'.sh')
