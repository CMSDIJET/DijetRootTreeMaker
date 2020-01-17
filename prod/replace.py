import os
FR = open('flat-MC-cfg_miniAOD_temp.py')

content=FR.read()
#path = '/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/BstarToJJ/'
#PATH =['1000/submit_20170925_100302/', '2000/submit_20170925_100318/', '3000/submit_20170925_100333/', '4000/submit_20170925_100350/', '5000/submit_20170925_100407/', '6000/submit_20170925_100423/', '7000/submit_20170925_100439/', '8000/submit_20170925_100453/', '9000/submit_20170925_100509/']
path = '/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/TylerW/Generation/MinBias/'

PATH = ['RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_miniAODv2_Tyler1000GeV_2MINIAOD/181126_041331/0000/','RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_miniAODv2_Tyler2000GeV_2MINIAOD/181126_041350/0000/','RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_miniAODv2_Tyler3000GeV_2MINIAOD/181126_041409/0000/','RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_miniAODv2_Tyler4000GeV_2MINIAOD/181126_041427/0000/']

for i in PATH:
  Rp = content
  F1 = ''
  F2 = ''
  F3 = ''
  rpath = path +i
  File =  os.listdir(rpath)
  index = 0
  for j in File :
    if not('RSGravitonToQuarkQuark' in j and 'root' in j):
      continue
    index += 1
    if index < 200:
      F1 += '\'file:'+rpath+j+'\',\n'
    if 200 > index > 400:
      F2 += '\'file:'+rpath+j+'\',\n'
    if index > 600:
      F3 += '\'file:'+rpath+j+'\',\n'
  Rp = Rp.replace('File111',F1).replace('File222',F2).replace('File333',F3).replace('THISROOTFILE','\'test_'+i.split('/')[0].split('_')[-2]+'.root\'')
  Fout=open('flat-MC-cfg_miniAOD_'+i.split('/')[0]+'.py','w+')
  Fout.write(Rp)
  Fout.close()
