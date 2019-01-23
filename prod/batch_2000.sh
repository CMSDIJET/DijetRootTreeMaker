#!/bin/sh

cd /afs/cern.ch/work/z/zhixing/private/CMSSW_9_4_0/src/CMSDIJET/DijetRootTreeMaker/prod/

eval `scramv1 runtime -sh`

cmsRun flat-MC-cfg_miniAOD_2000.py

mv /tmp/bStar_2000GeV.root /eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/TylerW/2017JetHT/Bstar_2000GeV/