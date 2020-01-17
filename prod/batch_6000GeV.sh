#!/bin/sh

cd /afs/cern.ch/work/z/zhixing/private/CMSSW_8_0_31/src/CMSDIJET/DijetRootTreeMaker/prod/

eval `scramv1 runtime -sh`

cmsRun flat-MC-cfg_miniAOD_6000.py
