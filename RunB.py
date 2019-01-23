from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = '/afs/cern.ch/work/z/zhixing/private/CMSSW_9_4_0/src/CMSDIJET/DijetRootTreeMaker/crab_output/' 
config.General.requestName = 'Oct_RunB'
config.JobType.psetName = '/afs/cern.ch/work/z/zhixing/private/CMSSW_9_4_0/src/CMSDIJET/DijetRootTreeMaker/prod/flat-data-cfg_miniAOD_B.py'
config.JobType.pluginName = 'Analysis'
config.Data.inputDataset = 'global'
config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.inputDataset = '/JetHT/Run2017B-17Nov2017-v1/MINIAOD'
config.Data.unitsPerJob = 5 #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.outputDatasetTag = 'analysis'
config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/TylerW/RunII2017' 
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
config.Site.blacklist = []
config.Site.storageSite = 'T2_CH_CERN'
