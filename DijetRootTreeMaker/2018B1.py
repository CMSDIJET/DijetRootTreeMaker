from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = '/afs/cern.ch/work/z/zhixing/private/CMSSW_10_2_5/src/CMSDIJET/DijetRootTreeMaker/crab_output/'
config.General.requestName = '2018B_Jun_v1'
config.JobType.psetName = '/afs/cern.ch/work/z/zhixing/private/CMSSW_10_2_5/src/CMSDIJET/DijetRootTreeMaker/prod/flat-data-cfg_miniAOD_B.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.Data.inputDataset = 'global'
#config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.inputDataset = '/JetHT/Run2018B-17Sep2018-v1/MINIAOD'
config.Data.unitsPerJob = 250 #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.outputDatasetTag = 'analysis'
config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/TylerW/2018JetHT'
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
config.Site.blacklist = []
config.Site.storageSite = 'T2_CH_CERN'

