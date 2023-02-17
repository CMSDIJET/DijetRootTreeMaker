import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = THISGLOBALTAG
#process.GlobalTag.globaltag = '106X_upgrade2018_realistic_v16_L1v1'


#--------------------- Report and output ---------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string(THISROOTFILE),
                                 #fileName=cms.string("test.root"),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

## --- suppress long output ---> wantSummary = cms.untracked.bool(False) 

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

############## output  edm format ###############
process.out = cms.OutputModule('PoolOutputModule',                                                                                                                  
                               fileName = cms.untracked.string('jettoolbox.root'),                                                                              
                               outputCommands = cms.untracked.vstring([
                                                                      'keep *_slimmedJets_*_*',                                                                  
                                                                      'keep *_slimmedJetsPuppi_*_*',                                                                  
                                                                       ])                                                                                           
                               )

# ----------------------- Jet Tool Box  -----------------
# ----- giulia test: do not recluster ak4 and ca8 jets to save time --------

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.slimmedGenJetsAK8 = ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)


#-------------------------------------------------------
# Gen Particles Pruner
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.prunedGenParticlesDijet = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
    "drop  *  ", # by default
    "keep ( status = 3 || (status>=21 && status<=29) )", # keep hard process particles
    )
)


#------------- Recluster Gen Jets to access the constituents -------
#already in toolbox, just add keep statements

process.out.outputCommands.append("keep *_slimmedGenJets_*_*")
process.out.outputCommands.append("keep *_slimmedGenJetsAK8_*_*")

##-------------------- Define the source  ----------------------------



process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Winter22MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13p6TeV_pythia8/MINIAODSIM/122X_mcRun3_2021_realistic_v9-v2/2430000/01796a5d-9b6a-46fc-a36f-150cb43af911.root')
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2520000/02DA51E2-0619-C045-877E-7235570345C6.root')
    )




##-------------------- User analyzer  --------------------------------

#Residue from deleted reco and AOD sequences
calo_collection=''
cluster_collection=''
pfcalo_collection=''
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',
  
  # There's no avoiding this in Consumes era
  isData          = cms.bool(False),
  
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsPuppi'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(0),
  ptMinAK8         = cms.double(0),

  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'),
  ptHat            = cms.untracked.InputTag('generator'),
  genParticles     = cms.InputTag('prunedGenParticles'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),   

  TriggerResultsTag	= cms.InputTag('TriggerResults','','HLT'),
  NoiseFilterResultsTag	= cms.InputTag('TriggerResults','','PAT'),
  l1tResults            = cms.InputTag(''),
  daqPartitions         = cms.uint32(1),
  l1tIgnoreMaskAndPrescale = cms.bool(False),
  throw                 = cms.bool(False),  


  triggerSelection = cms.vstring(

     ###
     #ATTENTION: The order below is NOT the order that they will appear in the "triggerResult" vector in the output tree.
     ###
     'HLT_PFHT780_v',  # it exists and it is the PFHT prescaled 
     'HLT_PFHT1050_v',  # it exists and it is the PFHT prescaled
     'HLT_PFHT890_v', # it exists and it is the PFHT unprescaled
     'HLT_PFJet400_v', # it exists and it is prescaled
     'HLT_PFJet450_v', # it exists and it is unprescaled
     'HLT_PFJet500_v', # it exists and it is unprescaled
     'HLT_PFJet550_v', # it exists and it is unprescaled
     'HLT_Mu50_v', 	# it exists and it is unprescaled
     'HLT_CaloJet500_NoJetID_v', # it exists and it is unprescaled
     'HLT_CaloJet550_NoJetID_v', # it exists and it is unprescaled
     'HLT_AK8PFJet320_v',       # it exists and it is unprescaled
     'HLT_AK8PFJet400_v',       # it exists and it is unprescaled
     'HLT_AK8PFJet450_v',       # it exists and it is unprescaled
     'HLT_AK8PFJet500_v',       # it exists and it is unprescaled
     'HLT_AK8PFJet550_v',       # it exists and it is unprescaled
     ###
  ),


  ## JECs ################
  redoJECs  = cms.bool(True),

  ## Version Summer15_25nsV3
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_RunCD_V1_DATA_L2L3Residual_AK4PFPuppi.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L2Relative_AK4PFPuppi.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L3Absolute_AK4PFPuppi.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_RunCD_V1_DATA_L2L3Residual_AK4PFPuppi.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L3Absolute_AK4PFchs.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L2Relative_AK4PFPuppi.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Winter22Run3_V1/Winter22Run3_V1_L3Absolute_AK4PFPuppi.txt')
)


# ------------------ path --------------------------

process.p = cms.Path()

process.p +=                     process.chs                      
process.p +=                     process.dijets 
