import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '130X_dataRun3_v2'
process.GlobalTag.globaltag = THISGLOBALTAG

#--------------------- Report and output ---------------------------
# Note: in grid runs this parameter is not used.
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TFileService=cms.Service("TFileService",
                                 #fileName=cms.string('test.root'),
                                 fileName=cms.string(THISROOTFILE),
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


# Added 'vertexRef().isNonnull() &&' check for 80X data compatibility. Juska
process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('vertexRef().isNonnull() && fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets


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

##-------------------- Define the source  ----------------------------



# Note: for grid running it does not matter what's here, as input data is
# handled separately there.

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/00DC6D7F-C300-504F-A3BC-497C7D8FDA98.root')
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2023C/JetMET0/MINIAOD/22Sep2023_v2-v1/2540000/01c06858-67d2-40e6-9ad3-9af02f98f951.root') 
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2022C/JetMET/MINIAOD/22Sep2023-v1/2530000/08f285fd-506b-40de-a0e2-28a611411ff4.root')  
)

##-------------------- User analyzer  --------------------------------
   

#UpdatePuppiTuneV15(process,False)

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData           = cms.bool(True),

  ## JETS/MET ########################################
  jetsAK4          = cms.InputTag('slimmedJets'), 
  jetsAK8          = cms.InputTag('slimmedJetsPuppi'),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  
  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'), # Updated from untracked to 80X by Juska
  ptHat            = cms.untracked.InputTag('generator'), # Why do these need to be 'untracked' anyway?
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),  

  TriggerResultsTag	= cms.InputTag('TriggerResults','','HLT'),
  NoiseFilterResultsTag	= cms.InputTag('TriggerResults','','RECO'),
  l1GtSrc               = cms.InputTag('gtStage2Digis','','PAT'),  #PAT for 19Dec2023, RECO for 22Sep2023 (do edmDumpEventContent to check) 
  UnprefirableEventToken= cms.InputTag('simGtExtUnprefireable','','PAT'),
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

  ## Version Summer15_25nsV3 ( https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/ )
  # Note that it hardly matters what is put in here, as these should be overriden in analysis step anyway. Juska.
  # That's also why these JEC's are greatly dated.
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
process.p +=                      process.chs
process.p +=                      process.dijets
