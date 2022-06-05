import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = THISGLOBALTAG
process.GlobalTag.globaltag = '122X_mcRun3_2021_realistic_v9'


#--------------------- Report and output ---------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService=cms.Service("TFileService",
                                 #fileName=cms.string(THISROOTFILE),
                                 fileName=cms.string("test.root"),
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
    #fileNames = cms.untracked.vstring('file:QstarToJJ_M_4000_TuneCUETP8M1_13TeV_pythia8__MINIAODSIM__Asympt50ns_MCRUN2_74_V9A-v1__70000__AA35D1E7-FEFE-E411-B1C5-0025905B858A.root')    
    #fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/QstarToJJ_M_1000_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/00F85752-BCFB-E411-A29A-000F5327349C.root')
    #fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0E4CEBFE-ECFB-E411-9F0C-842B2B29273C.root')
    #fileNames = cms.untracked.vstring('/store/mc/RunIIFall18MiniAOD/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v12-v2/100000/1DF8D896-EFA6-5440-B9DE-B9F2787BBF01.root')
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Winter22MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13p6TeV_pythia8/MINIAODSIM/122X_mcRun3_2021_realistic_v9-v2/2430000/004af1ea-531a-40ad-8859-bb0493b5c2fa.root')
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
  ptMinAK4         = cms.double(0),#edw
  ptMinAK8         = cms.double(0),#edw

  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'),
  ptHat            = cms.untracked.InputTag('generator'),
  genParticles     = cms.InputTag('prunedGenParticles'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),     


  ## trigger ################################### NOTE!!!!add PFHT1050 instead of PFHT650
  triggerAlias     = cms.vstring('PFHT780','PFHT890','PFHT1050',
                                 'PFJET200','PFJET400','PFJET450','PFJET500','PFJET550',
                                 'Mu50',
                                 'AK8PFJet320', 'AK8PFJet400','AK8PFJet450','AK8PFJet500','AK8PFJet550',
		                 'CaloJet500NoJetID','CaloJet550NoJetID',                               'HLT_PFHT500_PFMET100_PFMHT100_IDTight','HLT_PFHT500_PFMET110_PFMHT110_IDTight','HLT_PFHT700_PFMET85_PFMHT85_IDTight','HLT_PFHT700_PFMET95_PFMHT95_IDTight',
'HLT_PFHT800_PFMET75_PFMHT75_IDTight','HLT_PFHT800_PFMET85_PFMHT85_IDTight',
'PFHT380_SixJet32_DoubleBTagCSV_p075', 
'HLT_PFHT430_SixJet40_BTagCSV_p080'),                                 
#

  triggerSelection = cms.vstring(

     ###
     ### For JetHT PD ###
     ###
     'HLT_PFHT780_v*',  # it exists and it is the PFHT prescaled 
     'HLT_PFHT890_v*',  # it exists and it is the PFHT prescaled
     'HLT_PFHT1050_v*', # it exists and it is the PFHT unprescaled
     'HLT_PFJet200_v*',
     'HLT_PFJet400_v*', # it exists and it is prescaled
     'HLT_PFJet450_v*', # it exists and it is unprescaled
     'HLT_PFJet500_v*', # it exists and it is unprescaled
     'HLT_PFJet550_v*', # it exists and it is unprescaled
     'HLT_Mu50_v*', 	# it exists and it is unprescaled
     'HLT_AK8PFJet320_v*',# it exists and it is prerescaled
     'HLT_AK8PFJet400_v*',# it exists and it is unprescaled  
     'HLT_AK8PFJet450_v*',# it exists and it is unprescaled
     'HLT_AK8PFJet500_v*',# it exists and it is unprescaled
     'HLT_AK8PFJet550_v*',# it exists and it is unprescaled
     'HLT_CaloJet500_NoJetID_v*', # it exists and it is unprescaled
     'HLT_CaloJet550_NoJetID_v*', # it exists and it is unprescaled
     'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*', # it exists and it is unprescaled
     'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v*', # it exists and it is unprescaled 
     'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*', # it exists and it is unprescaled
     'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v*', # it exists and it is unprescaled
     'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*', # it exists and it is unprescaled
     'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v*', # it exists and it is unprescaled
     'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', # it exists and it is unprescaled
     'HLT_PFHT430_SixJet40_BTagCSV_p080_v*',       # it exists and it is unprescaled
     ###
  ),
  triggerConfiguration = cms.PSet(
    usePathStatus = cms.bool( True ),
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    #l1tIgnoreMask         = cms.bool(False),
   # l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),


  ## Noise Filters ###################################
  noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
  noiseFilterSelection_globalSuperTightHalo2016Filter = cms.string('Flag_globalSuperTightHalo2016Filter'),
  noiseFilterSelection_HBHENoiseIsoFilter = cms.string('Flag_HBHENoiseIsoFilter'),
  noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
  noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
  noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
  noiseFilterSelection_BadChargedCandidateFilter = cms.string('Flag_BadChargedCandidateFilter'),
  noiseFilterSelection_BadPFMuonFilter = cms.string('Flag_BadPFMuonFilter'),
  #noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
  #noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
  #_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
  #noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
  #noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
  #noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
  #noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
  #noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
  #noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
  # and the sub-filters
  #noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
  #noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
  #noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),

  noiseFilterConfiguration = cms.PSet(
    usePathStatus = cms.bool( True ),
    hltResults            = cms.InputTag('TriggerResults','','PAT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
  l1tIgnoreMaskAndPrescale = cms.bool(False),
    #l1tIgnoreMask         = cms.bool(False),
    #l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),


  ## JECs ################
  redoJECs  = cms.bool(True),

  ## Version Summer15_25nsV3
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK8PFchs.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK8PFchs.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK4PFchs.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK8PFchs.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK8PFchs.txt')
)


# ------------------ path --------------------------

process.p = cms.Path()

process.p +=                     process.chs                      
process.p +=                     process.dijets 
