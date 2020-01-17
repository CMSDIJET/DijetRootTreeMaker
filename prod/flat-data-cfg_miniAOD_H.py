import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')



## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
#process.GlobalTag.globaltag = THISGLOBALTAG

#--------------------- Report and output ---------------------------
# Note: in grid runs this parameter is not used.
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string('test.root'),
                                 #fileName=cms.string(THISROOTFILE),
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
                                                                                                                                     
                                                                       ])                                                                                           
                               )


# Added 'vertexRef().isNonnull() &&' check for 80X data compatibility. Juska
process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('vertexRef().isNonnull() && fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
bTagDiscriminators = [
         'softPFMuonBJetTags',
         'softPFElectronBJetTags',
         'pfJetBProbabilityBJetTags',
         'pfJetProbabilityBJetTags',
         'pfCombinedInclusiveSecondaryVertexV2BJetTags',
         'pfDeepCSVJetTags:probudsg',
         'pfDeepCSVJetTags:probb',
         'pfDeepCSVJetTags:probc',
         'pfDeepCSVJetTags:probbb',
         'pfDeepCSVJetTags:probcc',
         'pfCombinedMVAV2BJetTags',
     ]
bTagInfos = [
        'pfImpactParameterTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'pfDeepCSVTagInfos',
     ]
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
updateJetCollection(
        process,
        labelName = "DeepFlavour",
        jetSource = cms.InputTag('slimmedJets'),  # 'ak4Jets'
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

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
    #fileNames = cms.untracked.vstring("root://eoscms//eos/cms/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v2/000/273/411/00000/10CB3C59-721B-E611-AFB4-02163E012711.root")
    #fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/j/juska/eos/cms/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v2/000/273/411/00000/10CB3C59-721B-E611-AFB4-02163E012711.root")
    #fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/j/juska/eos/cms/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v2/000/273/730/00000/EA345ED4-B821-E611-BEA5-02163E0138E2.root")
    fileNames = cms.untracked.vstring(
#'/store/data/Run2016B/JetHT/MINIAOD/03Feb2017_ver1-v1/80000/F8F61751-A5EC-E611-92B0-02163E019D1E.root'
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/131/00000/3AEDC011-9D80-E611-915C-FA163EDE598F.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/189/00000/12C2DA56-0981-E611-88BE-FA163E298B60.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/167/00000/C29ED806-E380-E611-B40C-FA163E4FF519.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/085/00000/7E4416FF-887F-E611-9882-02163E0134A5.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/010/00000/2CEFC128-757E-E611-A6CB-02163E014679.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/130/00000/A4836AAA-7880-E611-9974-02163E01472A.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/117/00000/5CACF77F-0D80-E611-9F04-02163E01430B.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/156/00000/F4E4C8A6-D080-E611-95E8-FA163E8C4370.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/145/00000/38B8A77F-BB80-E611-8A23-02163E013579.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/178/00000/8E2B831A-F980-E611-9270-02163E011AEE.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/090/00000/CE81D20C-C37F-E611-ADA9-02163E01431B.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/165/00000/8C355254-D480-E611-BAF0-02163E0145B2.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/166/00000/CC52562E-DA80-E611-B5F0-FA163E720510.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/202/00000/20E08234-6B81-E611-A504-FA163E60DE95.root',
'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v1/000/281/201/00000/AC5F9427-6A81-E611-A899-02163E01420B.root'
)
    
)

##-------------------- User analyzer  --------------------------------


# Residue from AOD and RECO running
calo_collection=''
cluster_collection=''
pfcalo_collection=''
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData          = cms.bool(True),

  ## JETS/MET ########################################
  jets		   = cms.InputTag('selectedUpdatedPatJetsDeepFlavour'),
  jetsAK4             = cms.InputTag('slimmedJets'), 
  #jetsAK8             = cms.InputTag('slimmedJetsAK8'),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
 # ptMinAK8         = cms.double(10),
  
  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'), # Updated from untracked to 80X by Juska
  ptHat            = cms.untracked.InputTag('generator'), # Why do these need to be 'untracked' anyway?
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
 # genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),  


  ## trigger ###################################
  #triggerAlias     = cms.vstring('Fat','PFHT650','PFNoPUHT650','HT750','HT550'),
  ##### For 0T data  #####
  #triggerAlias     = cms.vstring('L1Jet68','L1Jet36','L1Jet16','L1EG20','L1EG5'),
  ##### For JetHT PD ##### 
  triggerAlias     = cms.vstring('PFHT780','PFHT890','PFHT1050',
                                 'PFJET400','PFJET450','PFJET500','PFJET550',
                                 'Mu50',
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
     'HLT_PFJet400_v*', # it exists and it is prescaled
     'HLT_PFJet450_v*', # it exists and it is unprescaled
     'HLT_PFJet500_v*', # it exists and it is unprescaled
     'HLT_PFJet550_v*', # it exists and it is unprescaled
     'HLT_Mu50_v*', 	# it exists and it is unprescaled
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

  noiseFilterConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','RECO'), #for prompt reco
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    #l1tIgnoreMask         = cms.bool(False),
    #l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),


  ## JECs ################
  redoJECs  = cms.bool(True),

  ## Version Summer15_25nsV3 ( https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/ )
  # Note that it hardly matters what is put in here, as these should be overriden in analysis step anyway. Juska.
  # That's also why these JEC's are greatly dated.
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK4PFchs.txt'),



)


# ------------------ path --------------------------

process.p = cms.Path()
#process.p += process.patJetCorrFactorsDeepFlavour*
process.p += process.updatedPatJetsDeepFlavour
#process.p += process.patJetCorrFactorsTransientCorrectedDeepFlavour*
#process.p += rocess.softPFMuonsTagInfosDeepFlavour*
#process.p += process.softPFMuonBJetTagsDeepFlavour*
#process.p += process.softPFElectronsTagInfosDeepFlavour*
#process.p += process.softPFElectronBJetTagsDeepFlavour*
#process.p += process.pfImpactParameterTagInfosDeepFlavour*
#process.p += process.pfJetBProbabilityBJetTagsDeepFlavour*
#process.p += process.pfJetProbabilityBJetTagsDeepFlavour*
#process.p += process.pfInclusiveSecondaryVertexFinderTagInfosDeepFlavour*
#process.p += process.pfCombinedInclusiveSecondaryVertexV2BJetTagsDeepFlavour*
#process.p += process.pfDeepCSVTagInfosDeepFlavour*
#process.p += process.pfDeepCSVJetTagsDeepFlavour*
#process.p += process.pfSecondaryVertexTagInfosDeepFlavour*
#process.p += process.pfCombinedMVAV2BJetTagsDeepFlavour*
#process.p += process.updatedPatJetsTransientCorrectedDeepFlavour*
#process.p += process.selectedUpdatedPatJetsDeepFlavour*
process.p += process.chs
process.p += process.dijets
