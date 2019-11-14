import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = THISGLOBALTAG
process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'


#--------------------- Report and output ---------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService=cms.Service("TFileService",
                                 #fileName=cms.string(THISROOTFILE),
                                 fileName=cms.string('RSG2bb_9000.root'),
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


# ----------------------- Jet Tool Box  -----------------
# ----- giulia test: do not recluster ak4 and ca8 jets to save time --------

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

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



process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:QstarToJJ_M_4000_TuneCUETP8M1_13TeV_pythia8__MINIAODSIM__Asympt50ns_MCRUN2_74_V9A-v1__70000__AA35D1E7-FEFE-E411-B1C5-0025905B858A.root')    
    #fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/QstarToJJ_M_1000_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/00F85752-BCFB-E411-A29A-000F5327349C.root')
    #fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0E4CEBFE-ECFB-E411-9F0C-842B2B29273C.root')
    fileNames = cms.untracked.vstring(
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_10.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_11.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_12.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_13.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_14.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_15.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_16.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_17.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_18.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_19.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_1.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_20.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_21.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_22.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_23.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_24.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_25.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_26.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_27.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_28.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_29.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_2.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_30.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_31.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_32.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_33.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_34.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_35.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_36.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_37.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_38.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_39.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_3.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_40.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_41.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_42.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_43.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_44.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_45.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_46.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_47.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_48.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_49.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_50.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_51.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_52.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_53.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_54.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_55.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_56.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_57.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_58.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_59.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_5.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_60.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_61.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_62.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_63.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_64.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_65.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_66.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_67.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_68.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_69.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_6.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_70.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_71.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_72.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_73.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_74.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_75.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_76.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_77.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_78.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_79.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_80.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_81.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_82.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_83.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_84.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_85.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_86.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_87.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_88.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_89.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_8.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_90.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_91.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_92.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_93.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_94.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_95.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_96.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_97.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_98.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_99.root',
'file:/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/deguio/sampleGen_2018/RSGravitonToBBbar/9000/submit_20190508_181049/step3_9.root',

)
#    fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v11_ext1-v1/70000/041AA818-0820-E811-8B46-02163E01A0A4.root')

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
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),

  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'),
  ptHat            = cms.untracked.InputTag('generator'),
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 


  ## trigger ###################################
  triggerAlias     = cms.vstring('PFHT900','PFHT1050','PFHT600','PFHT350'
                                 ,'PFHT650MJJ950','PFHT650MJJ900'
                                 ,'PFJET500','PFJET450','PFJET200'),
  triggerSelection = cms.vstring(
     'HLT_PFHT900_v*',
     'HLT_PFHT1050_v*',
     'HLT_PFHT600_v*',
     'HLT_PFHT350_v*',
     'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*',
     'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*',
     'HLT_PFJet500_v*',
     'HLT_PFJet450_v*',
     'HLT_PFJet200_v*',
     
  ),
  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    #l1tIgnoreMask         = cms.bool(False),
    #l1techIgnorePrescales = cms.bool(False),
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
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Autumn18_V8_MC/Autumn18_V8_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Autumn18_V8_MC/Autumn18_V8_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Autumn18_V8_MC/Autumn18_V8_MC_L3Absolute_AK4PFchs.txt'),
)


# ------------------ path --------------------------

process.p = cms.Path()

#process.p +=                     process.prunedGenParticlesDijet
process.p +=                     process.chs 
process.p +=                     process.dijets 
