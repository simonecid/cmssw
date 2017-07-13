import FWCore.ParameterSet.Config as cms

process = cms.Process("MatchGenJetWithL1Objects")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.TFileService = cms.Service("TFileService", fileName = cms.string("l1tObjectGenJetMatching.root") )

#process.source = cms.Source("PoolSource",
#  # replace 'myfile.root' with the source file you want to use
#  fileNames = cms.untracked.vstring(
#    'file:/afs/cern.ch/user/s/sbologna/CMSSW/CMSSW_9_0_0/src/step2.root'
#  )
#)
readFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [
  '/store/relval/CMSSW_9_0_0/RelValADDMonoJet_d3MD3_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_cc7_rsb-v1/00000/0421E5E3-5C26-E711-A4FC-001E67457E36.root',
  '/store/relval/CMSSW_9_0_0/RelValADDMonoJet_d3MD3_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_cc7_rsb-v1/00000/04EB2B69-5826-E711-9302-A4BF0103375E.root',
  '/store/relval/CMSSW_9_0_0/RelValADDMonoJet_d3MD3_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_cc7_rsb-v1/00000/1044C16D-5926-E711-BA0C-0242AC130002.root',
  '/store/relval/CMSSW_9_0_0/RelValADDMonoJet_d3MD3_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_cc7_rsb-v1/00000/56BB0CED-5826-E711-8436-0242AC130006.root',
  '/store/relval/CMSSW_9_0_0/RelValADDMonoJet_d3MD3_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_cc7_rsb-v1/00000/E2B0E749-5E26-E711-87FE-A4BF0103325A.root',
  '/store/relval/CMSSW_9_0_0/RelValADDMonoJet_d3MD3_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_cc7_rsb-v1/00000/F65D0477-5826-E711-907E-001E674FAF23.root' 
] )

process.source = source

process.demo = cms.EDAnalyzer('MatchGenJetWithL1Objects',
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  l1tJetCollectionTag = cms.InputTag("simCaloStage2Digis"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
  l1tTauCollectionTag = cms.InputTag("simCaloStage2Digis"),
  l1tEGammaCollectionTag = cms.InputTag("simCaloStage2Digis"),
)


process.p = cms.Path(process.demo)
