import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# fileList = FileUtils.loadListFromFile('ttbar.list')
# readFiles = cms.untracked.vstring(*fileList)

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring("file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/QCD_Pt-15To7000_PU200/inputs104X_71.root"),
  #fileNames = cms.untracked.vstring(
  #  "file:pf500.root",
  #)
)

process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetCalibrator_*_*",
    "keep *_Phase1L1TSumsProducer_*_*",
  ),
)

# process.SaveSums = cms.EDProducer("SaveGenSumsAndL1Sums",
#   genMETCollectionTag = cms.InputTag("", "", "")
#   l1tMETCollectionTag = cms.InputTag("", "", "")
#   genJetCollectionTag = cms.InputTag("", "", "")
#   l1tHTCollectionTag = cms.InputTag("", "", "")
# )

process.p = cms.Path(process.Phase1L1TJetsSequence )

process.e = cms.EndPath(process.out)


# process.out = cms.OutputModule("PoolOutputModule",
#   fileName = cms.untracked.string('myOutputFile.root'),
#   outputCommands = cms.untracked.vstring(
#     "drop *",
#     "keep *_Phase1L1TJetProducer_*_*",
#     "keep *_ak4GenJetsNoNu_*_*",
#     "keep *_Phase1L1TJetCalibrator_*_*",
#   ),
# )

# from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer

# Phase1L1TJetsSequence = cms.Sequence(
#   Phase1L1TJetProducer
# )

# process.p = cms.Path(process.Phase1L1TJetsSequence)

# process.e = cms.EndPath(process.out)