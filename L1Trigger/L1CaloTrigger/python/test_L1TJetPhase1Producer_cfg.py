import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TJetPhase1Producer")
#process = cms.Process("SaveEvent")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/six/sb17498/CMSSW_10_1_5/src/TTBar_300events.root')
)
#process.source.skipEvents = cms.untracked.uint32(1)

process.L1TJetPhase1Producer = cms.EDProducer('L1TJetPhase1Producer',
  pfCandidateCollectionTag = cms.InputTag("l1pfProducer", "Puppi", "IN"),
  pfClusterCollectionTag = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN")
)

process.p = cms.Path(process.L1TJetPhase1Producer)

#process.e = cms.EndPath(process.out)
