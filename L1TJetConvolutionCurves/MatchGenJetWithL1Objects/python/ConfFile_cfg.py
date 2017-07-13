import FWCore.ParameterSet.Config as cms

process = cms.Process("MatchGenJetWithL1Objects")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
  # replace 'myfile.root' with the source file you want to use
  fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/user/s/sbologna/CMSSW/CMSSW_9_0_0/src/step2.root'
  )
)

process.demo = cms.EDAnalyzer('MatchGenJetWithL1Objects',
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  l1tJetCollectionTag = cms.InputTag("hltGtStage2Digis")
)


process.p = cms.Path(process.demo)
