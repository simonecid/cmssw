import FWCore.ParameterSet.Config as cms

process = cms.Process("MatchGenJetWithL1Objects")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.TFileService = cms.Service("TFileService", fileName = cms.string("l1tObjectGenJetMatching.root") )

process.MatchGenJetWithL1Objects = cms.EDAnalyzer('MatchGenJetWithL1Objects',
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  l1tJetCollectionTag = cms.InputTag("simCaloStage2Digis"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
  l1tTauCollectionTag = cms.InputTag("simCaloStage2Digis"),
  l1tEGammaCollectionTag = cms.InputTag("simCaloStage2Digis"),
)

process.source = source_QCD_Pt_15_3000

process.p = cms.Path(process.MatchGenJetWithL1Objects)
