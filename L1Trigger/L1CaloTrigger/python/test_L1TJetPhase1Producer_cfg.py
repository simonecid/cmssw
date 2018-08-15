import FWCore.ParameterSet.Config as cms
from L1Trigger.L1CaloTrigger.caloEtaSegmentation import caloEtaSegmentation
from math import pi


import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process = cms.Process("L1TJetPhase1Producer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('file:/six/sb17498/CMSSW_10_1_5/src/TTBar_300events.root')
)

process.TFileService = cms.Service('TFileService', fileName = cms.string("Histograms.root"))
#process.source.skipEvents = cms.untracked.uint32(1)

GenJetParameters.src = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN")
#GenJetParameters.src = cms.InputTag("genParticles", "", "HLT")

process.FastjetJetProducer = cms.EDProducer(
    "FastjetJetProducer",
    GenJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4)
)

process.L1TJetPhase1Producer = cms.EDProducer('L1TJetPhase1Producer',
  pfCandidateCollectionTag = cms.InputTag("l1pfProducerTightTK", "Puppi", "IN"),
  pfClusterCollectionTag = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN"),
  etaBinning = caloEtaSegmentation,
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-pi),
  phiUp = cms.double(pi),
  jetIEtaSize = cms.uint32(9),
  jetIPhiSize = cms.uint32(9),
  seedPtThreshold = cms.double(4) # GeV
)

process.PrintMomentum = cms.EDAnalyzer('PrintMomentum',
  genJetCollectionTag = cms.InputTag("FastjetJetProducer", "", "L1TJetPhase1Producer"),
  phase1L1TJetFromPfClustersTag = cms.InputTag("L1TJetPhase1Producer", "Phase1L1TJetFromPfClusters", "L1TJetPhase1Producer"),
)

process.MatchGenJetWithL1Objects = cms.EDAnalyzer('MatchGenJetWithL1Objects',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  l1tJetCollectionTag = cms.InputTag("L1TJetPhase1Producer", "Phase1L1TJetFromPfClusters", "L1TJetPhase1Producer"),
)

#process.p = cms.Path(process.FastjetJetProducer + process.L1TJetPhase1Producer + process.MatchGenJetWithL1Objects)
process.p = cms.Path(process.FastjetJetProducer + process.L1TJetPhase1Producer + process.PrintMomentum)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root')
)

process.e = cms.EndPath(process.out)
