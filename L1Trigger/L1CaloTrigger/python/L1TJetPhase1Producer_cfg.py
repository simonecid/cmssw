import FWCore.ParameterSet.Config as cms
from L1Trigger.L1CaloTrigger.caloEtaSegmentation import caloEtaSegmentation
from math import pi
from copy import deepcopy

import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process = cms.Process("L1TJetPhase1Producer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD/GeneratePfClustersAndCandidatesFromQCD_3715975.0.root')
  #fileNames = cms.untracked.vstring('file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.0.root')
)

process.TFileService = cms.Service('TFileService', fileName = cms.string("Histograms.root"))
#process.source.skipEvents = cms.untracked.uint32(1)

ak4JetFromPfClustersParameters = deepcopy(GenJetParameters)
ak4JetFromPfCandidatesParameters = deepcopy(GenJetParameters)

#GenJetParameters.src = cms.InputTag("genParticles", "", "HLT")

ak4JetFromPfClustersParameters.src = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN")
process.ak4GenJetFromPfClusters = cms.EDProducer(
    "FastjetJetProducer",
    ak4JetFromPfClustersParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    jetCollInstanceName = cms.string("ak4GenJetFromPfClusters")

)

ak4JetFromPfCandidatesParameters.src = cms.InputTag("l1pfProducer", "PF", "IN")
process.ak4GenJetFromPfCandidates = cms.EDProducer(
    "FastjetJetProducer",
    ak4JetFromPfCandidatesParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    jetCollInstanceName = cms.string("ak4GenJetFromPfCandidates")
)

process.ConvertGenJetToL1TJet = cms.EDProducer('ConvertGenJetToL1TJet',
  ak4GenJetFromPfClustersCollectionTag = cms.InputTag("ak4GenJetFromPfClusters", "ak4GenJetFromPfClusters", "L1TJetPhase1Producer"),
  ak4GenJetFromPfCandidatesCollectionTag = cms.InputTag("ak4GenJetFromPfCandidates", "ak4GenJetFromPfCandidates", "L1TJetPhase1Producer"),
)

process.L1TJetPhase1Producer = cms.EDProducer('L1TJetPhase1Producer',
  pfCandidateCollectionTag = ak4JetFromPfCandidatesParameters.src,
  pfClusterCollectionTag = ak4JetFromPfClustersParameters.src,
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

process.MatchAK4GenJetWithPhase1L1TJetFromPfClusters = cms.EDAnalyzer('MatchGenJetToL1Jet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  l1tJetCollectionTag = cms.InputTag("L1TJetPhase1Producer", "Phase1L1TJetFromPfClusters", "L1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithPhase1L1TJetFromPfCandidates = cms.EDAnalyzer('MatchGenJetToL1Jet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  l1tJetCollectionTag = cms.InputTag("L1TJetPhase1Producer", "Phase1L1TJetFromPfCandidates", "L1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithAK4JetFromPfClusters = cms.EDAnalyzer('MatchGenJetToL1Jet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  l1tJetCollectionTag = cms.InputTag("ConvertGenJetToL1TJet", "ak4L1TJetFromPfClusters", "L1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithAK4JetFromPfCandidates = cms.EDAnalyzer('MatchGenJetToL1Jet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  l1tJetCollectionTag = cms.InputTag("ConvertGenJetToL1TJet", "ak4L1TJetFromPfCandidates", "L1TJetPhase1Producer"),
)

#process.p = cms.Path(process.ak4JetFromPfClusters + process.ak4JetFromPfCandidates + process.L1TJetPhase1Producer + process.MatchGenJetToL1Jet)
process.p = cms.Path(
  process.ak4GenJetFromPfClusters +
  process.ak4GenJetFromPfCandidates +
  process.ConvertGenJetToL1TJet +
  process.L1TJetPhase1Producer +
  process.MatchAK4GenJetWithPhase1L1TJetFromPfClusters +
  process.MatchAK4GenJetWithPhase1L1TJetFromPfCandidates +
  process.MatchAK4GenJetWithAK4JetFromPfClusters +
  process.MatchAK4GenJetWithAK4JetFromPfCandidates

)
#process.p = cms.Path(process.FastjetJetProducer + process.L1TJetPhase1Producer + process.PrintMomentum)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring(
    "keep *",
    "drop *_ak4GenJetFromPfClusters_*_*",
    "drop *_ak4GenJetFromPfCandidates_*_*",    
  ),
)

process.e = cms.EndPath(process.out)
