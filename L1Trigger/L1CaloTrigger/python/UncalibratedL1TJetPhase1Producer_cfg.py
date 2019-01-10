import FWCore.ParameterSet.Config as cms
from L1Trigger.L1CaloTrigger.caloEtaSegmentation import caloEtaSegmentation
from JetCalibration.ApplyCalibrationFactors.BuildCalibrationVPSet import BuildCalibrationVPSet
from math import pi
from copy import deepcopy

import pickle as pkl

import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process = cms.Process("UncalibratedL1TJetPhase1Producer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = process.source = cms.Source("PoolSource",
#  fileNames = cms.untracked.vstring('file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD/GeneratePfClustersAndCandidatesFromQCD_3715975.0.root')
#)

process.source = process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.0.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.100.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.101.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.102.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.103.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.104.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.105.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/GeneratePfClustersAndCandidatesFromQCD_PU200/GeneratePfClustersAndCandidatesFromQCD_PU200_3720957.106.root",
  )
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

process.ConvertGenJetToCaloJet = cms.EDProducer('ConvertGenJetToCaloJet',
  ak4GenJetFromPfClustersCollectionTag = cms.InputTag("ak4GenJetFromPfClusters", "ak4GenJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
  ak4GenJetFromPfCandidatesCollectionTag = cms.InputTag("ak4GenJetFromPfCandidates", "ak4GenJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
)

process.Phase1L1TJetFromPfCandidatesProducer = cms.EDProducer('L1TJetPhase1Producer',
  inputCollectionTag = ak4JetFromPfCandidatesParameters.src,
  etaBinning = caloEtaSegmentation,
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-pi),
  phiUp = cms.double(pi),
  jetIEtaSize = cms.uint32(9),
  jetIPhiSize = cms.uint32(9),
  seedPtThreshold = cms.double(4), # GeV
  puSubtraction = cms.bool(True),
  outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates"),
  vetoZeroPt = cms.bool(True)
)

process.Phase1L1TJetFromPfClustersProducer = cms.EDProducer('L1TJetPhase1Producer',
  inputCollectionTag = ak4JetFromPfClustersParameters.src,
  etaBinning = caloEtaSegmentation,
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-pi),
  phiUp = cms.double(pi),
  jetIEtaSize = cms.uint32(9),
  jetIPhiSize = cms.uint32(9),
  seedPtThreshold = cms.double(4), # GeV
  puSubtraction = cms.bool(True),
  outputCollectionName = cms.string("Phase1L1TJetFromPfClusters"),
  vetoZeroPt = cms.bool(True)
)

with open('JetCalibration/Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_PU200_PUSubtraction_NoZeroPtJets.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters = pkl.load(f)

process.CalibratePhase1L1TJetFromPfClusters = cms.EDProducer('ApplyCalibrationFactors',
  inputCollectionTag = cms.InputTag("Phase1L1TJetFromPfClustersProducer", "Phase1L1TJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
  absEtaBinning = cms.vdouble(0, 1.4, 3, 6),
  calibration = Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters,
  outputCollectionName = cms.string("CalibratedPhase1L1TJetFromPfClusters")
)

process.PrintMomentum = cms.EDAnalyzer('PrintMomentum',
  genJetCollectionTag = cms.InputTag("FastjetJetProducer", "", "UncalibratedL1TJetPhase1Producer"),
  phase1L1TJetFromPfClustersTag = cms.InputTag("L1TJetPhase1Producer", "Phase1L1TJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithPhase1L1TJetFromPfClusters = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  #caloJetCollectionTag = cms.InputTag("Phase1L1TJetFromPfClustersProducer", "Phase1L1TJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
  caloJetCollectionTag = cms.InputTag("CalibratePhase1L1TJetFromPfClusters", "CalibratedPhase1L1TJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithPhase1L1TJetFromPfCandidates = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("Phase1L1TJetFromPfCandidatesProducer", "Phase1L1TJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithAK4JetFromPfClusters = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("ConvertGenJetToCaloJet", "ak4CaloJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithAK4JetFromPfCandidates = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("ConvertGenJetToCaloJet", "ak4CaloJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
)

#process.p = cms.Path(process.ak4JetFromPfClusters + process.ak4JetFromPfCandidates + process.L1TJetPhase1Producer + process.MatchGenJetToRecoCaloJet)

process.p = cms.Path(
  process.ak4GenJetFromPfClusters +
  process.ak4GenJetFromPfCandidates +
  process.ConvertGenJetToCaloJet +
  process.Phase1L1TJetFromPfCandidatesProducer +
  process.Phase1L1TJetFromPfClustersProducer +
  process.CalibratePhase1L1TJetFromPfClusters +
  process.MatchAK4GenJetWithPhase1L1TJetFromPfCandidates +
  process.MatchAK4GenJetWithPhase1L1TJetFromPfClusters +
  process.MatchAK4GenJetWithAK4JetFromPfClusters +
  process.MatchAK4GenJetWithAK4JetFromPfCandidates
)

#process.p = cms.Path(process.FastjetJetProducer + process.L1TJetPhase1Producer + process.PrintMomentum)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring(
    "keep *",
  ),
)

process.e = cms.EndPath(process.out)
