import FWCore.ParameterSet.Config as cms
from L1Trigger.L1CaloTrigger.caloEtaSegmentation import caloEtaSegmentation
from math import pi
from copy import deepcopy

import pickle as pkl

import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process = cms.Process("CalibratedL1TJetPhase1Producer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.337.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.32.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.46.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.475.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.457.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.382.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.374.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.184.root",
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_3738062.63.root",
  )
)

process.TFileService = cms.Service('TFileService', fileName = cms.string("Histograms.root"))
#process.source.skipEvents = cms.untracked.uint32(1)

with open('JetCalibration/Calibration_MatchAK4GenJetWithAK4JetFromPfClusters_PU200_PUSubtraction_NoZeroPtJets.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithAK4JetFromPfClusters = pkl.load(f)
with open('JetCalibration/Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates_PU200_PUSubtraction_NoZeroPtJets.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates = pkl.load(f)
with open('JetCalibration/Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_PU200_PUSubtraction_NoZeroPtJets.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters = pkl.load(f)
with open('JetCalibration/Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_PU200_PUSubtraction_NoZeroPtJets.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates = pkl.load(f)

process.CalibratePhase1L1TJetFromPfClusters = cms.EDProducer('ApplyCalibrationFactors',
  inputCollectionTag = cms.InputTag("Phase1L1TJetFromPfClustersProducer", "Phase1L1TJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
  absEtaBinning = cms.vdouble(0, 1.4, 3, 6),
  calibration = Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters,
  outputCollectionName = cms.string("CalibratedPhase1L1TJetFromPfClusters")
)

process.CalibratePhase1L1TJetFromPfCandidates = cms.EDProducer('ApplyCalibrationFactors',
   inputCollectionTag = cms.InputTag("Phase1L1TJetFromPfCandidatesProducer", "Phase1L1TJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
   absEtaBinning = cms.vdouble(0, 1.4, 3, 6),
   calibration = Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates,
   outputCollectionName = cms.string("CalibratedPhase1L1TJetFromPfCandidates")
)

process.CalibrateAK4JetFromPfClusters = cms.EDProducer('ApplyCalibrationFactors',
   inputCollectionTag = cms.InputTag("ConvertGenJetToCaloJet", "ak4CaloJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
   absEtaBinning = cms.vdouble(0, 1.4, 3, 6),
   calibration = Calibration_MatchAK4GenJetWithAK4JetFromPfClusters,
   outputCollectionName = cms.string("CalibratedAK4JetFromPfClusters")
)

process.CalibrateAK4JetFromPfCandidates = cms.EDProducer('ApplyCalibrationFactors',
   inputCollectionTag = cms.InputTag("ConvertGenJetToCaloJet", "ak4CaloJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
   absEtaBinning = cms.vdouble(0, 1.4, 3, 6),
   calibration = Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates,
   outputCollectionName = cms.string("CalibratedAK4JetFromPfCandidates")
)

process.MatchAK4GenJetWithPhase1L1TJetFromPfClusters = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("CalibratePhase1L1TJetFromPfClusters", "CalibratedPhase1L1TJetFromPfClusters", "CalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithPhase1L1TJetFromPfCandidates = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("CalibratePhase1L1TJetFromPfCandidates", "CalibratedPhase1L1TJetFromPfCandidates", "CalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithAK4JetFromPfClusters = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("CalibrateAK4JetFromPfClusters", "CalibratedAK4JetFromPfClusters", "CalibratedL1TJetPhase1Producer"),
)

process.MatchAK4GenJetWithAK4JetFromPfCandidates = cms.EDAnalyzer('MatchGenJetToRecoCaloJet',
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
  caloJetCollectionTag = cms.InputTag("CalibrateAK4JetFromPfCandidates", "CalibratedAK4JetFromPfCandidates", "CalibratedL1TJetPhase1Producer"),
)

#process.p = cms.Path(process.ak4JetFromPfClusters + process.ak4JetFromPfCandidates + process.UncalibratedL1TJetPhase1Producer + process.MatchGenJetToRecoCaloJet)

process.p = cms.Path(
  process.CalibratePhase1L1TJetFromPfClusters +
  process.CalibratePhase1L1TJetFromPfCandidates +
  process.CalibrateAK4JetFromPfClusters +
  process.CalibrateAK4JetFromPfCandidates +
  process.MatchAK4GenJetWithPhase1L1TJetFromPfClusters +
  process.MatchAK4GenJetWithPhase1L1TJetFromPfCandidates +
  process.MatchAK4GenJetWithAK4JetFromPfClusters +
  process.MatchAK4GenJetWithAK4JetFromPfCandidates
)

#process.p = cms.Path(process.FastjetJetProducer + process.UncalibratedL1TJetPhase1Producer + process.PrintMomentum)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring(
    "keep *",
  ),
)

process.e = cms.EndPath(process.out)
