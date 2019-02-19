import FWCore.ParameterSet.Config as cms
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
    "file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_0.4Square/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_0.8Square_3751805.0.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.32.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.46.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.475.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.457.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.382.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.374.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.184.root",
    #"file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_FinerGranularity_3742478.63.root",
  )
)

process.TFileService = cms.Service('TFileService', fileName = cms.string("Histograms.root"))
#process.source.skipEvents = cms.untracked.uint32(1)

with open('/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PUSubtractionBranch_PU200_Puppi_NoZeroPtJets_MatchAK4GenJetWithAK4JetFromPfClusters.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithAK4JetFromPfClusters = pkl.load(f)
with open('/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_MatchAK4GenJetWithAK4JetFromPfCandidates.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates = pkl.load(f)
with open('/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_MatchAK4GenJetWithPhase1L1TJetFromPfClusters.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters = pkl.load(f)
with open('/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates.pickle', 'rb') as f:
  Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates = pkl.load(f)

process.CalibratePhase1L1TJetFromPfClusters = cms.EDProducer('ApplyCalibrationFactors',
  inputCollectionTag = cms.InputTag("Phase1L1TJetFromPfClustersProducer", "Phase1L1TJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
  absEtaBinning = cms.vdouble([p.etaMin.value() for p in Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters] + [Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters[-1].etaMax.value()]),
  calibration = Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters,
  outputCollectionName = cms.string("CalibratedPhase1L1TJetFromPfClusters")
)

process.CalibratePhase1L1TJetFromPfCandidates = cms.EDProducer('ApplyCalibrationFactors',
   inputCollectionTag = cms.InputTag("Phase1L1TJetFromPfCandidatesProducer", "Phase1L1TJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
   absEtaBinning = cms.vdouble([p.etaMin.value() for p in Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates] + [Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates[-1].etaMax.value()]),
   calibration = Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates,
   outputCollectionName = cms.string("CalibratedPhase1L1TJetFromPfCandidates")
)

process.CalibrateAK4JetFromPfClusters = cms.EDProducer('ApplyCalibrationFactors',
   inputCollectionTag = cms.InputTag("ConvertGenJetToCaloJet", "ak4CaloJetFromPfClusters", "UncalibratedL1TJetPhase1Producer"),
   absEtaBinning = cms.vdouble([p.etaMin.value() for p in Calibration_MatchAK4GenJetWithAK4JetFromPfClusters] + [Calibration_MatchAK4GenJetWithAK4JetFromPfClusters[-1].etaMax.value()]),
   calibration = Calibration_MatchAK4GenJetWithAK4JetFromPfClusters,
   outputCollectionName = cms.string("CalibratedAK4JetFromPfClusters")
)

process.CalibrateAK4JetFromPfCandidates = cms.EDProducer('ApplyCalibrationFactors',
   inputCollectionTag = cms.InputTag("ConvertGenJetToCaloJet", "ak4CaloJetFromPfCandidates", "UncalibratedL1TJetPhase1Producer"),
   absEtaBinning = cms.vdouble([p.etaMin.value() for p in Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates] + [Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates[-1].etaMax.value()]),
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

CalibrateAK4JetFromPfClusters_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfClusters_3752474.0.pickle"
CalibrateAK4JetFromPfCandidates_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfCandidates_3752475.0.pickle"
CalibratePhase1L1TJetFromPfClusters_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_3752476.0.pickle"
CalibratePhase1L1TJetFromPfCandidates_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_3752477.0.pickle"

with open(CalibratePhase1L1TJetFromPfClusters_path, 'rb') as f:
  process.CalibratePhase1L1TJetFromPfClusters.calibration = pkl.load(f)
with open(CalibratePhase1L1TJetFromPfCandidates_path, 'rb') as f:
  process.CalibratePhase1L1TJetFromPfCandidates.calibration = pkl.load(f)
with open(CalibrateAK4JetFromPfClusters_path, 'rb') as f:
  process.CalibrateAK4JetFromPfClusters.calibration = pkl.load(f)
with open(CalibrateAK4JetFromPfCandidates_path, 'rb') as f:
  process.CalibrateAK4JetFromPfCandidates.calibration = pkl.load(f)

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
