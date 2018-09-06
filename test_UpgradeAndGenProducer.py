
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1NTUPLE")

# conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond_condDBv2 import autoCond
process.GlobalTag.globaltag = cms.string( autoCond['run2_mc'] )

# input file
from L1Trigger.L1TNtuples.RelValInputFiles import *
process.source = cms.Source (
  "PoolSource",
  fileNames = cms.untracked.vstring( "file:/six/sb17498/CMSSW_10_1_5/src/myOutputFile.root" )
)

# N events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Tree.root')
)

# producer under test
process.load("L1Trigger.L1TNtuples.l1UpgradeTree_cfi")
process.load("L1Trigger.L1TNtuples.l1GeneratorTree_cfi")

l1GeneratorTree = cms.EDAnalyzer(
    "L1GenTreeProducer",
    genJetToken     = cms.untracked.InputTag("ak4GenJetsNoNu"),
    genParticleToken = cms.untracked.InputTag("genParticles"),
    pileupInfoToken     = cms.untracked.InputTag("addPileupInfo")
)


process.l1UpgradeTree = cms.EDAnalyzer(
    "L1UpgradeTreeProducer",
    egToken = cms.untracked.InputTag("caloStage2Digis","EGamma"),
    tauTokens = cms.untracked.VInputTag(cms.InputTag("caloStage2Digis","Tau")),
    muonToken = cms.untracked.InputTag("gmtStage2Digis","Muon"),
    muonLegacyToken = cms.untracked.InputTag("muonLegacyInStage2FormatDigis","legacyMuon"),
    sumToken = cms.untracked.InputTag("caloStage2Digis","EtSum"),
    jetToken = cms.untracked.InputTag("L1TJetPhase1Producer", "Phase1L1TJetFromPfCandidates", "L1TJetPhase1Producer"),
    maxL1Upgrade = cms.uint32(60)
)

process.p = cms.Path(
  process.l1GeneratorTree +
  process.l1UpgradeTree
)

