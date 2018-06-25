import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from L1TJetConvolutionCurves.MatchL1TMuonWithGenLevelMuons.source_SingleMu import *

process = cms.Process("MatchL1TMuonWithGenLevelMuons")

process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(
    process.GlobalTag, '90X_upgrade2023_realistic_v9', '')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
options.register ('source',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Source sample")
options.outputFile = 'MatchL1TMuonWithGenLevelMuons_SingleMu.root'
options.source = "source"
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if options.source:
  process.source = locals()[options.source]

process.TFileService = cms.Service('TFileService', fileName = cms.string(options.outputFile))

process.MatchL1TMuonWithGenLevelMuons = cms.EDAnalyzer('MatchL1TMuonWithGenLevelMuons',
  genParticleCollectionTag = cms.InputTag( "PropagateGenMuonAndGenJetsTo2ndMuonStation", "PropagatedGenMuons", "MatchL1TMuonWithGenLevelMuons"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
)

process.PropagateGenMuonAndGenJetsTo2ndMuonStation = cms.EDProducer('PropagateGenMuonAndGenJetsTo2ndMuonStation',
  genParticleCollectionTag = cms.InputTag("genParticles"),
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  prop2nd = cms.PSet(
    useTrack = cms.string("none"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
    useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
    useSimpleGeometry = cms.bool(True),
    useStation2 = cms.bool(True),
    fallbackToME1 = cms.bool(False)
  )
)

process.p = cms.Path(process.PropagateGenMuonAndGenJetsTo2ndMuonStation *
                     process.MatchL1TMuonWithGenLevelMuons)
