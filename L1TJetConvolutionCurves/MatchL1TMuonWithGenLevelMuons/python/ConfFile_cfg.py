import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from L1TJetConvolutionCurves.MatchL1TMuonWithGenLevelMuons.source_SingleNeutrinoPU140_splitted import *

process = cms.Process("MatchL1TMuonWithGenLevelMuons")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
options.register ('source',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Source sample")
options.outputFile = 'MatchL1TMuonWithGenLevelMuons.root'
options.source = "source_0"
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

if options.source:
  process.source = locals()[options.source]

process.TFileService = cms.Service('TFileService', fileName = cms.string(options.outputFile))

process.MatchL1TMuonWithGenLevelMuons = cms.EDAnalyzer('MatchL1TMuonWithGenLevelMuons',
  genParticleCollectionTag = cms.InputTag("genParticles"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
)


process.p = cms.Path(process.MatchL1TMuonWithGenLevelMuons)