import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from importlib import import_module

process = cms.Process("MatchL1TMuonWithGenLevelMuons")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
options.register ('source',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Source sample")
options.register ('sourceFile',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "File containing the splitted sources")
options.outputFile = 'MatchL1TMuonWithGenLevelMuons_GluGlu_HToMuMu.root'
options.source = "source_0"
options.sourceFile = "source_SingleNeutrinoPU140_splitted"

options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = getattr(import_module("L1TJetConvolutionCurves.MatchGenJetWithL1Objects." + options.sourceFile), options.source)

process.TFileService = cms.Service('TFileService', fileName = cms.string(options.outputFile))

process.MatchL1TMuonWithGenLevelMuons = cms.EDAnalyzer('MatchL1TMuonWithGenLevelMuons',
  genParticleCollectionTag = cms.InputTag("genParticles"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
)


process.p = cms.Path(process.MatchL1TMuonWithGenLevelMuons)
