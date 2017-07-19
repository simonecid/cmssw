import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from L1TJetConvolutionCurves.MatchGenJetWithL1Objects.source_QCD_Pt_15_3000_splitted import *

process = cms.Process("MatchGenJetWithL1Objects")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

options = VarParsing.VarParsing ('analysis')
options.register ('source',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Source sample")
options.outputFile = 'l1tObjectGenJetMatching.root'
options.source = source_0
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if options.source:
  process.source = locals()[options.source]

process.TFileService = cms.Service('TFileService', fileName = cms.string(options.outputFile))

process.MatchGenJetWithL1Objects = cms.EDAnalyzer('MatchGenJetWithL1Objects',
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  l1tJetCollectionTag = cms.InputTag("simCaloStage2Digis"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
  l1tTauCollectionTag = cms.InputTag("simCaloStage2Digis"),
  l1tEGammaCollectionTag = cms.InputTag("simCaloStage2Digis"),
)

process.p = cms.Path(process.MatchGenJetWithL1Objects)  