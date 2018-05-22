import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from importlib import import_module
#from L1TJetConvolutionCurves.MatchGenJetWithL1Objects.source_SingleNeutrinoPU140_splitted import *

#BXVector<l1t::EtSum> "simCaloStage2Digis"        ""                "HLT"
#BXVector<l1t::EtSum> "hltGtStage2Digis"          "EtSum"           "HLT"
#BXVector<l1t::EtSum> "simCaloStage2Digis"        "MP"              "HLT"


process = cms.Process("MatchGenJetWithL1Objects")
#process = cms.Process("SaveEvent")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
options.register ('source',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Source sample")
options.outputFile = 'MHTTree.root'
options.source = "source_0"
options.register ('sourceFile',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "File containing the splitted sources")
options.sourceFile = "source_SingleNeutrinoPU140_splitted"
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = getattr(import_module("L1TJetConvolutionCurves.MatchGenJetWithL1Objects." + options.sourceFile), options.source)
#process.source.skipEvents = cms.untracked.uint32(1)

process.TFileService = cms.Service('TFileService', fileName = cms.string(options.outputFile))

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring("p")
    )
)

process.ComputeMHT = cms.EDAnalyzer('ComputeMHT',
  l1tEtSumCollectionTag = cms.InputTag("simCaloStage2Digis", "", "HLT"),
)

process.p = cms.Path(process.ComputeMHT)
