import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from importlib import import_module
#from L1TJetConvolutionCurves.SaveGenLevelInfo.source_SingleNeutrinoPU140_splitted import *

process = cms.Process("SaveGenLevelInfo")
#process = cms.Process("SaveEvent")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

options = VarParsing.VarParsing('analysis')
options.register('source',
                 "",  # default value
                 VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                 VarParsing.VarParsing.varType.string,          # string, int, or float
                 "Source sample")
options.outputFile = 'genInfo.root'
options.source = "source_MinBias_TuneCUETP8M1_14TeV"
options.register('sourceFile',
                 "",  # default value
                 VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                 VarParsing.VarParsing.varType.string,          # string, int, or float
                 "File containing the splitted sources")
options.sourceFile = "source_MinBias_TuneCUETP8M1_14TeV"
options.parseArguments()

process.source = getattr(import_module(
    "L1TJetConvolutionCurves.MatchGenJetWithL1Objects." + options.sourceFile), options.source)

# UNCOMMENT ME!
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(5000))

process.TFileService = cms.Service(
    'TFileService', fileName=cms.string(options.outputFile))

process.out = cms.OutputModule("PoolOutputModule",
                               fileName=cms.untracked.string(
                                   'myOutputFile.root'),
                               SelectEvents=cms.untracked.PSet(
                                   SelectEvents=cms.vstring("p")
                               )
                               )


#tempFiles = cms.untracked.vstring()
#tempSource = cms.Source("PoolSource", fileNames=tempFiles)
#tempFiles.extend([
#    '/store/mc/PhaseIITDRFall17GS/MinBias_TuneCUETP8M1_14TeV-v1-pythia8/GEN-SIM/93X_upgrade2023_realistic_v2-v1/10000/000B06D1-769F-E711-AF73-02163E011F7D.root'
#])
#process.source = tempSource

process.SaveGenLevelInfo = cms.EDAnalyzer('SaveGenLevelInfo',
                                            genParticleCollectionTag=cms.InputTag(
                                              "genParticles"),
                                            genJetCollectionTag=cms.InputTag(
                                              "ak4GenJets"),
                                            genJetNoNuCollectionTag=cms.InputTag(
                                              "ak4GenJetsNoNu"),
                                          )


process.SaveEvent = cms.EDProducer('SaveEvent'
                                   )

process.EventNumberFilter = cms.EDFilter('EventNumberFilter'
                                         )

process.p = cms.Path(process.SaveGenLevelInfo)

#process.e = cms.EndPath(process.out)
