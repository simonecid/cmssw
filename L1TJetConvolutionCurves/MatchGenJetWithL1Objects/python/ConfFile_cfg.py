import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from importlib import import_module
#from L1TJetConvolutionCurves.MatchGenJetWithL1Objects.source_SingleNeutrinoPU140_splitted import *

process = cms.Process("MatchGenJetWithL1Objects")

#process.Tracer = cms.Service("Tracer")
#process = cms.Process("SaveEvent")
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi')
#process.load('RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
options.register ('source',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Source sample")
options.outputFile = 'l1tObjectGenJetMatching.root'
options.source = "source_0"
options.register ('sourceFile',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "File containing the splitted sources")
options.sourceFile = "source_QCD_Pt_15_3000_splitted"
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = getattr(import_module("L1TJetConvolutionCurves.MatchGenJetWithL1Objects." + options.sourceFile), options.source)
#process.source.skipEvents = cms.untracked.uint32(1)

process.TFileService = cms.Service('TFileService', fileName = cms.string(options.outputFile))

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring("p")
    )
)

process.GenMuonCollectionProducer = cms.EDProducer('GenMuonCollectionProducer',
  genParticleCollectionTag = cms.InputTag("genParticles"),
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

process.MatchGenJetWithL1Objects = cms.EDAnalyzer('MatchGenJetWithL1Objects',
  genParticleCollectionTag = cms.InputTag( "GenMuonCollectionProducer", "GenMuons", "MatchGenJetWithL1Objects"),
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
)

process.MatchLeadingGenJetWithL1Objects = cms.EDAnalyzer('MatchLeadingGenJetWithL1Objects',
  genParticleCollectionTag = cms.InputTag("genParticles"),
  genJetCollectionTag = cms.InputTag("ak4GenJets"),
  l1tMuonCollectionTag = cms.InputTag("simGmtStage2Digis"),
)


process.SaveEvent = cms.EDProducer('SaveEvent'
)

process.EventNumberFilter = cms.EDFilter('EventNumberFilter'
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(
    process.GlobalTag, '90X_upgrade2023_realistic_v9', '')

#process.p = cms.Path(process.PropagateGenMuonAndGenJetsTo2ndMuonStation *
#                     process.MatchGenJetWithL1Objects)
process.p = cms.Path(process.GenMuonCollectionProducer *
                     process.MatchGenJetWithL1Objects)

#process.e = cms.EndPath(process.out)
