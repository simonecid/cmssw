import FWCore.ParameterSet.Config as cms

process = cms.Process("PtPlotter")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
  # replace 'myfile.root' with the source file you want to use
  fileNames = cms.untracked.vstring(
    '/store/mc/PhaseIISpring17GS/MinBias_TuneCUETP8M1_14TeV-pythia8/GEN-SIM/90X_upgrade2023_realistic_v9-v1/00000/002B5E76-6025-E711-81F3-5065F381A2F1.root'
  )
)

process.demo = cms.EDAnalyzer('PtPlotter',
  GenJetTag=cms.InputTag("ak4GenJets")
)


process.p = cms.Path(process.demo)
