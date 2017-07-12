import FWCore.ParameterSet.Config as cms

process = cms.Process("MatchGenJetL1Object")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/s/sbologna/CMSSW/CMSSW_9_0_0/src/step2.root'
    )
)

process.demo = cms.EDAnalyzer('MatchGenJetWithL1Object'
)


process.p = cms.Path(process.demo)
