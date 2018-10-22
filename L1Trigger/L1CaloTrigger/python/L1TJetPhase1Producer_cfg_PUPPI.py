from L1Trigger.L1CaloTrigger.L1TJetPhase1Producer_cfg import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

ak4JetFromPfClustersParameters.src = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN")
