from L1Trigger.L1CaloTrigger.UncalibratedL1TJetPhase1Producer_FinerGranularity_cfg import *

ak4JetFromPfClustersParameters.src = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN")
ak4JetFromPfCandidatesParameters.src = cms.InputTag("l1pfProducer", "Puppi", "IN")

process.Phase1L1TJetFromPfCandidatesProducer.seedPtThreshold = cms.double(10)
process.Phase1L1TJetFromPfClustersProducer.seedPtThreshold = cms.double(10)
