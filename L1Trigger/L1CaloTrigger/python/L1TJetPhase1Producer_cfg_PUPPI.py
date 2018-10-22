from L1Trigger.L1CaloTrigger.L1TJetPhase1Producer_cfg import *

ak4JetFromPfClustersParameters.src = cms.InputTag("pfClustersFromCombinedCalo", "calibrated", "IN")
ak4JetFromPfCandidatesParameters.src = cms.InputTag("l1pfProducer", "Puppi", "IN")

