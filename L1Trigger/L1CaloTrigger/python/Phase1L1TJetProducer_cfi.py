import FWCore.ParameterSet.Config as cms
from math import pi

caloEtaSegmentation = cms.vdouble(
  0.0,
  0.0833,
  0.1666,
  0.2499,
  0.3332,
  0.4165,
  0.4998,
  0.5831,
  0.6664,
  0.7497,
  0.833,
  0.9163,
  0.9996,
  1.0829,
  1.1662,
  1.2495,
  1.3328,
  1.4161,
  1.5,
)

Phase1L1TJetProducer = cms.EDProducer('Phase1L1TJetProducer',
  inputCollectionTag = cms.InputTag("l1pfCandidates", "Puppi", "REPR"),
  etaBinning = caloEtaSegmentation,
#  nBinsPhi = cms.uint32(72),
#  phiLow = cms.double(-pi),
#  phiUp = cms.double(pi),
  nBinsPhi = cms.uint32(8),
  phiLow = cms.double(0),
  phiUp = cms.double(0.7),
  jetIEtaSize = cms.uint32(5),
  jetIPhiSize = cms.uint32(5),
  seedPtThreshold = cms.double(5), # GeV
  puSubtraction = cms.bool(False),
  outputCollectionName = cms.string("UncalibratedPhase1L1TJetFromPfCandidates"),
  vetoZeroPt = cms.bool(True)
)
