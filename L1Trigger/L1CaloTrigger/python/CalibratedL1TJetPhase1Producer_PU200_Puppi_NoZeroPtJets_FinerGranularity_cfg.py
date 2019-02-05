from L1Trigger.L1CaloTrigger.CalibratedL1TJetPhase1Producer_cfg import *

CalibrateAK4JetFromPfClusters_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_FinerGranularity_MatchAK4GenJetWithAK4JetFromPfClusters.pickle"
CalibrateAK4JetFromPfCandidates_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_FinerGranularity_MatchAK4GenJetWithAK4JetFromPfCandidates.pickle"
CalibratePhase1L1TJetFromPfClusters_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_FinerGranularity_MatchAK4GenJetWithPhase1L1TJetFromPfClusters.pickle"
CalibratePhase1L1TJetFromPfCandidates_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Trees_Calibration_PU200_Puppi_NoZeroPtJets_FinerGranularity_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates.pickle"

with open(CalibratePhase1L1TJetFromPfClusters_path, 'rb') as f:
  process.CalibratePhase1L1TJetFromPfClusters.calibration = pkl.load(f)
with open(CalibratePhase1L1TJetFromPfCandidates_path, 'rb') as f:
  process.CalibratePhase1L1TJetFromPfCandidates.calibration = pkl.load(f)
with open(CalibrateAK4JetFromPfClusters_path, 'rb') as f:
  process.CalibrateAK4JetFromPfClusters.calibration = pkl.load(f)
with open(CalibrateAK4JetFromPfCandidates_path, 'rb') as f:
  process.CalibrateAK4JetFromPfCandidates.calibration = pkl.load(f)
