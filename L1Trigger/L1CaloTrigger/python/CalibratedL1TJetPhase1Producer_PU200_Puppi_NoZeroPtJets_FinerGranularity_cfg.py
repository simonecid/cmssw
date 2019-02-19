from L1Trigger.L1CaloTrigger.CalibratedL1TJetPhase1Producer_cfg import *

CalibrateAK4JetFromPfClusters_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfClusters_3752474.0.pickle"
CalibrateAK4JetFromPfCandidates_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfCandidates_3752475.0.pickle"
CalibratePhase1L1TJetFromPfClusters_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_3752476.0.pickle"
CalibratePhase1L1TJetFromPfCandidates_path = "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_3752477.0.pickle"


with open(CalibratePhase1L1TJetFromPfClusters_path, 'rb') as f:
  process.CalibratePhase1L1TJetFromPfClusters.calibration = pkl.load(f)
with open(CalibratePhase1L1TJetFromPfCandidates_path, 'rb') as f:
  process.CalibratePhase1L1TJetFromPfCandidates.calibration = pkl.load(f)
with open(CalibrateAK4JetFromPfClusters_path, 'rb') as f:
  process.CalibrateAK4JetFromPfClusters.calibration = pkl.load(f)
with open(CalibrateAK4JetFromPfCandidates_path, 'rb') as f:
  process.CalibrateAK4JetFromPfCandidates.calibration = pkl.load(f)
