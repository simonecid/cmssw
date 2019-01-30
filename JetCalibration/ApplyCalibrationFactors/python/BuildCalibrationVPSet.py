import FWCore.ParameterSet.Config as cms
import os
import pickle as pkl
from sys import argv

def BuildCalibrationVPSet(rootfile_path, etaBinning):

  from ROOT import TFile
  parameters = cms.VPSet()
  
  tfile = TFile(rootfile_path)

  ptBinning = []

  for binIdx in xrange(0, len(etaBinning) -1 ):

    parameter = cms.PSet()
    
    graph = tfile.Get("l1corr_eta_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]))
    
    if graph == None:
      raise ValueError("Object " + "l1corr_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]) + " not found in file " + rootfile_path)
    
    l1tPtCentres = [pt for pt in graph.GetX()]

    parameter.l1tPtBins = cms.vdouble([float("-inf")] + [(l1tPtCentres[i] + l1tPtCentres[i+1])/2. for i in xrange(0, len(l1tPtCentres)-1)] + [float("+inf")])
    parameter.l1tCalibrationFactors = cms.vdouble([y for y in graph.GetY()])
    parameter.etaMin = cms.double(etaBinning[binIdx])
    parameter.etaMax = cms.double(etaBinning[binIdx + 1])
    
    parameters.append(parameter)

  tfile.Close()

  with open(os.path.splitext(rootfile_path)[0] + '.pickle', 'wb') as f:
    pkl.dump(parameters, f, pkl.HIGHEST_PROTOCOL)

  return


if __name__ == "__main__":
    BuildCalibrationVPSet("JetCalibration/Calibration_MatchAK4GenJetWithAK4JetFromPfCandidates_PU200_NoZeroPtJets_FinerGranularity.root", [0, 1.4, 3, 5.191])
    BuildCalibrationVPSet("JetCalibration/Calibration_MatchAK4GenJetWithAK4JetFromPfClusters_PU200_NoZeroPtJets_FinerGranularity.root", [0, 1.4, 3, 5.191])
    BuildCalibrationVPSet("JetCalibration/Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_PU200_NoZeroPtJets_FinerGranularity.root", [0, 1.4, 3, 5.191])
    BuildCalibrationVPSet("JetCalibration/Calibration_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_PU200_NoZeroPtJets_FinerGranularity.root", [0, 1.4, 3, 5.191])