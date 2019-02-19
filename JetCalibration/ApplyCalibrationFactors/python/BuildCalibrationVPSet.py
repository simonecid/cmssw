#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
import os
import pickle as pkl
from sys import argv

def BuildCalibrationVPSet(rootfile_path, etaBinning, outputfile_path):

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

  with open(outputfile_path, 'wb') as f:
    pkl.dump(parameters, f, pkl.HIGHEST_PROTOCOL)

  return


if __name__ == "__main__":
  BuildCalibrationVPSet("Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfCandidates_3752475.0.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964, 3.489, 4.191, 5.191],
  "Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfCandidates_3752475.0.pickle"
  )
  BuildCalibrationVPSet("Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfClusters_3752474.0.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964, 3.489, 4.191, 5.191],
  "Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithAK4JetFromPfClusters_3752474.0.pickle"
  )
  BuildCalibrationVPSet("Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_3752477.0.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964, 3.489, 4.191, 5.191],
  "Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_3752477.0.pickle"
  )
  BuildCalibrationVPSet("Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_3752476.0.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964, 3.489, 4.191, 5.191],
  "Calibration_Tree_Calibration_QCD_PU200_Puppi_NoZeroPTJets_0.4Square_MatchAK4GenJetWithPhase1L1TJetFromPfClusters_3752476.0.pickle"
  )