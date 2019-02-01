#!/usr/bin/env python

import os

TreeNames = [
             "MatchAK4GenJetWithAK4JetFromPfClusters",
             "MatchAK4GenJetWithAK4JetFromPfCandidates",
             "MatchAK4GenJetWithPhase1L1TJetFromPfClusters",
             "MatchAK4GenJetWithPhase1L1TJetFromPfCandidates",
]
InputFiles = [
            "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets/Trees_Calibration_PU200_Puppi_NoZeroPtJets.root",
            "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PU200_Puppi_NoZeroPTJets_FinerGranularity/Trees_Calibration_PU200_Puppi_NoZeroPtJets_FinerGranularity.root",
            "/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/ComputeUncalibratedPhase1AndAK4L1TJetsFromPfClustersAndCandidates_QCD_PUSubtractionBranch_PU200_Puppi_NoZeroPTJets/Trees_Calibration_PUSubtractionBranch_PU200_Puppi_NoZeroPtJets.root",
]

tmpFileName = "/tmp/sb17498_tmpjob.job"

for InputFile in InputFiles:
  for TreeName in TreeNames:
  
    JobName = "Calibration_" + os.path.splitext(os.path.basename(InputFile))[0] + "_" + TreeName

    tmpJob = open(tmpFileName, "w")

    jobDescription =  \
"""Universe = vanilla
job = {JobName}
cmd = JetCalibration/ApplyCalibrationFactors/python/runCalibration.py
args= --no-correction-fit --stage2 --treeName {TreeName} {InputFile} {OutputFile}
transfer_input_files = JetCalibration/ApplyCalibrationFactors/python/finer_binning.py,JetCalibration/ApplyCalibrationFactors/python/common_utils.py,JetCalibration/ApplyCalibrationFactors/python/__init__.py
# Better not to put logs in the same folder, as they are plenty
output = /six/sb17498/logs/CMSSW/$(job)_$(cluster).$(process).out
error = /six/sb17498/logs/CMSSW/$(job)_$(cluster).$(process).err
log = /six/sb17498/logs/CMSSW/$(job)_$(cluster).$(process).log
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
#+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/kreczko/workernode:centos7"
getenv=True

request_cpus = 1
request_memory = 2000
request_disk = 3000000

queue 1""".format(JobName=JobName, TreeName=TreeName, InputFile=InputFile, OutputFile=JobName+".root")

    print "Submitting", JobName + ":\n" + jobDescription + "\n-----------------------"
    
    tmpJob.write(jobDescription)
    tmpJob.close()
    os.system("condor_submit " + tmpFileName)


