#! /bin/bash

set -o xtrace

while getopts "j:c:p:s:i:" o; do
  case "${o}" in
    j)
      jobName=${OPTARG}
      ;;
    c)
      clusterId=${OPTARG}
      ;;
    p)
      processId=${OPTARG}
      ;;  
    s)
      sourceFile=${OPTARG}
      ;;
    i)
      inputFile=${OPTARG}
      ;;
    esac
done

echo "Dumping sysinfo"

echo "> lsb_release -a"
lsb_release -a

echo "> Platform name"
python /cvmfs/fcc.cern.ch/sw/0.8.1/tools/hsf_get_platform.py --get=os 

echo "I am running on" $HOSTNAME

echo "Running heppy job"

source /cvmfs/cms.cern.ch/cmsset_default.sh

set -o xtrace
cp -r /software/sb17498/CMSSW_9_0_0 .
cd CMSSW_9_0_0/src
cmsenv
scramv1 b ProjectRename
scram b
cmsenv

mkdir ${SAVE_DESTINATION}

heppy ${HOME_FOLDER}/${SAVE_DESTINATION} ${inputFile} --option=sample=source_${processId} --option=sourceFile=${sourceFile} --option=outputFile=output.root -f

# Zip file
cd ${HOME_FOLDER}
tar -czvf ${jobName}_${sampleName}_${clusterId}.${processId}.tar.gz ${SAVE_DESTINATION}

set +o xtrace
