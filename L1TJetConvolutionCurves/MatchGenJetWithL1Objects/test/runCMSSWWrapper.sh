#! /bin/bash

#set -o xtrace

while getopts "j:c:p:s:i:t:d:" o; do
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
      sampleName=${OPTARG}
      ;;
    i)
      inputFile=${OPTARG}
      ;;
    d)
      HDFS_DEST=${OPTARG}
      ;;
    esac
done

echo "I am running on" $HOSTNAME

echo "Running CMSSW job"

HOME_FOLDER="$(pwd)"
OUTPUT_FILENAME=${jobName}_${sampleName}_${clusterId}.${processId}.root

source /cvmfs/cms.cern.ch/cmsset_default.sh
cp -r /software/sb17498/CMSSW/CMSSW_9_0_0 .
cd CMSSW_9_0_0/src
cmsenv
scram b

set -o xtrace
cmsRun ${inputFile} source=${sampleName}_${processId} outputFile=${OUTPUT_FILENAME}

echo "Will save on" /FCC-hh/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /FCC-hh/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal ${OUTPUT_FILENAME} /FCC-hh/${HDFS_DEST}

set +o xtrace