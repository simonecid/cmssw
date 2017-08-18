#! /bin/sh

#set -o xtrace

while getopts "j:c:p:i:d:" o; do
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
OUTPUT_FILENAME=${jobName}_${clusterId}.${processId}.root
source /cvmfs/cms.cern.ch/cmsset_default.sh

set -o xtrace
cp -r /software/sb17498/CMSSW_9_0_0 .
cd CMSSW_9_0_0/src
cmsenv
scramv1 b ProjectRename
scram b

cmsRun ${inputFile} source=source_${processId} outputFile=${OUTPUT_FILENAME}

echo "Will save on" /FCC-hh/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /FCC-hh/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal ${OUTPUT_FILENAME} /FCC-hh/${HDFS_DEST}

set +o xtrace