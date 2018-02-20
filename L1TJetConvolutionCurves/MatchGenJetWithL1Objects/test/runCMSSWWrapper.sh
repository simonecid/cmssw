#! /bin/sh

#set -o xtrace

export HOME=/users/sb17498

while getopts "j:c:p:i:d:s:b:" o; do
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
    s)
      SAMPLE_FILE=${OPTARG}
      ;;
    b)
      BRANCH=${OPTARG}
    esac
done

export SCRAM_ARCH=slc6_amd64_gcc630

echo "I am running on" $HOSTNAME

echo "Running CMSSW job"

HOME_FOLDER="$(pwd)"
OUTPUT_FILENAME=${jobName}_${clusterId}.${processId}.root
source /cvmfs/cms.cern.ch/cmsset_default.sh

set -o xtrace
cmsrel CMSSW_9_3_0
cd CMSSW_9_3_0/src
cmsenv
git cms-merge-topic simonecid:${BRANCH}
scram b

cmsRun ${inputFile} source=source_${processId} outputFile=${OUTPUT_FILENAME} sourceFile=${SAMPLE_FILE}

echo "Will save on" /FCC-hh/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /FCC-hh/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal ${OUTPUT_FILENAME} /FCC-hh/${HDFS_DEST}

set +o xtrace