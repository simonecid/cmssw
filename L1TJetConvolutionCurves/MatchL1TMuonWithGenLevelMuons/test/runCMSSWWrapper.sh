#! /bin/sh

#set -o xtrace

while getopts "j:c:p:i:d:b:" o; do
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
    b)
      BRANCH=${OPTARG}
    esac
done

echo "I am running on" $HOSTNAME

echo "Running CMSSW job"

HOME_FOLDER="$(pwd)"
OUTPUT_FILENAME=${jobName}_${clusterId}.${processId}.root
source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc6_amd64_gcc630

set -o xtrace
cmsrel CMSSW_9_0_0
cd CMSSW_9_0_0/src
set +o xtrace
cmsenv
set -o xtrace
git cms-merge-topic simonecid:${BRANCH}
scram b

cmsRun ${inputFile} source=source_${processId} outputFile=${OUTPUT_FILENAME}

echo "Will save on" /FCC-hh/${HDFS_DEST}

mv ${OUTPUT_FILENAME} ${HOME_FOLDER}

cd ${HOME_FOLDER}

/usr/bin/hdfs dfs -mkdir -p /FCC-hh/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal ${OUTPUT_FILENAME} /FCC-hh/${HDFS_DEST}

set +o xtrace