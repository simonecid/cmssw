#! /bin/sh

#set -o xtrace

export HOME=/users/sb17498

while getopts "j:c:p:i:d:b:" o; do
  case "${o}" in
    j)
      jobName=${OPTARG}
      echo jobName=${OPTARG}
      ;;
    c)
      clusterId=${OPTARG}
      echo clusterId=${OPTARG}
      ;;
    p)
      processId=${OPTARG}
      echo processId=${OPTARG}
      ;;
    i)
      inputFile=${OPTARG}
      echo inputFile=${OPTARG}
      ;;
    d)
      HDFS_DEST=${OPTARG}
      echo HDFS_DEST=${OPTARG}
      ;;
    b)
      BRANCH=${OPTARG}
      echo BRANCH=${OPTARG}
    esac
done

export SCRAM_ARCH=slc6_amd64_gcc630

echo "I am running on" $HOSTNAME

echo "Running CMSSW job"

HOME_FOLDER="$(pwd)"
OUTPUT_FILENAME=${jobName}_${clusterId}.${processId}.root
source /cvmfs/cms.cern.ch/cmsset_default.sh

set -o xtrace

cmsrel CMSSW_9_0_0
cd CMSSW_9_0_0/src
cmsenv 
git cms-merge-topic simonecid:${BRANCH}
git cms-addpkg MuonGenerator

scram b

mkdir __output

cmsRun MuonGenerator/CMSSWJobDriver.py clusterId=${clusterId} processId=${processId} outputFile=__output/${OUTPUT_FILENAME} cfgFile=${inputFile}

echo "Will save on" /FCC-hh/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /FCC-hh/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal __output/* /FCC-hh/${HDFS_DEST}

set +o xtrace
