#! /bin/sh

#set -o xtrace

export HOME=/users/sb17498

while getopts "j:c:p:in:d:s:b:" o; do
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
    n) 
      numberOfJobs=${OPTARG}
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

cmsrel CMSSW_9_0_0
cd CMSSW_9_0_0/src
cmsenv 
git cms-merge-topic simonecid:${BRANCH}
scram b

mkdir __output

cmsRun CMSSWJobDriver.py numberOfBlocks=${numberOfJobs} outputFile=__output/${OUTPUT_FILENAME} cfgFile=${inputFile}

echo "Will save on" /user/sb17498/CMS_Phase_2/jetMETStudies/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /user/sb17498/CMS_Phase_2/jetMETStudies/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal __output/* /user/sb17498/CMS_Phase_2/jetMETStudies/${HDFS_DEST}

set +o xtrace
