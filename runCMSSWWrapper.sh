#! /bin/sh

#set -o xtrace

export HOME=/users/sb17498

while getopts "j:c:p:i:n:d:s:b:" o; do
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
cmsrel CMSSW_10_1_5
rm -r CMSSW_10_1_5/src
cp -r /software/sb17498/CMSSW_10_1_5/src CMSSW_10_1_5
cd CMSSW_10_1_5/src
cmsenv
scramv1 b ProjectRename
#git cms-merge-topic simonecid:${BRANCH}
scram b

mkdir __output

cmsRun CMSSWJobDriver.py blockIndex=${processId} numberOfBlocks=${numberOfJobs} outputFile=__output/${OUTPUT_FILENAME} cfgFile=${inputFile} inputFiles=${SAMPLE_FILE}

echo "Will save on" /user/sb17498/CMS_Phase_2/jetMETStudies/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /user/sb17498/CMS_Phase_2/jetMETStudies/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal __output/* /user/sb17498/CMS_Phase_2/jetMETStudies/${HDFS_DEST}

set +o xtrace