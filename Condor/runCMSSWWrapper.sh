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

cmsrel CMSSW_9_0_0
#rm -r CMSSW_10_1_5/src
#mkdir CMSSW_tmp 

#/usr/bin/hdfs dfs -copyToLocal /user/sb17498/software/CMSSW_10_1_5.tar.gz CMSSW_tmp
#cd CMSSW_tmp
#tar xzvf CMSSW_10_1_5.tar.gz
#mv CMSSW_10_1_5/src ../CMSSW_10_1_5/src
#cd ..


cd CMSSW_9_0_0/src
cmsenv
git cms-merge-topic simonecid:convolutionCurves_MuonMatch_GenMuonPropagation
git cms-addpkg Condor
git cms-addpkg L1TJetConvolutionCurves
scram b

mkdir __output

cmsRun Condor/CMSSWJobDriver.py blockIndex=${processId} numberOfBlocks=${numberOfJobs} outputFile=__output/${OUTPUT_FILENAME} cfgFile=${inputFile} inputFiles=${SAMPLE_FILE}

echo "Will save on" /FCC-hh/${HDFS_DEST}

/usr/bin/hdfs dfs -mkdir -p /FCC-hh/${HDFS_DEST}
/usr/bin/hdfs dfs -moveFromLocal __output/* /FCC-hh/${HDFS_DEST}

set +o xtrace