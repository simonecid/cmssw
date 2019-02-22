from imp import load_source
from copy import deepcopy
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.ParameterSet.Config as cms

def cantor(k1, k2):
  return ((k1 + k2) * (k1 + k2 + 1)) / 2 + k2

'''
#DEBUG
cfgFile = "FastPUPPI/NtupleProducer/python/runNewInputs.py"
sourceFilePath = "qcdfiles.py"
blockIndex = 1
numberOfBlocks = 2
'''

'''split the list aList in numberOfBlocks blocks of equal size'''
def splitInBlocks (aList, numberOfBlocks):
  k = len(aList) / numberOfBlocks
  r = len(aList) % numberOfBlocks

  i = 0
  blocks = []
  while i < len(aList):
    if len(blocks)<r:
      blocks.append(aList[i:i+k+1])
      i += k+1
    else:
      blocks.append(aList[i:i+k])
      i += k

  return blocks

options = VarParsing.VarParsing('analysis')
options.register('clusterId',
                 0,  # default value
                 VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                 VarParsing.VarParsing.varType.int,          # string, int, or float
                 "Cluster index")
options.register('processId',
                 0,  # default value
                 VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                 VarParsing.VarParsing.varType.int,          # string, int, or float
                 "Process index")
options.outputFile = 'output.root'
options.register('cfgFile',
                 "",  # default value
                 VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                 VarParsing.VarParsing.varType.string,          # string, int, or float
                 "Path to CMSSW cfg file to run")
options.parseArguments()

#process = load_source("process", "L1Trigger/L1CaloTrigger/python/test_L1TJetPhase1Producer_cfg.py")
# First we load the process from the file path specified as argument
'''Contains the process we want to execute'''
processModule = load_source(
                       "processModule",
                       options.cfgFile
                     )

process = processModule.process                     

# Loading the python module containing the source files that will be splitted
'''CMSSW source to be splitted'''
if len(options.inputFiles) > 0:
  inputFiles = open(options.inputFiles[0])
  process.source.fileNames = cms.untracked.vstring(inputFiles.readlines())
  inputFiles.close()
  # Splitting the source and reassigning the correct chunk to source
  process.source.fileNames = splitInBlocks(process.source.fileNames, options.numberOfBlocks)[options.blockIndex]


#Setting output file name 
if options.outputFile != "":
  process.out.fileName = options.outputFile
  if hasattr(process, "TFileService"):
    #if it has a fileservice we use the output file name without extension as a prefix
    process.TFileService.fileName = cms.string(options.outputFile.replace(".root", "") + "_" + process.TFileService.fileName.value())

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
k1 = options.clusterId
k2 = options.processId
seed = cantor(k1, k2) % (2**32)
seedList = [seed] * randSvc.countSeeds()
randSvc.insertSeeds(*seedList)