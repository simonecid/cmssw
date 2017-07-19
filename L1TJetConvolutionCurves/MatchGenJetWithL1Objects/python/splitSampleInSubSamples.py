from L1TJetConvolutionCurves.MatchGenJetWithL1Objects.MySamples_cfg import *

def splitInBlocks (l, n):
  """split the list l in n blocks of equal size"""
  k = len(l) / n
  r = len(l) % n

  i = 0
  blocks = []
  while i < len(l):
    if len(blocks)<r:
      blocks.append(l[i:i+k+1])
      i += k+1
    else:
      blocks.append(l[i:i+k])
      i += k

  return blocks

import argparse
parser = argparse.ArgumentParser()

from L1TJetConvolutionCurves.MatchGenJetWithL1Objects.MySamples_cfg import *

parser.add_argument('--numberOfBlocks',
    type=int, dest="numberOfBlocks")
parser.add_argument('--sample',
    dest="sample")

args = parser.parse_args()

sample = locals()[args.sample]

blocks = splitInBlocks(sample.fileNames, args.numberOfBlocks)

filelist = open(args.sample + "_splitted.py", 'w')
filelist.write("import FWCore.ParameterSet.Config as cms \n")
for idx, block in enumerate(blocks):
  filelist.write("files_" + str(idx) + " = cms.untracked.vstring() \n")
  filelist.write("files_" + str(idx) + ".extend(" + str(block) + ") \n")
  filelist.write("source_" + str(idx) + " = cms.Source (\"PoolSource\", fileNames = files_" + str(idx) + ") \n")
filelist.close()

