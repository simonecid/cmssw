import os
import copy
import PhysicsTools.HeppyCore.framework.config as cfg
import logging
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
from PhysicsTools.HeppyCore.framework.looper import Looper
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.analyzers.examples.Printer import Printer
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
from PhysicsTools.HeppyCore.framework.heppy_loop import _heppyGlobalOptions
from PhysicsTools.Heppy.analyzers.Matcher import Matcher
from PhysicsTools.Heppy.analyzers.Selector import Selector
from PhysicsTools.Heppy.analyzers.CMSReader import CMSReader
from PhysicsTools.Heppy.analyzers.MatchedParticlesTreeProducer import MatchedParticlesTreeProducer

from importlib import import_module


sourceFileName = _heppyGlobalOptions["sourceFile"]
sampleName = _heppyGlobalOptions["sample"]
outputFile = _heppyGlobalOptions["outputFile"]

sample = getattr(import_module("L1TJetConvolutionCurves.MatchGenJetWithL1Objects." + sourceFileName), sampleName)

# next 2 lines necessary to deal with reimports from ipython
logging.shutdown()
reload(logging)
logging.basicConfig(level=logging.WARNING)

selectedComponents = [cfg.MCComponent(
  sampleName,
  files = [
#    "root://xrootd-cms.infn.it/" + path for path in sample.fileNames
#    "/hdfs/dpm/phy.bris.ac.uk/home/cms/store/mc/PhaseIISpring17D/SingleNeutrino/GEN-SIM-DIGI-RAW/PU140_90X_upgrade2023_realistic_v9-v1/70000/061FB060-8D26-E711-9DF5-0242AC130002.root"
#    "root://xrootd-cms.infn.it//store/mc/PhaseIISpring17D/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/100000/0013DDCD-8F46-E711-B875-0CC47A537688.root"
    "/storage/sb17498/CMSSW/CMSSW_9_0_0/test/0013DDCD-8F46-E711-B875-0CC47A537688.root"
  ]
)]

pdgIds = {
  'electron-': 11,
  'muon-': 13,
  'tau-': 15,
  'photon': 22,
  'pion+': 211,
  'kaon+': 321,
  'kaon_long': 130,
  'kaon_short': 310,
  'bottom': 5
}

tfile_service_1 = cfg.Service(
  TFileService,
  'tfile1',
  fname=outputFile,
  option='recreate'
)


def pt (ptc):
  return ptc.pt()

def eta (ptc):
  return ptc.eta()

def phi (ptc):
  return ptc.phi()

def matchedParticlePt (ptc):
  return ptc.match.pt()

def matchedParticleEta (ptc):
  return ptc.match.eta()

def ptRatioWithMatched (ptc):
  return ptc.pt()/ptc.match.pt()

def deltaR (ptc):
  return ptc.dr

def isMatched(ptc):
  return ptc.match is not None

def hasMatches(ptc):
  return len(ptc.matches) > 0

def particleCheckerFactory (ptcName):
  def particleChecker (ptc):
    return (abs(ptc.pdgId()) == pdgIds[ptcName])
  return particleChecker

def getFinalStateBQuark (ptc):
  return (abs(ptc.pdgid()) == 5) and (ptc.status() == 71)

def genParticleStatus (ptc):
  return ptc.status()

def pt (ptc):
  return ptc.pt()

cmsReader = cfg.Analyzer(
  CMSReader,
  reads = {
    "genParticles": {
      "handle": AutoHandle( 'genParticles', 'vector<reco::GenParticle>' ),
      "type": "GenParticle"
    },
    "l1tMuons": {
      "handle": AutoHandle( 'simGmtStage2Digis', 'BXVector<l1t::Muon>' ),
      "type": "PhysicsObject"
    }
  }
)

examplePrinter = cfg.Analyzer(
  Printer,
  'examplePrinter',
  input_objects = "genMuons",
  value_func = genParticleStatus
)

muonGenMuonMatcher = cfg.Analyzer(
  Matcher,
  instance_label = 'muonGenMuonMatcher',
  delta_r = 0.5,
  match_particles = "genMuons",
  particles = "l1tMuons"
)

genMuonSelector = cfg.Analyzer(
  Selector,
  'genMuonSelector',
  output = 'genMuons',
  input_objects = 'genParticles',
  filter_func = particleCheckerFactory("muon-")
)

matchedMuonSelector = cfg.Analyzer(
  Selector,
  instance_label = 'matchedMuonSelector',
  input_objects = 'l1tMuons',
  output = 'matched_l1tMuons',
  filter_func = isMatched
)

matchedL1TMuonGenMuonTreeProducer = cfg.Analyzer(
  MatchedParticlesTreeProducer,
  file_label = "tfile1",
  tree_name = 'matchedL1TMuonGenMuon',
  tree_title = 'Tree containing info about matched L1T muon and gen muons',
  particle_collection = 'matched_l1tMuons',
  particle_name = "l1tMuon",
  matched_particle_name = "genMuon"
)

# definition of a sequence of analyzers,
# the analyzers will process each event in this order
sequence = cfg.Sequence( [
  cmsReader,
  genMuonSelector,
  muonGenMuonMatcher,
  matchedMuonSelector,
  matchedL1TMuonGenMuonTreeProducer
] )

config = cfg.Config(
  components = selectedComponents,
  sequence = sequence,
  services = [tfile_service_1],
  events_class = Events
)

if __name__ == '__main__':

  def next():
      loop.process(loop.iEvent+1)

  loop = Looper( 'looper', config,
                 nEvents=100,
                 nPrint=0,
                 timeReport=True)
  loop.process()
  print loop.event
