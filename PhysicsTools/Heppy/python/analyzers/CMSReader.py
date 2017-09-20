from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from importlib import import_module

import collections

'''Reads from the root files and creates iterable arrays.
'''

class CMSReader(Analyzer):
  '''
  
  Simple example configuration:: 
  
    from heppy.analyzers.CMSReader import CMSReader
    reader = cfg.Analyzer(
      CMSReader,
      reads = {
        "genParticles": {
          "handle": AutoHandle( 'genParticles', 'vector<reco::GenParticle>' ),
          "type": "GenParticle"
        },
        "l1tMuons": {
          "handle": AutoHandle( 'simGmtStage2Digis', 'BXVector<l1t::Muon>' ),
          "type": "TriggerObject"
        }
      }
    )
    
  * reads: dictionary of string -> handle, particle type

  For each key-handle pair, the reader creates a property "key" in the events and fills it with an array of object pointed by the handle.abs
  This workaround allows to easily uses every other heppy analyser.
  '''
  
  def declareHandles(self):
    super(CMSReader, self).declareHandles()
    for key, read in self.cfg_ana.reads.iteritems():
      self.handles[key] = read["handle"]
    
  def process(self, event):
    super(CMSReader, self).readCollections(event.input)
    for key, read in self.cfg_ana.reads.iteritems():
      module = import_module("PhysicsTools.Heppy.physicsobjects." + read["type"])
      object_model = getattr(module, read["type"])
      setattr(event,  key, map(object_model, read["handle"].product()))