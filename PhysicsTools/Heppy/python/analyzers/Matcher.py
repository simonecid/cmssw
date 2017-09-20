from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from importlib import import_module

from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaR2

import collections

'''Particle matcher.
'''

def bestMatch( object, matchCollection):
  '''Return the best match to object in matchCollection, which is the closest object in delta R'''
  deltaR2Min = float('+inf')
  bm = None
  for match in matchCollection:
    dR2 = deltaR2( object, match )
    if dR2 < deltaR2Min:
      deltaR2Min = dR2
      bm = match
  return bm, deltaR2Min

def matchObjectCollection( objects, matchCollection, deltaR2Max, filter = lambda x,y : True):
  pairs = {}
  if len(objects)==0:
    return pairs
  x = 0
  for object in objects:
    if len(matchCollection)==0:
      pairs[x] = None
    else:
      bm, dr2 = bestMatch( object, [mob for mob in matchCollection if filter(object,mob)] )
      if dr2<deltaR2Max:
        import pdb; pdb.set_trace()        
        pairs[x] = bm
      else:
        pairs[x] = None            
    
    x += 1
        
  return pairs


class Matcher(Analyzer):
  '''Particle matcher. 
  
  Works with any kind of object with a p4 method. 
  
  Simple example configuration:: 
  
    from heppy.analyzers.Matcher import Matcher
    papas_jet_match = cfg.Analyzer(
      Matcher,
      instance_label = 'papas',   
      delta_r = 0.3,
      match_particles = "genJets",
      particles = "jets",
    )

  @param particles: handle containing the particles to be matched. 
  @param match_particles: handle of the collection containing the particles where a match 
         is to be found. 
  '''
  
  def process(self, event):
    '''process event
    
    The event must contain:
     - self.cfg_ana.particles: the particles to be matched
     - self.cfg_ana.match_particles: the particles in which the match has to be found
     
    Modifies the particles in event.<self.cfg_ana.particles>
    '''
    particles = getattr(event, self.cfg_ana.particles)
    match_particles = getattr(event, self.cfg_ana.match_particles)
    pairs = matchObjectCollection(particles, match_particles,
                      self.cfg_ana.delta_r**2)

    x = 0
    for ptc in particles:
      matchname = 'match'
      match = pairs[x]
      setattr(ptc, matchname, match)
      if match:
        drname = 'dr'
        dr = deltaR(ptc, match)
        import pdb; pdb.set_trace()
        setattr(ptc, drname, dr)
      x += 1
