from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from importlib import import_module

'''
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle

Prints a single quantity.

Example:

def pt (ptc):
  return ptc.pt()

examplePrinter = cfg.Analyzer(
  Printer,
  'examplePrinter',
  input_objects = "jets"
  value_func = pt
)

 * input_objects: name of object collection
 * value_func: function that returns the variable to be returned. A particle object is passed to the function.

'''


class Printer(Analyzer):

  def process(self, event):
    self.input_objects = getattr(event, self.cfg_ana.input_objects)
    if len(self.input_objects) > 0: 
      print "Printing for", event.iEv 
    for x in xrange(0,len(self.input_objects)):
      ptc = self.input_objects[x]
      print x, self.cfg_ana.value_func(ptc)