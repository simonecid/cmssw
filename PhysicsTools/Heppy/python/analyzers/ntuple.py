#!/bin/env python

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )

# simple p4

def bookP4( tree, pName ):
    var(tree, '{pName}_e'.format(pName=pName))
    var(tree, '{pName}_pt'.format(pName=pName))
    var(tree, '{pName}_theta'.format(pName=pName))
    var(tree, '{pName}_eta'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))
    var(tree, '{pName}_m'.format(pName=pName))

def fillP4( tree, pName, p4 ):
    fill(tree, '{pName}_e'.format(pName=pName), p4.e() )
    fill(tree, '{pName}_pt'.format(pName=pName), p4.pt() )
    fill(tree, '{pName}_theta'.format(pName=pName), p4.theta() )
    fill(tree, '{pName}_eta'.format(pName=pName), p4.eta() )
    fill(tree, '{pName}_phi'.format(pName=pName), p4.phi() )
    fill(tree, '{pName}_m'.format(pName=pName), p4.M() )

# simple particle

def bookParticle( tree, pName ):
    var(tree, '{pName}_pdgId'.format(pName=pName))
    var(tree, '{pName}_ip'.format(pName=pName)) #TODO Colin clean up hierarchy
    var(tree, '{pName}_ip_signif'.format(pName=pName))
    bookP4(tree, pName)
    
def fillParticle( tree, pName, particle ):
    fill(tree, '{pName}_pdgId'.format(pName=pName), particle.pdgId() )
    ip = -99
    ip_signif = -1e9
    if hasattr(particle, 'path'):
        path = particle.path
        if hasattr(path, 'IP'):
            ip = path.IP
        if hasattr(path, 'IP_signif'):
            ip_signif = path.IP_signif
    fill(tree, '{pName}_ip'.format(pName=pName), ip )
    fill(tree, '{pName}_ip_signif'.format(pName=pName), ip_signif )
    fillP4(tree, pName, particle.p4() )


def bookCluster( tree, name ):
    var(tree, '{name}_e'.format(name=name))
    var(tree, '{name}_layer'.format(name=name))

layers = dict(
    ecal_in = 0,
    hcal_in = 1
)
    
def fillCluster( tree, name, cluster ):
    fill(tree, '{name}_e'.format(name=name), cluster.energy )
    fill(tree, '{name}_layer'.format(name=name), layers[cluster.layer] )
    
# jet

def bookComponent( tree, pName ):
    var(tree, '{pName}_e'.format(pName=pName))
    var(tree, '{pName}_pt'.format(pName=pName))
    var(tree, '{pName}_num'.format(pName=pName))

def fillComponent(tree, pName, component):
    fill(tree, '{pName}_e'.format(pName=pName), component.e() )
    fill(tree, '{pName}_pt'.format(pName=pName), component.pt() )
    fill(tree, '{pName}_num'.format(pName=pName), component.num() )
    
    
pdgIds = [211, 22, 130, 11, 13]
    
def bookJet( tree, pName, taggers=None):
    bookP4(tree, pName )
    for pdgId in pdgIds:
        bookComponent(tree, '{pName}_{pdgId:d}'.format(pName=pName, pdgId=pdgId))
    if taggers:
        for tagger in taggers:
            var(tree, '{pName}_{tagger}'.format(pName=pName, tagger=tagger))


def fillJet( tree, pName, jet, taggers=None):
    fillP4(tree, pName, jet.p4() )
    if taggers:
        for tagger in taggers:
            if tagger in jet.tags:
                fill(tree,
                     '{pName}_{tagger}'.format(pName=pName, tagger=tagger),
                     jet.tags.get(tagger, None))
            else:   
                fill(tree, '{pName}_{tagger}'.format(pName=pName, tagger=tagger), -99)
                
    for pdgId in pdgIds:
        component = jet.constituents.get(pdgId, None)
        if component is not None:
            fillComponent(tree,
                          '{pName}_{pdgId:d}'.format(pName=pName, pdgId=pdgId),
                          component )
        else:
            import pdb; pdb.set_trace()
            print jet
    
def bookMet(tree, pName):
    var(tree, '{pName}_pt'.format(pName=pName)  )
    var(tree, '{pName}_sumet'.format(pName=pName)  )    
    var(tree, '{pName}_phi'.format(pName=pName)  )

def fillMet(tree, pName, met):
    fill(tree, '{pName}_pt'.format(pName=pName), met.pt() )
    fill(tree, '{pName}_sumet'.format(pName=pName), met.sum_et() )
    fill(tree, '{pName}_phi'.format(pName=pName), met.phi() )


