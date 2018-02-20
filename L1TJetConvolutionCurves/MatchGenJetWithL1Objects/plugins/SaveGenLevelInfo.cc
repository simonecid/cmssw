// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/SaveGenLevelInfo
// Class:      SaveGenLevelInfo
// 
/**\class SaveGenLevelInfo SaveGenLevelInfo.cc L1TJetConvolutionCurves/SaveGenLevelInfo/plugins/SaveGenLevelInfo.cc

 Description: Matches a gen jet with a L1T object and creates a tree

 Implementation:
    Receives tags of the various L1T objects
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 12 Jul 2017 14:21:56 GMT
//
//


// system include files
#include <memory>
#include <signal.h>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <utility>


struct Particle {
  unsigned int id;
  float pt, eta, phi;
};

class SaveGenLevelInfo : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit SaveGenLevelInfo(const edm::ParameterSet&);
    ~SaveGenLevelInfo();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetNoNuCollectionTag;
    edm::EDGetTokenT< std::vector< reco::GenParticle > > *_genParticleCollectionTag;

    TTree * _genJetTree;
    TTree * _leadingGenJetTree;
    TTree * _genMuonTree;
    TTree * _leadingGenMuonTree;
    TTree * _genJetNoNuTree;
    TTree * _leadingGenJetNoNuTree;
    TTree * _numberOfGenJetInEventTree;
    TTree * _numberOfGenJetNoNuInEventTree;
    TTree * _numberOfGenMuonInEventTree;

    Particle _genJetParticle;

    unsigned int _numberOfGenJetInEvent;
    unsigned int _numberOfGenMuonInEvent;

};

SaveGenLevelInfo::SaveGenLevelInfo(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  this -> _genJetTree = fs -> make<TTree>("genJetTree", "TTree with generator-level jet information");
  this -> _leadingGenJetTree = fs -> make<TTree>("leadingGenJetTree", "TTree with generator-level jet information");
  this -> _genMuonTree = fs -> make<TTree>("genMuonTree", "TTree with generator-level muon information");
  this -> _leadingGenMuonTree = fs -> make<TTree>("leadingGenMuonTree", "TTree with generator-level muon information");
  this -> _genJetNoNuTree = fs -> make<TTree>("genJetNoNuTree", "TTree with generator-level jet (no nu) information");
  this -> _leadingGenJetNoNuTree = fs -> make<TTree>("leadingGenJetNoNuTree", "TTree with generator-level jet (no nu) information");
  this -> _numberOfGenJetInEventTree = fs -> make<TTree>("numberOfGenJetInEventTree", "TTree with number of gen jets per event");
  this -> _numberOfGenJetNoNuInEventTree = fs -> make<TTree>("numberOfGenJetNoNuInEventTree", "TTree with number of gen jets (no nu) per event");
  this -> _numberOfGenMuonInEventTree = fs -> make<TTree>("numberOfGenMuonInEventTree", "TTree with number of gen muons (no nu) per event");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _genJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _genJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _genJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

  this -> _leadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _leadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _leadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _leadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  
  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genMuonTree -> Branch("genMuon_id", &(this -> _genJetParticle.id), "genMuon_id/i");
  this -> _genMuonTree -> Branch("genMuon_pt", &(this -> _genJetParticle.pt), "genMuon_pt/F");
  this -> _genMuonTree -> Branch("genMuon_eta", &(this -> _genJetParticle.eta), "genMuon_eta/F");
  this -> _genMuonTree -> Branch("genMuon_phi", &(this -> _genJetParticle.phi), "genMuon_phi/F");

  this -> _leadingGenMuonTree -> Branch("genMuon_id", &(this -> _genJetParticle.id), "genMuon_id/i");
  this -> _leadingGenMuonTree -> Branch("genMuon_pt", &(this -> _genJetParticle.pt), "genMuon_pt/F");
  this -> _leadingGenMuonTree -> Branch("genMuon_eta", &(this -> _genJetParticle.eta), "genMuon_eta/F");
  this -> _leadingGenMuonTree -> Branch("genMuon_phi", &(this -> _genJetParticle.phi), "genMuon_phi/F");
  
  this -> _genJetNoNuTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _genJetNoNuTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _genJetNoNuTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _genJetNoNuTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

  this -> _leadingGenJetNoNuTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _leadingGenJetNoNuTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _leadingGenJetNoNuTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _leadingGenJetNoNuTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

  this -> _numberOfGenJetInEventTree -> Branch("numberOfGenJetInEvent", &(this -> _numberOfGenJetInEvent), "numberOfGenJetInEvent/i");
  this -> _numberOfGenJetNoNuInEventTree -> Branch("numberOfGenJetInEvent", &(this -> _numberOfGenJetInEvent), "numberOfGenJetInEvent/i");
  this -> _numberOfGenMuonInEventTree -> Branch("numberOfGenMuonInEvent", &(this -> _numberOfGenMuonInEvent), "numberOfGenMuonInEvent/i");

}

void SaveGenLevelInfo::_getTokens(const edm::ParameterSet& iConfig)
{
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  this -> _genJetNoNuCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetNoNuCollectionTag")));
  this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
}

void SaveGenLevelInfo::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _genJetNoNuCollectionTag) delete this -> _genJetNoNuCollectionTag;
  if (this -> _genParticleCollectionTag) delete this -> _genParticleCollectionTag;
}

SaveGenLevelInfo::~SaveGenLevelInfo()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
SaveGenLevelInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  edm::Handle < std::vector< reco::GenJet > > genJetNoNuCollectionHandle;
  edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);
  iEvent.getByToken(*(this -> _genJetNoNuCollectionTag), genJetNoNuCollectionHandle);
  iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);

  // Saving every genJet to get the misidentification probability
  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    this -> _genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
    this -> _genJetParticle.pt = genJetIterator -> pt();
    //std::cout << "Jet pt " << genJetIterator -> pt() << " eta " << genJetIterator -> eta() << std::endl;
    this -> _genJetParticle.eta = genJetIterator -> eta();
    this -> _genJetParticle.phi = genJetIterator -> phi();
    this -> _genJetTree -> Fill();
  }

  // Saving every genJet to get the misidentification probability
  float maxPt = 0;
  bool save = false;
  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    if (genJetIterator->pt() > maxPt) 
    {      
      this -> _genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
      this -> _genJetParticle.pt = genJetIterator -> pt();
      this -> _genJetParticle.eta = genJetIterator -> eta();
      this -> _genJetParticle.phi = genJetIterator -> phi();
      save = true;
      maxPt = genJetIterator -> pt();
    }
  }
  if (save)
  {
    //std::cout << "Saving jet pt " << this -> _genJetParticle.pt << " eta " << this -> _genJetParticle.eta << std::endl;
    this -> _leadingGenJetTree -> Fill();
  }

  this -> _numberOfGenJetInEvent = genJetCollectionHandle -> size();
  this -> _numberOfGenJetInEventTree -> Fill();
  
  // Saving every genJet to get the misidentification probability
  for (auto genJetIterator = genJetNoNuCollectionHandle -> begin(); genJetIterator != genJetNoNuCollectionHandle -> end(); genJetIterator++ )
  {
    this -> _genJetParticle.id = (genJetIterator - genJetNoNuCollectionHandle->begin());
    this -> _genJetParticle.pt = genJetIterator -> pt();
    //std::cout << "Jet pt " << genJetIterator -> pt() << " eta " << genJetIterator -> eta() << std::endl;
    this -> _genJetParticle.eta = genJetIterator -> eta();
    this -> _genJetParticle.phi = genJetIterator -> phi();
    this -> _genJetNoNuTree -> Fill();
  }

  maxPt = 0;
  save = false;
  // Saving every genJet to get the misidentification probability
  for (auto genJetIterator = genJetNoNuCollectionHandle -> begin(); genJetIterator != genJetNoNuCollectionHandle -> end(); genJetIterator++ )
  {
    if (genJetIterator->pt() > maxPt) 
    {      
      this -> _genJetParticle.id = (genJetIterator - genJetNoNuCollectionHandle->begin());
      this -> _genJetParticle.pt = genJetIterator -> pt();
      this -> _genJetParticle.eta = genJetIterator -> eta();
      this -> _genJetParticle.phi = genJetIterator -> phi();
      save = true;
      maxPt = genJetIterator -> pt();
    }
  }
  if (save)
  {
    //std::cout << "Saving jet pt " << this -> _genJetParticle.pt << " eta " << this -> _genJetParticle.eta << std::endl;
    this -> _leadingGenJetNoNuTree -> Fill();
  }

  this -> _numberOfGenJetInEvent = genJetNoNuCollectionHandle -> size();
  this -> _numberOfGenJetNoNuInEventTree -> Fill();

  const reco::GenParticle* leadingGenMuon = NULL;
  maxPt = 0;
  save = false;
  this -> _numberOfGenMuonInEvent = 0;
  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    if (abs(genParticleIterator->pdgId()) == 13) {

      this -> _genJetParticle.id = this -> _numberOfGenMuonInEvent;
      (this -> _numberOfGenMuonInEvent)++;
      this -> _genJetParticle.pt = genParticleIterator -> pt();
      this -> _genJetParticle.eta = genParticleIterator -> eta();
      this -> _genJetParticle.phi = genParticleIterator -> phi();
      //std::cout << "Muon pt " << genParticleIterator -> pt() << " eta " << genParticleIterator -> eta() << std::endl;
      this -> _genMuonTree -> Fill();

      if (this -> _genJetParticle.pt > maxPt) 
      {
        maxPt = this -> _genJetParticle.pt;
        leadingGenMuon = &(*(genParticleIterator));
        save = true;
      }
    }
      
  }

  if (save)
  {
    this -> _genJetParticle.id = 0;
    this -> _genJetParticle.pt = leadingGenMuon -> pt();
    this -> _genJetParticle.eta = leadingGenMuon -> eta();
    this -> _genJetParticle.phi = leadingGenMuon -> phi();
    //std::cout << "Saving muon pt " << this -> _genJetParticle.pt << " eta " << this -> _genJetParticle.eta << std::endl;
    this -> _leadingGenMuonTree -> Fill();    
  }

  this -> _numberOfGenMuonInEventTree -> Fill();
  

}

// ------------ method called once each job just before starting event loop  ------------
void SaveGenLevelInfo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void SaveGenLevelInfo::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SaveGenLevelInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(SaveGenLevelInfo);
