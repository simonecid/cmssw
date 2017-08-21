// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchL1TMuonWithGenLevelMuons
// Class:      MatchL1TMuonWithGenLevelMuons
// 
/**\class MatchL1TMuonWithGenLevelMuons MatchL1TMuonWithGenLevelMuons.cc L1TJetConvolutionCurves/MatchL1TMuonWithGenLevelMuons/plugins/MatchL1TMuonWithGenLevelMuons.cc

 Description: [one line class summary]

 Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Simone Bologna
//         Created:  Fri, 18 Aug 2017 14:10:47 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <cmath>
#include "TTree.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

struct Particle {
  unsigned int id;
  float pt, eta, phi;
};

class MatchL1TMuonWithGenLevelMuons : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit MatchL1TMuonWithGenLevelMuons(const edm::ParameterSet&);
    ~MatchL1TMuonWithGenLevelMuons();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector < std::tuple < Particle, Particle, float > > 
    _matchL1TMuonWithGenLevelMuons
    (
      const edm::Handle< std::vector< reco::GenParticle > > &,
      const edm::Handle < BXVector < T > > & 
    );

    edm::EDGetTokenT< std::vector< reco::GenParticle > > *_genParticleCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
    Particle _genParticle;
    Particle _l1tObjectParticle;
    float _deltaR2;

    TTree * _l1tMuonGenParticleTree;
    TTree * _genParticleTree;
    
    // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MatchL1TMuonWithGenLevelMuons::MatchL1TMuonWithGenLevelMuons(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  this -> _l1tMuonCollectionTag = new edm::EDGetTokenT< BXVector< l1t::Muon > >(consumes< BXVector< l1t::Muon > > (iConfig.getParameter< edm::InputTag >("l1tMuonCollectionTag")));  

  this -> _l1tMuonGenParticleTree = fs -> make<TTree>("matchedL1TMuonGenParticleTree", "TTree with generator-level Muon / L1T Muon information");
  this -> _genParticleTree = fs -> make<TTree>("genParticleTree", "TTree with generator-level muon information");

  this -> _l1tMuonGenParticleTree -> Branch("genParticle_id", &(this -> _genParticle.id), "genParticle_id/i");
  this -> _l1tMuonGenParticleTree -> Branch("genParticle_pt", &(this -> _genParticle.pt), "genParticle_pt/F");
  this -> _l1tMuonGenParticleTree -> Branch("genParticle_eta", &(this -> _genParticle.eta), "genParticle_eta/F");
  this -> _l1tMuonGenParticleTree -> Branch("genParticle_phi", &(this -> _genParticle.phi), "genParticle_phi/F");
  this -> _l1tMuonGenParticleTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tMuonGenParticleTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonGenParticleTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonGenParticleTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonGenParticleTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genParticleTree -> Branch("genParticle_id", &(this -> _genParticle.id), "genParticle_id/i");
  this -> _genParticleTree -> Branch("genParticle_pt", &(this -> _genParticle.pt), "genParticle_pt/F");
  this -> _genParticleTree -> Branch("genParticle_eta", &(this -> _genParticle.eta), "genParticle_eta/F");
  this -> _genParticleTree -> Branch("genParticle_phi", &(this -> _genParticle.phi), "genParticle_phi/F");
}


MatchL1TMuonWithGenLevelMuons::~MatchL1TMuonWithGenLevelMuons()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  if (this -> _genParticleCollectionTag) delete this -> _genParticleCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
  
}

template <class T> // <3
const std::vector <std::tuple < Particle, Particle, float > > 
MatchL1TMuonWithGenLevelMuons::_matchL1TMuonWithGenLevelMuons
(
  const edm::Handle < std::vector< reco::GenParticle > > & genParticleCollectionHandle,
  const edm::Handle < BXVector < T > > & l1tMuonCollectionHandle 
)
{

  std::vector< std::tuple < Particle , Particle , float > > l1tObjectGenParticlePairs;
  // for each object in the l1t collection we look for the closest jet in a wide range
  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    // Only the muonssss
    if (abs(genParticleIterator -> pdgId()) != 13) continue;

    bool foundMatch = false;
    float dr2Min = 25; // i.e. dr = 5
    std::tuple<Particle, Particle, float > l1tObjectGenParticlePair;
    
    Particle & matchedGenParticle = std::get<1>(l1tObjectGenParticlePair);
    matchedGenParticle.id = (genParticleIterator - genParticleCollectionHandle -> begin());
    matchedGenParticle.pt = genParticleIterator -> pt();
    matchedGenParticle.phi = genParticleIterator -> phi();
    matchedGenParticle.eta = genParticleIterator -> eta();
    
    for (
      typename BXVector<T>::const_iterator bx0Iterator = l1tMuonCollectionHandle -> begin(0);
      bx0Iterator != l1tMuonCollectionHandle -> end(0);
      bx0Iterator++
    )
    {
      float dr2 = reco::deltaR2(*bx0Iterator, *genParticleIterator);
      if (dr2 < dr2Min)
      {
        std::get<2>(l1tObjectGenParticlePair) = dr2;
        dr2Min = dr2;
        foundMatch = true;
        Particle & matchedL1TObject = std::get<0>(l1tObjectGenParticlePair);
        matchedL1TObject.id = (bx0Iterator - l1tMuonCollectionHandle -> begin(0));
        matchedL1TObject.pt = bx0Iterator -> pt();
        matchedL1TObject.phi = bx0Iterator -> phi();
        matchedL1TObject.eta = bx0Iterator -> eta();
      }
    }
    
    //if we have found a compatible gen jet we push the object-jet pair in a vector which will be our result
    if (foundMatch)
      l1tObjectGenParticlePairs.push_back(l1tObjectGenParticlePair);
  }

  return l1tObjectGenParticlePairs;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchL1TMuonWithGenLevelMuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //Retrieving muon and gen particle data
  edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
  iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);
  edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
  iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);

  // Matching
  auto l1tMuonGenMuonPairs = this -> _matchL1TMuonWithGenLevelMuons<>(genParticleCollectionHandle, l1tMuonCollectionHandle);

  //Debug printout

    for (const auto & matchTuple : l1tMuonGenMuonPairs) 
  {
    const Particle & l1tMuon = std::get<0>(matchTuple);
    const Particle & genMuon = std::get<1>(matchTuple);
    float deltaR2 = std::get<2>(matchTuple);
    std::cout << "MATCHED MUON:" << std::endl;
    std::cout << "Gen-level id: " << genMuon.id << "\t L1T-level id: " << l1tMuon.id << std::endl;
    std::cout << "Gen-level pt: " << genMuon.pt << "\t L1T-level pt: " << l1tMuon.pt << std::endl;
    std::cout << "Gen-level eta: " << genMuon.eta << "\t L1T-level eta: " << l1tMuon.eta << std::endl;
    std::cout << "Gen-level phi: " << genMuon.phi << "\t L1T-level phi: " << l1tMuon.phi << std::endl;
    std::cout << "deltaR2: " << deltaR2 << std::endl;
  }

  //Saving all the gen muons

  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    this -> _genParticle.id = (genParticleIterator - genParticleCollectionHandle->begin());
    this -> _genParticle.pt = genParticleIterator -> pt();
    this -> _genParticle.eta = genParticleIterator -> eta();
    this -> _genParticle.phi = genParticleIterator -> phi();
    this -> _genParticleTree -> Fill();
  }

  //Saving all the matched muon-gen muon pairs
  for (const auto & matchTuple : l1tMuonGenMuonPairs)
  {
    const Particle & l1tObject = std::get<0>(matchTuple);
    const Particle & genMuon = std::get<1>(matchTuple);
    float deltaR2 = std::get<2>(matchTuple);
    this -> _genParticle.id = genMuon.id;
    this -> _genParticle.pt = genMuon.pt;
    this -> _genParticle.eta = genMuon.eta;
    this -> _genParticle.phi = genMuon.phi;
    this -> _l1tObjectParticle.id = l1tObject.id;
    this -> _l1tObjectParticle.pt = l1tObject.pt;
    this -> _l1tObjectParticle.eta = l1tObject.eta;
    this -> _l1tObjectParticle.phi = l1tObject.phi;
    this -> _deltaR2 = deltaR2;
    this -> _l1tMuonGenParticleTree -> Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MatchL1TMuonWithGenLevelMuons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MatchL1TMuonWithGenLevelMuons::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MatchL1TMuonWithGenLevelMuons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchL1TMuonWithGenLevelMuons);
