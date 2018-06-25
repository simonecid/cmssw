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

#include "L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchingAlgorithms.h"

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

struct TriggerObject {
  unsigned int id;
  float pt, eta, phi;
  int hwQual;
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

    template <class TParticle, class TTrigger> // <3
    const std::vector <std::tuple<const TTrigger*, const TParticle*, float, int> > 
    _matchL1TMuonWithGenLevelMuons
    (
      const edm::Handle<std::vector<TParticle>>& particleCollectionHandle,
      const edm::Handle<BXVector<TTrigger>>& l1tObjectCollectionHandle,
      float dr2Min = 25,
      bool crossMatch = true
    );

    edm::EDGetTokenT< std::vector< reco::GenParticle > > *_genParticleCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
    Particle _genParticle;
    TriggerObject _l1tObjectParticle;
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
  this -> _l1tMuonGenParticleTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
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

template <class TParticle, class TTrigger> // <3
const std::vector <std::tuple<const TTrigger*, const TParticle*, float, int> > 
MatchL1TMuonWithGenLevelMuons::_matchL1TMuonWithGenLevelMuons
(
  const edm::Handle<std::vector<TParticle>>& particleCollectionHandle,
  const edm::Handle<BXVector<TTrigger>>& l1tObjectCollectionHandle,
  float dr2Min,
  bool crossMatch
)
{
  std::vector <std::tuple<const TTrigger*, const TParticle*, float, int> > l1tObjectParticlePairs;
  // for each object in the particle collection we look for the closest l1tobject in a wide range
  for (auto particleIterator = particleCollectionHandle -> begin(); particleIterator != particleCollectionHandle -> end(); particleIterator++ )
  {
    std::tuple<const TTrigger *, const TParticle *, float, int> l1tObjectParticlePair = 
      MatchingAlgorithms::matchParticleWithL1Object<>
      (
        *particleIterator,
        particleCollectionHandle,
        l1tObjectCollectionHandle,
        dr2Min,
        crossMatch
      );
    // if we have found a match let's add it.
    if (std::get<0>(l1tObjectParticlePair) != NULL) l1tObjectParticlePairs.push_back(l1tObjectParticlePair);
  }
  // return the matched stuff
  return l1tObjectParticlePairs;
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
  auto l1tMuonGenMuonPairs = this -> _matchL1TMuonWithGenLevelMuons<>(genParticleCollectionHandle, l1tMuonCollectionHandle, 25, true);

  for (const auto & matchTuple : l1tMuonGenMuonPairs)
  {

    const l1t::Muon* matchedL1TMuon = std::get<0>(matchTuple);
    const reco::GenParticle* matchedGenMuon = std::get<1>(matchTuple);
    float deltaR2 = std::get<2>(matchTuple);

    this -> _genParticle.id = 0;
    this -> _genParticle.pt = matchedGenMuon -> pt();
    this -> _genParticle.eta = matchedGenMuon -> eta();
    this -> _genParticle.phi = matchedGenMuon -> phi();
    this -> _l1tObjectParticle.id = 0;
    this -> _l1tObjectParticle.pt = matchedL1TMuon -> pt();
    this -> _l1tObjectParticle.eta = matchedL1TMuon -> eta();
    this -> _l1tObjectParticle.phi = matchedL1TMuon -> phi();
    this -> _l1tObjectParticle.hwQual = matchedL1TMuon -> hwQual();
    this -> _deltaR2 = deltaR2;
    this -> _l1tMuonGenParticleTree -> Fill();
  }

  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    // Only the muonssss
    if (abs(genParticleIterator -> pdgId()) != 13) continue;
    std::cout << "Gen-Muon object" << std::endl;
    this -> _genParticle.id = (genParticleIterator - genParticleCollectionHandle -> begin());
    this -> _genParticle.pt = genParticleIterator -> pt();
    this -> _genParticle.phi = genParticleIterator -> phi();
    this -> _genParticle.eta = genParticleIterator -> eta();
    this -> _genParticleTree -> Fill();
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
