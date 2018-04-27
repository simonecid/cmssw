// -*- C++ -*-
//
// Package:    l1tMuonConvolutionCurves/MatchLeadingGenJetWithL1Objects
// Class:      MatchLeadingGenJetWithL1Objects
// 
/**\class MatchLeadingGenJetWithL1Objects MatchLeadingGenJetWithL1Objects.cc l1tMuonConvolutionCurves/MatchLeadingGenJetWithL1Objects/plugins/MatchLeadingGenJetWithL1Objects.cc

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
#include <csignal>

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

#include "L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchingAlgorithms.h"

struct Particle {
  unsigned int id;
  float pt, eta, phi;
};

struct TriggerObject {
  unsigned int id;
  float pt, eta, phi;
  int hwQual;
};

class MatchLeadingGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MatchLeadingGenJetWithL1Objects(const edm::ParameterSet&);
    ~MatchLeadingGenJetWithL1Objects();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector< std::tuple < TriggerObject, Particle, float, int> >
    _matchGenJetWithL1Object(
      const reco::GenJet &,
      const edm::Handle<BXVector<T>> &,
      float = 0.25,
      bool = true
    );

    template <class TParticle, class TTrigger> // <3
    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector<std::tuple<const TTrigger*, const TParticle*, float, int> > &
    );

    edm::EDGetTokenT<std::vector<reco::GenParticle>> *_genParticleCollectionTag;
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;

    TTree * _l1tMuonLeadingGenJetTree;
    TTree * _leadingGenJetTree;

    Particle _genJetParticle;
    TriggerObject _l1tObjectParticle;
    int _hwQual;
    float _deltaR2;
    int _matchingQuality;
};

MatchLeadingGenJetWithL1Objects::MatchLeadingGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tMuonLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TMuonTree", "TTree with generator-level jet / L1T Muon information");
  this -> _leadingGenJetTree = fs -> make<TTree>("leadingGenJetTree", "TTree with generator-level jet information in detector space");

  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _leadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _leadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _leadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _leadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

}

void MatchLeadingGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> genParticleCollectionTag not found." << std::endl;
    this -> _genParticleCollectionTag = NULL;
  }

  try
  {
    this -> _l1tMuonCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Muon > >(consumes< BXVector< l1t::Muon > > (iConfig.getParameter< edm::InputTag >("l1tMuonCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tMuonCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tMuonCollectionTag = NULL;
  }
  
  return;
}

void MatchLeadingGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _genParticleCollectionTag) delete this -> _genParticleCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
}

MatchLeadingGenJetWithL1Objects::~MatchLeadingGenJetWithL1Objects()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchLeadingGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);

  //Looking for the leading gen jet

  float maxPt = 0;
  bool save = false;
  const reco::GenJet* leadingGenJetInDetector = NULL;

  for (auto genJetIterator = genJetCollectionHandle->begin(); genJetIterator != genJetCollectionHandle->end(); genJetIterator++)
  {
    // Only detector leading l1t jets
    if ((genJetIterator->pt() > maxPt) && (genJetIterator->eta() < 2.44) && (genJetIterator->eta() > -2.44))
    {
      save = true;
      maxPt = genJetIterator->pt();
      leadingGenJetInDetector = &(*genJetIterator);
      this->_genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
      this->_genJetParticle.pt = genJetIterator->pt();
      this->_genJetParticle.eta = genJetIterator->eta();
      this->_genJetParticle.phi = genJetIterator->phi();
    }
  }
  if (save)
    this->_leadingGenJetTree->Fill();

  maxPt = 0;
  save = false;

  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes
  // Muon veto. For muon matching
  edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
  iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);
  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    if (abs(genParticleIterator->pdgId()) == 13)
      return;
  }

  if (leadingGenJetInDetector != NULL){
    if (this -> _l1tMuonCollectionTag)
    {

      edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
      iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);
      std::tuple<const l1t::Muon *, const reco::GenJet *, float, int> l1tMuonGenJetPair = MatchingAlgorithms::matchParticleWithL1Object<>(
          *leadingGenJetInDetector,
          genJetCollectionHandle,
          l1tMuonCollectionHandle,
          5,
          true);
      std::vector<std::tuple<const l1t::Muon *, const reco::GenJet *, float, int>> l1tMuonGenJetPairs;
      if (std::get<0>(l1tMuonGenJetPair) != NULL) l1tMuonGenJetPairs.push_back(l1tMuonGenJetPair);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tMuonLeadingGenJetTree), l1tMuonGenJetPairs);
    }
    
  }
}

template <class TParticle, class TTrigger> // <3
void MatchLeadingGenJetWithL1Objects::_fillTreeWithMatchedPairs(
    TTree &aTree,
    const std::vector< std::tuple<const TTrigger*,const  TParticle*, float, int>> &matchedL1TObjectJetPairs)
{
  for (const auto & matchTuple : matchedL1TObjectJetPairs)
  {
    const TTrigger & l1tObject = *(std::get<0>(matchTuple));
    const TParticle & genJet = *(std::get<1>(matchTuple));
    float deltaR2 = std::get<2>(matchTuple);
    int matchingQuality = std::get<3>(matchTuple);
    this -> _genJetParticle.id = 0;
    this -> _genJetParticle.pt = genJet.pt();
    this -> _genJetParticle.eta = genJet.eta();
    this -> _genJetParticle.phi = genJet.phi();
    this -> _l1tObjectParticle.id = 0;
    this -> _l1tObjectParticle.pt = l1tObject.pt();
    this -> _l1tObjectParticle.eta = l1tObject.eta();
    this -> _l1tObjectParticle.phi = l1tObject.phi();
    this -> _l1tObjectParticle.hwQual = l1tObject.hwQual();
    this -> _deltaR2 = deltaR2;
    this -> _matchingQuality = matchingQuality;
    aTree.Fill();
  }
}

// Considers a genJet and looks for the closest l1tObject within that radius (1 genjet can have only one l1tobject)

template <class T> // <3
const std::vector <std::tuple < TriggerObject, Particle, float, int > > 
MatchLeadingGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const reco::GenJet & genJet,
  const edm::Handle < BXVector < T > > & l1tObjectCollectionHandle,
  float dr2Min,
  bool crossMatch
)
{

  std::vector< std::tuple < TriggerObject, Particle , float, int > > l1tObjectGenJetPairs;
  // for each object in the genJet collection we look for the closest l1tobject in a wide range
  
  bool foundMatch = false;
  std::tuple<TriggerObject, Particle, float, int > l1tObjectGenJetPair;
  
  Particle & matchedGenJet = std::get<1>(l1tObjectGenJetPair);
  matchedGenJet.id = 0;
  matchedGenJet.pt = genJet.pt();
  matchedGenJet.phi = genJet.phi();
  matchedGenJet.eta = genJet.eta();

  for (
        typename BXVector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle->begin(0);
        bx0Iterator != l1tObjectCollectionHandle->end(0);
        bx0Iterator++
  )
  {
    float dr2 = reco::deltaR2(*bx0Iterator, genJet);

    if ((dr2 < dr2Min))
    {
      TriggerObject &matchedL1TObject = std::get<0>(l1tObjectGenJetPair);
      matchedL1TObject.id = (bx0Iterator - l1tObjectCollectionHandle->begin(0));
      matchedL1TObject.pt = bx0Iterator->pt();
      matchedL1TObject.phi = bx0Iterator->phi();
      matchedL1TObject.eta = bx0Iterator->eta();
      matchedL1TObject.hwQual = bx0Iterator->hwQual();

      std::get<2>(l1tObjectGenJetPair) = dr2;
      dr2Min = dr2;
      foundMatch = true;
    }
  
  } 
  //if we have found a compatible l1t object we push the object-jet pair in a vector which will be our result
  if (foundMatch)
    l1tObjectGenJetPairs.push_back(l1tObjectGenJetPair);
  

  return l1tObjectGenJetPairs;
}

// ------------ method called once each job just before starting event loop  ------------
void MatchLeadingGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchLeadingGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchLeadingGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchLeadingGenJetWithL1Objects);
