// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchGenJetWithL1Objects
// Class:      MatchGenJetWithL1Objects
// 
/**\class MatchGenJetWithL1Objects MatchGenJetWithL1Objects.cc L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchGenJetWithL1Objects.cc

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

class MatchGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MatchGenJetWithL1Objects(const edm::ParameterSet&);
    ~MatchGenJetWithL1Objects();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector< std::tuple < TriggerObject, Particle, float> >
    _matchGenJetWithL1Object
    (
      const edm::Handle<std::vector<reco::GenJet>> &,
      const edm::Handle<BXVector<T>> &,
      float = 0.25
    );

    template <class TParticle, class TTrigger> // <3
    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector < std::tuple < const TTrigger*, const TParticle*, float, int > > &
    );

    template <class TParticle, class TTrigger> // <3
    const std::vector< std::tuple<const TTrigger*, const TParticle*, float, int> >
    _matchParticleWithL1Object
    (
      const edm::Handle<std::vector<TParticle>> &,
      const edm::Handle<BXVector<TTrigger>> &,
      float = 0.25,
      bool = true
    );

    edm::EDGetTokenT<std::vector<reco::GenParticle>> *_genMuonCollectionTag;
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
 
    TTree * _l1tMuonGenJetGenMuonTree;
    Particle _genJetParticle;
    Particle _genMuonParticle;
    TriggerObject _l1tObjectParticle;
    float _deltaR2_genJet_l1tMuon;
    float _deltaR2_genMuon_l1tMuon;
};

MatchGenJetWithL1Objects::MatchGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tMuonGenJetGenMuonTree = fs -> make<TTree>("matchedL1TMuonGenJetTree", "TTree with generator-level jet & muon + L1T Muon information");

  this -> _l1tMuonGenJetGenMuonTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/I");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genMuon_id", &(this -> _genMuonParticle.id), "genMuon_id/I");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genMuon_pt", &(this -> _genMuonParticle.pt), "genMuon_pt/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genMuon_eta", &(this -> _genMuonParticle.eta), "genMuon_eta/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("genMuon_phi", &(this -> _genMuonParticle.phi), "genMuon_phi/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/I");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("deltaR2_genJet_l1tMuon", &(this -> _deltaR2_genJet_l1tMuon), "_deltaR2_genJet_l1tMuon/F");
  this -> _l1tMuonGenJetGenMuonTree -> Branch("deltaR2_genMuon_l1tMuon", &(this -> _deltaR2_genMuon_l1tMuon), "_deltaR2_genMuon_l1tMuon/F");


}

void MatchGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _genMuonCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> genParticleCollectionTag not found." << std::endl;
    this -> _genMuonCollectionTag = NULL;
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

void MatchGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _genMuonCollectionTag) delete this -> _genMuonCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
}

MatchGenJetWithL1Objects::~MatchGenJetWithL1Objects()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);
  
  edm::Handle < std::vector< reco::GenParticle > > genMuonCollectionHandle;
  iEvent.getByToken(*(this -> _genMuonCollectionTag), genMuonCollectionHandle);

  edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
  iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);

  // Stop if no l1t muon

  if (l1tMuonCollectionHandle -> begin(0) == l1tMuonCollectionHandle -> end(0))
    return;

  std::vector<std::tuple<const l1t::Muon*,const  reco::GenParticle*, float, int> > 
    l1tMuonGenMuonPairs;

  std::vector<std::tuple<const l1t::Muon*,const  reco::GenJet*, float, int> > 
    l1tMuonGenJetPairs;

  for (auto l1tMuonIterator = l1tMuonCollectionHandle -> begin(); l1tMuonIterator != l1tMuonCollectionHandle -> end(); l1tMuonIterator++ )
  {
    std::tuple<const l1t::Muon *, const reco::GenParticle *, float, int> l1tMuonGenMuonPair = 
      MatchingAlgorithms::matchL1ObjectWithParticle<>
      (
        *l1tMuonIterator,
        l1tMuonCollectionHandle,
        genMuonCollectionHandle,
        5,
        true
      );

    std::tuple<const l1t::Muon *, const reco::GenJet *, float, int> l1tMuonGenJetPair = 
      MatchingAlgorithms::matchL1ObjectWithParticle<>
      (
        *l1tMuonIterator,
        l1tMuonCollectionHandle,
        genJetCollectionHandle,
        5,
        true
      );

    // if we have found a match let's add it.
    if (std::get<1>(l1tMuonGenMuonPair) != NULL) 
    {
      const reco::GenParticle* genMuon = std::get<1>(l1tMuonGenMuonPair);
      this -> _genMuonParticle.id = int(genMuon - &(*genMuonCollectionHandle -> begin()));
      this -> _genMuonParticle.pt = genMuon -> pt();
      this -> _genMuonParticle.eta = genMuon -> eta();
      this -> _genMuonParticle.phi = genMuon -> phi();
    }
    else 
    {
      this -> _genMuonParticle.id = -1;
      this -> _genMuonParticle.pt = -1;
      this -> _genMuonParticle.eta = -1;
      this -> _genMuonParticle.phi = -1;
    }
      
    if (std::get<1>(l1tMuonGenJetPair) != NULL) 
    {
      const reco::GenJet* genJet = std::get<1>(l1tMuonGenJetPair);
      this -> _genJetParticle.id = int(genJet - &(*genJetCollectionHandle -> begin()));
      this -> _genJetParticle.pt = genJet -> pt();
      this -> _genJetParticle.eta = genJet -> eta();
      this -> _genJetParticle.phi = genJet -> phi();
    }
    else
    {
      this -> _genJetParticle.id = -1;
      this -> _genJetParticle.pt = -1;
      this -> _genJetParticle.eta = -1;
      this -> _genJetParticle.phi = -1;
    }

    this -> _l1tObjectParticle.id = int(l1tMuonIterator - l1tMuonCollectionHandle -> begin(0));
    this -> _l1tObjectParticle.pt = l1tMuonIterator -> pt();
    this -> _l1tObjectParticle.eta = l1tMuonIterator -> eta();
    this -> _l1tObjectParticle.phi = l1tMuonIterator -> phi();
    this -> _l1tObjectParticle.hwQual = l1tMuonIterator -> hwQual();

    this -> _deltaR2_genJet_l1tMuon = std::get<2>(l1tMuonGenJetPair);
    this -> _deltaR2_genMuon_l1tMuon = std::get<2>(l1tMuonGenMuonPair);

    this -> _l1tMuonGenJetGenMuonTree -> Fill();

    return;

  }  
}
  
  
template<class TParticle, class TTrigger>
void
MatchGenJetWithL1Objects::_fillTreeWithMatchedPairs
(
  TTree & aTree,
  const std::vector < std::tuple <const TTrigger*, const TParticle*, float, int > > & matchedL1TObjectJetPairs
)
{
  for (const auto & matchTuple : matchedL1TObjectJetPairs)
  {
    const TTrigger & l1tObject = *(std::get<0>(matchTuple));
    const TParticle & genJet = *(std::get<1>(matchTuple));
    float deltaR2 = std::get<2>(matchTuple);
    float matchingQuality = std::get<3>(matchTuple);
    this -> _genJetParticle.id = 0;
    this -> _genJetParticle.pt = genJet.pt();
    this -> _genJetParticle.eta = genJet.eta();
    this -> _genJetParticle.phi = genJet.phi();
    this -> _l1tObjectParticle.id = 0;
    this -> _l1tObjectParticle.pt = l1tObject.pt();
    this -> _l1tObjectParticle.eta = l1tObject.eta();
    this -> _l1tObjectParticle.phi = l1tObject.phi();
    this -> _l1tObjectParticle.hwQual = l1tObject.hwQual();
    //this -> _deltaR2 = deltaR2;
    //this -> _matchingQuality = matchingQuality;
    
    aTree.Fill();
    
  }
}

template <class TParticle, class TTrigger> // <3
const std::vector <std::tuple<const TTrigger*, const TParticle*, float, int> > 
MatchGenJetWithL1Objects::_matchParticleWithL1Object
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

// Considers a genJet and looks for the closest l1tObject within that radius (1 genjet can have only one l1tobject)

template <class T> // <3
const std::vector <std::tuple < TriggerObject, Particle, float > > 
MatchGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const edm::Handle < std::vector< reco::GenJet > > & genJetCollectionHandle,
  const edm::Handle < BXVector < T > > & l1tObjectCollectionHandle,
  float dr2Min
)
{

  std::vector< std::tuple < TriggerObject, Particle , float > > l1tObjectGenJetPairs;
  // for each object in the genJet collection we look for the closest l1tobject in a wide range
  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {    
    bool foundMatch = false;
    std::tuple<TriggerObject, Particle, float > l1tObjectGenJetPair;
    
    Particle & matchedGenJet = std::get<1>(l1tObjectGenJetPair);
    matchedGenJet.id = (genJetIterator - genJetCollectionHandle -> begin());
    matchedGenJet.pt = genJetIterator -> pt();
    matchedGenJet.phi = genJetIterator -> phi();
    matchedGenJet.eta = genJetIterator -> eta();

    for (
          typename BXVector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle->begin(0);
          bx0Iterator != l1tObjectCollectionHandle->end(0);
          bx0Iterator++
    )
    {
      float dr2 = reco::deltaR2(*bx0Iterator, *genJetIterator);

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
  }

  return l1tObjectGenJetPairs;
}

// ------------ method called once each job just before starting event loop  ------------
void MatchGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchGenJetWithL1Objects);
