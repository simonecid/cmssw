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

    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Jet > > *_l1tJetCollectionTag;
    //!TODO: Like this

    TTree * _l1tJetGenJetTree;
    TTree * _l1tJetTree;
    TTree * _genJetTree;

    Particle _genJetParticle;
    TriggerObject _l1tObjectParticle;
    int _hwQual;
    float _deltaR2;
    int _matchingQuality;
};

MatchGenJetWithL1Objects::MatchGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tJetGenJetTree = fs -> make<TTree>("matchedL1TJetGenJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _genJetTree = fs -> make<TTree>("genJetTree", "TTree with generator-level jet information");

  this -> _l1tJetGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tJetGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tJetGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tJetGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");
  this -> _l1tJetGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _l1tJetGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _genJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _genJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _genJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

  this -> _l1tJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");

}

void MatchGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> genJetCollectionTag not found." << std::endl;
    this -> _genJetCollectionTag = NULL;
  }

  try
  {
    this -> _l1tJetCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Jet > >(consumes< BXVector< l1t::Jet > > (iConfig.getParameter< edm::InputTag >("l1tJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tJetCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tJetCollectionTag = NULL;
  }
  
  return;
}

void MatchGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _genParticleCollectionTag) delete this -> _genParticleCollectionTag;
  if (this -> _l1tCaloTowerCollectionTag) delete this -> _l1tCaloTowerCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
  if (this -> _l1tJetCollectionTag) delete this -> _l1tJetCollectionTag;
  if (this -> _l1tEGammaCollectionTag) delete this -> _l1tEGammaCollectionTag;
  if (this -> _l1tTauCollectionTag) delete this -> _l1tTauCollectionTag;
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

  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes
  
  if (this -> _l1tMuonCollectionTag)
  {
    // Muon veto. For muon matching
    edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
    iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);
    bool hasMuons = false;
    for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
    {
      if (abs(genParticleIterator->pdgId()) == 13)
        hasMuons = true;
    }

    if (!hasMuons)
    {
      edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
      iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);
      const std::vector<std::tuple<const l1t::Muon*,const  reco::GenJet*, float, int> > 
        l1tMuonGenJetPairs = 
          this -> _matchParticleWithL1Object<>(genJetCollectionHandle, l1tMuonCollectionHandle, 5, true);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tMuonGenJetTree), l1tMuonGenJetPairs);
      float maxPt = 0;
      bool save = false;
      for (auto l1tMuonIterator = l1tMuonCollectionHandle -> begin(0); l1tMuonIterator != l1tMuonCollectionHandle -> end(0); l1tMuonIterator++ )
      {
        if (l1tMuonIterator -> hwQual() < 8) continue;
        if (l1tMuonIterator -> pt() > maxPt)
        {
          save = true;
          maxPt = l1tMuonIterator -> pt();
          this -> _l1tObjectParticle.hwQual = l1tMuonIterator -> hwQual();      
          this -> _l1tObjectParticle.id = (l1tMuonIterator - l1tMuonCollectionHandle->begin(0));
          this -> _l1tObjectParticle.pt = l1tMuonIterator -> pt();
          this -> _l1tObjectParticle.eta = l1tMuonIterator -> eta();
          this -> _l1tObjectParticle.phi = l1tMuonIterator -> phi();
        }
      }
      if (save) {
        //std::cout << "I am saving a muon object w/ qual " << this -> _l1tObjectParticle.hwQual << std::endl;      
        this -> _l1tLeadingMuonTree -> Fill();
      }
      
      for (auto l1tMuonIterator = l1tMuonCollectionHandle -> begin(0); l1tMuonIterator != l1tMuonCollectionHandle -> end(0); l1tMuonIterator++ )
      {
        this -> _l1tObjectParticle.id = (l1tMuonIterator - l1tMuonCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tMuonIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tMuonIterator -> eta();
        this -> _l1tObjectParticle.phi = l1tMuonIterator -> phi();
        this -> _l1tObjectParticle.hwQual = l1tMuonIterator -> hwQual();      
        this -> _l1tMuonTree -> Fill();
      }
      
    }
  }
  
  if (this -> _l1tJetCollectionTag)
  {
    edm::Handle < BXVector< l1t::Jet > > l1tJetCollectionHandle;
    iEvent.getByToken(*(this -> _l1tJetCollectionTag), l1tJetCollectionHandle);
    const auto l1tJetGenJetPairs = this -> _matchParticleWithL1Object<>(genJetCollectionHandle, l1tJetCollectionHandle, 1);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tJetGenJetTree), l1tJetGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tJetIterator = l1tJetCollectionHandle -> begin(0); l1tJetIterator != l1tJetCollectionHandle -> end(0); l1tJetIterator++ )
    {
      // Only barrel leading l1t jets
      if ((l1tJetIterator->pt() > maxPt) && (l1tJetIterator->eta() < 1.44) && (l1tJetIterator->eta() > -1.44))
      {
        save = true;        
        maxPt = l1tJetIterator -> pt();
        this -> _l1tObjectParticle.hwQual = l1tJetIterator -> hwQual();      
        this -> _l1tObjectParticle.id = (l1tJetIterator - l1tJetCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tJetIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tJetIterator -> eta();
        this -> _l1tObjectParticle.phi = l1tJetIterator -> phi();
      }
    }
    if (save) this -> _l1tLeadingJetTree -> Fill();
    
    for (auto l1tJetIterator = l1tJetCollectionHandle -> begin(0); l1tJetIterator != l1tJetCollectionHandle -> end(0); l1tJetIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tJetIterator - l1tJetCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tJetIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tJetIterator -> eta();
      this -> _l1tObjectParticle.phi = l1tJetIterator -> phi();
      this -> _l1tObjectParticle.hwQual = l1tJetIterator -> hwQual();            
      this -> _l1tJetTree -> Fill();
    }
    
  }

//  if (this->_l1tCaloTowerCollectionTag)
//  {
//    edm::Handle<BXVector<l1t::CaloTower>> l1tCaloTowerCollectionHandle;
//    iEvent.getByToken(*(this->_l1tCaloTowerCollectionTag), l1tCaloTowerCollectionHandle);
//    auto l1tCaloTowerGenJetPairs = this->_matchGenJetWithL1Object<>(genJetCollectionHandle, l1tCaloTowerCollectionHandle);
//    this->_fillTreeWithMatchedPairs(*(this->_l1tCaloTowerGenJetTree), l1tCaloTowerGenJetPairs);
//    float maxPt = 0;
//    bool save = false;
//    for (auto l1tCaloTowerIterator = l1tCaloTowerCollectionHandle->begin(0); l1tCaloTowerIterator != l1tCaloTowerCollectionHandle->end(0); l1tCaloTowerIterator++)
//    {
//      if ((l1tCaloTowerIterator->pt() > maxPt) && (l1tCaloTowerIterator->eta() < 1.44) && (l1tCaloTowerIterator->eta() > -1.44))
//      {
//        save = true;
//        maxPt = l1tCaloTowerIterator->pt();
//        this->_l1tObjectParticle.hwQual = l1tCaloTowerIterator->hwQual();
//        this->_l1tObjectParticle.id = (l1tCaloTowerIterator - l1tCaloTowerCollectionHandle->begin(0));
//        this->_l1tObjectParticle.pt = l1tCaloTowerIterator->pt();
//        this->_l1tObjectParticle.eta = l1tCaloTowerIterator->eta();
//        this->_l1tObjectParticle.phi = l1tCaloTowerIterator->phi();
//      }
//    }
//    if (save) this->_l1tLeadingCaloTowerTree->Fill();
//
//    for (auto l1tCaloTowerIterator = l1tCaloTowerCollectionHandle->begin(0); l1tCaloTowerIterator != l1tCaloTowerCollectionHandle->end(0); l1tCaloTowerIterator++)
//    {
//      this->_l1tObjectParticle.id = (l1tCaloTowerIterator - l1tCaloTowerCollectionHandle->begin(0));
//      this->_l1tObjectParticle.pt = l1tCaloTowerIterator->pt();
//      this->_l1tObjectParticle.eta = l1tCaloTowerIterator->eta();
//      this->_l1tObjectParticle.phi = l1tCaloTowerIterator->phi();
//      this->_l1tObjectParticle.hwQual = l1tCaloTowerIterator->hwQual();
//      this->_l1tCaloTowerTree->Fill();
//    }
//  }

  if (this -> _l1tEGammaCollectionTag)
  {
    edm::Handle < BXVector< l1t::EGamma > > l1tEGammaCollectionHandle;
    iEvent.getByToken(*(this -> _l1tEGammaCollectionTag), l1tEGammaCollectionHandle);
    const auto l1tEGammaGenJetPairs = this -> _matchParticleWithL1Object<>(genJetCollectionHandle, l1tEGammaCollectionHandle, 1);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tEGammaGenJetTree), l1tEGammaGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tEGammaIterator = l1tEGammaCollectionHandle -> begin(0); l1tEGammaIterator != l1tEGammaCollectionHandle -> end(0); l1tEGammaIterator++ )
    {
      if (l1tEGammaIterator -> pt() > maxPt)
      {
        save = true;
        maxPt = l1tEGammaIterator -> pt();
        this -> _l1tObjectParticle.hwQual = l1tEGammaIterator -> hwQual();        
        this -> _l1tObjectParticle.id = (l1tEGammaIterator - l1tEGammaCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tEGammaIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tEGammaIterator -> eta();
        this -> _l1tObjectParticle.phi = l1tEGammaIterator -> phi();
      }
    }
    if (save) this -> _l1tLeadingEGammaTree -> Fill();
    
    for (auto l1tEGammaIterator = l1tEGammaCollectionHandle -> begin(0); l1tEGammaIterator != l1tEGammaCollectionHandle -> end(0); l1tEGammaIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tEGammaIterator - l1tEGammaCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tEGammaIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tEGammaIterator -> eta();
      this -> _l1tObjectParticle.phi = l1tEGammaIterator -> phi();
      this -> _l1tObjectParticle.hwQual = l1tEGammaIterator -> hwQual();                  
      this -> _l1tEGammaTree -> Fill();
    }
    
  }
  
  if (this -> _l1tTauCollectionTag)
  {
    edm::Handle < BXVector< l1t::Tau > > l1tTauCollectionHandle;
    iEvent.getByToken(*(this -> _l1tTauCollectionTag), l1tTauCollectionHandle);
    const auto l1tTauGenJetPairs = this -> _matchParticleWithL1Object<>(genJetCollectionHandle, l1tTauCollectionHandle, 1);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tTauGenJetTree), l1tTauGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tTauIterator = l1tTauCollectionHandle -> begin(0); l1tTauIterator != l1tTauCollectionHandle -> end(0); l1tTauIterator++ )
    {
      if (l1tTauIterator -> pt() > maxPt)
      {
        save = true;
        maxPt = l1tTauIterator -> pt();
        this -> _l1tObjectParticle.id = (l1tTauIterator - l1tTauCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tTauIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tTauIterator -> eta();
        this -> _l1tObjectParticle.hwQual = l1tTauIterator -> hwQual();                        
        this -> _l1tObjectParticle.phi = l1tTauIterator -> phi();
      }
    }
    if (save) this -> _l1tLeadingTauTree -> Fill();
    
    for (auto l1tTauIterator = l1tTauCollectionHandle -> begin(0); l1tTauIterator != l1tTauCollectionHandle -> end(0); l1tTauIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tTauIterator - l1tTauCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tTauIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tTauIterator -> eta();
      this -> _l1tObjectParticle.hwQual = l1tTauIterator -> hwQual();                        
      this -> _l1tObjectParticle.phi = l1tTauIterator -> phi();
      this -> _l1tTauTree -> Fill();
    }
    
  }

  // Saving every genJet to get the misidentification probability

  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    this -> _genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
    this -> _genJetParticle.pt = genJetIterator -> pt();
    //std::cout << "Jet pt " << genJetIterator -> pt() << std::endl;
    this -> _genJetParticle.eta = genJetIterator -> eta();
    this -> _genJetParticle.phi = genJetIterator -> phi();
    this -> _genJetTree -> Fill();
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
    this -> _deltaR2 = deltaR2;
    this -> _matchingQuality = matchingQuality;
    
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
